// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential correctness: the new `ppvm-lossy-pauli-word-2::LossyPauliWord`
//! must agree with the old `ppvm-pauli-word::loss::LossyPauliWord` reference on
//! every *observable* operation, driven by the shared seeded generators of
//! `ppvm-conformance-2` (`seeded_rng` + `random_lossy_pauli_string`, which marks
//! ≈1/5 of the sites `Lost`).
//!
//! What is compared (design: `traits-2-implementation-plan.md` Phase 2;
//! `word-data-structures.md` §"Lossy Pauli word"):
//! construction from a shared lossy string, per-site `get(i)` (`LossySite`),
//! `weight()` (non-identity present sites **and** lost sites), `loss_weight()`,
//! `iter()`, `x_bit`/`z_bit`/`is_lost` reads, the Pauli product on the
//! present-site projection (result bits *and* the emitted `iᵏ` phase), and
//! Clifford conjugation for H/S/CNOT/CZ on present sites (resulting word).
//!
//! We compare **observable algebra, not raw hash digests** — the two crates
//! finalize the structural hash with different folds by design (the new one
//! splits X/Z from loss into two `OnceLock` components), so digest equality is
//! neither expected nor required. The hash *contract* is checked separately at
//! the bottom of this file.
//!
//! ## Clifford: bit-only, loss-preserving
//!
//! Both the old and new lossy words realize the *bit-only* (`Sp(2n,2)`) Clifford
//! map and additionally **no-op on any lost qubit** (the design's loss guard:
//! `word-data-structures.md` §"Loss-specific behavior"). The differential
//! therefore compares the resulting lossy word symbol-for-symbol; the
//! conjugation *sign* a phased word would pick up is dropped by both and is
//! validated against the ℤ[i] oracle on the present projection in
//! `lossy_pauli_word_lean.rs`.

use std::hash::BuildHasher;

use ppvm_conformance_2::{GateOp, random_circuit, random_lossy_pauli_string, seeded_rng};

// New crate under test.
use ppvm_lossy_pauli_word_2::LossyPauliWord as NewLossyWord;
use ppvm_pauli_word_2::PauliWord as NewWord;
use ppvm_traits_2::{
    Clifford as NewClifford, IdentityBuildHasher, Indexable, KeyProduct, LossySite,
    Pauli as NewPauli, PauliBits as NewPauliBits, Phase, Word as NewWordTrait,
};

// Old reference crate.
use ppvm_pauli_word::loss::LossyPauliWord as OldLossyWord;
use ppvm_pauli_word::phase::PhasedPauliWord;
use ppvm_traits::char::Pauli as OldPauli;
use ppvm_traits::traits::{Clifford as OldClifford, PauliIter, PauliWordTrait};

/// A 64-qubit-capacity new lossy word with the default fxhash finalizer.
type NewLossy = NewLossyWord<u64>;
/// A 64-qubit-capacity old lossy word.
type OldLossy = OldLossyWord<u64>;
/// The new *ordinary* word, used to multiply the present-site projection.
type New = NewWord<u64>;
/// The old *phased* word, used to multiply the present-site projection.
type OldPhased = PhasedPauliWord<u64>;

/// Seeds every property sweep replays under.
const SEEDS: [u64; 12] = [1, 2, 3, 7, 11, 42, 99, 123, 777, 2024, 31337, 88888];
/// Qubit widths exercised (≤ 64, the `u64` backing capacity).
const WIDTHS: [usize; 6] = [1, 2, 3, 5, 16, 60];

/// Render the new `LossySite<Pauli>` the way the old `Pauli` renders (`L` for
/// loss), so the two words can be diffed symbol-for-symbol.
fn new_site_char(s: LossySite<NewPauli>) -> char {
    match s {
        LossySite::Present(NewPauli::I) => 'I',
        LossySite::Present(NewPauli::X) => 'X',
        LossySite::Present(NewPauli::Y) => 'Y',
        LossySite::Present(NewPauli::Z) => 'Z',
        LossySite::Lost => 'L',
    }
}

// ---------------------------------------------------------------------------
// 1. Construction, inspection: get / weight / loss_weight / iter / bit+loss reads
// ---------------------------------------------------------------------------

#[test]
fn construction_and_inspection_match_old() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let s = random_lossy_pauli_string(&mut rng, n);
            let new: NewLossy = s.as_str().into();
            let old: OldLossy = OldLossy::from(s.as_str());

            // n_sites / n_qubits.
            assert_eq!(new.n_sites(), old.n_qubits(), "width {n} seed {seed}: {s}");

            // Whole-word rendering (per-site `get` composed).
            assert_eq!(new.to_string(), old.to_string(), "display {s}");

            // Per-site get(i) (LossySite vs old L-carrying Pauli), and
            // x_bit / z_bit / is_lost reads.
            for i in 0..n {
                assert_eq!(
                    new_site_char(new.get(i)).to_string(),
                    old.get(i).to_string(),
                    "get({i}) of {s}"
                );
                assert_eq!(new.x_bit(i), old.get_xbit(i), "x_bit({i}) of {s}");
                assert_eq!(new.z_bit(i), old.get_zbit(i), "z_bit({i}) of {s}");
                assert_eq!(
                    NewLossyWord::is_lost(&new, i),
                    old.get_lbit(i),
                    "is_lost({i}) of {s}"
                );
                // The `PauliBits::is_lost` override must agree with the inherent one.
                assert_eq!(
                    NewPauliBits::is_lost(&new, i),
                    NewLossyWord::is_lost(&new, i),
                    "PauliBits::is_lost({i}) of {s}"
                );
            }

            // weight() — non-identity present sites AND lost sites.
            assert_eq!(new.weight(), old.weight(), "weight of {s}");
            // loss_weight() — only the lost sites.
            assert_eq!(new.loss_weight(), old.loss_weight(), "loss_weight of {s}");

            // iter() sequence.
            let new_iter: Vec<String> = new.iter().map(|p| new_site_char(p).to_string()).collect();
            let old_iter: Vec<String> = PauliIter::iter(&old).map(|p| p.to_string()).collect();
            assert_eq!(new_iter, old_iter, "iter of {s}");
        }
    }
}

// ---------------------------------------------------------------------------
// 2. Pauli product on the present-site projection: result bits AND phase iᵏ
// ---------------------------------------------------------------------------

/// The present-site projection of a *new* lossy word as a plain Pauli string:
/// every `Lost` site becomes identity `I`, present sites keep their Pauli. This
/// is exactly the ordinary Pauli word the loss-agnostic product acts on.
fn new_present_projection(w: &NewLossy) -> String {
    (0..w.n_sites())
        .map(|i| match w.get(i) {
            LossySite::Lost => 'I',
            other => new_site_char(other),
        })
        .collect()
}

/// The present-site projection of an *old* lossy word (`L ↦ I`), read through the
/// old word's own `get`.
fn old_present_projection(w: &OldLossy) -> String {
    (0..w.n_qubits())
        .map(|i| match w.get(i) {
            OldPauli::L => 'I',
            p => p.to_string().chars().next().unwrap(),
        })
        .collect()
}

/// Old product of two present-site projections: the phased product of two
/// phaseless (`+`) words yields exactly the phase exponent `k` in its `phase`
/// field (both start at phase `0`), so it equals what the new `key_mul` emits.
fn old_product(v: &str, w: &str) -> (String, u8) {
    let a: OldPhased = format!("+{v}").as_str().into();
    let b: OldPhased = format!("+{w}").as_str().into();
    let prod = a * b;
    (prod.word.to_string(), prod.phase)
}

fn new_product(v: &str, w: &str) -> (String, u8) {
    let a: New = v.into();
    let b: New = w.into();
    let (word, phase): (New, Phase) = a.key_mul(&b);
    (word.to_string(), phase.exponent())
}

#[test]
fn present_projection_product_matches_old() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let sv = random_lossy_pauli_string(&mut rng, n);
            let sw = random_lossy_pauli_string(&mut rng, n);

            let nv: NewLossy = sv.as_str().into();
            let ov: OldLossy = OldLossy::from(sv.as_str());
            let nw: NewLossy = sw.as_str().into();
            let ow: OldLossy = OldLossy::from(sw.as_str());

            // The projections read through each crate's own lossy `get` must agree.
            let pv_new = new_present_projection(&nv);
            let pv_old = old_present_projection(&ov);
            let pw_new = new_present_projection(&nw);
            let pw_old = old_present_projection(&ow);
            assert_eq!(pv_new, pv_old, "present projection v of {sv}");
            assert_eq!(pw_new, pw_old, "present projection w of {sw}");

            // The product on the present sites (bits + emitted phase) must agree.
            assert_eq!(
                new_product(&pv_new, &pw_new),
                old_product(&pv_old, &pw_old),
                "present-projection product {sv} * {sw}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// 3. Clifford conjugation: H / S / CNOT / CZ (resulting lossy word, loss guard)
// ---------------------------------------------------------------------------

fn replay_clifford_new(word: &mut NewLossy, circuit: &[GateOp]) {
    for &op in circuit {
        match op {
            GateOp::H(q) => word.h(q),
            GateOp::S(q) => word.s(q),
            GateOp::Cnot(c, t) => word.cnot(c, t),
            GateOp::Rx(..) | GateOp::Rz(..) => {}
        }
    }
}

fn replay_clifford_old(word: &mut OldLossy, circuit: &[GateOp]) {
    for &op in circuit {
        match op {
            GateOp::H(q) => word.h(q),
            GateOp::S(q) => word.s(q),
            GateOp::Cnot(c, t) => word.cnot(c, t),
            GateOp::Rx(..) | GateOp::Rz(..) => {}
        }
    }
}

#[test]
fn clifford_replay_matches_old() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[2usize, 3, 5, 16] {
            let s = random_lossy_pauli_string(&mut rng, n);
            let circuit = random_circuit(&mut rng, n, 200);

            let mut new: NewLossy = s.as_str().into();
            let mut old: OldLossy = OldLossy::from(s.as_str());

            replay_clifford_new(&mut new, &circuit);
            replay_clifford_old(&mut old, &circuit);

            assert_eq!(new.to_string(), old.to_string(), "H/S/CNOT {s} seed {seed}");
            // Loss is untouched by the Clifford bit ops (same lost sites survive).
            assert_eq!(
                new.loss_weight(),
                old.loss_weight(),
                "loss preserved by Clifford {s} seed {seed}"
            );
        }
    }
}

/// CNOT and CZ are exercised exhaustively over the full two-qubit lossy alphabet
/// (`I/X/Y/Z/L`² = 25 words), including every present/lost combination, so the
/// loss guard on a mixed present-control/lost-target pair is covered directly.
#[test]
fn two_qubit_gates_match_old_exhaustive() {
    const A: [char; 5] = ['I', 'X', 'Y', 'Z', 'L'];
    for a in A {
        for b in A {
            let s: String = [a, b].iter().collect();

            let mut new_cnot: NewLossy = s.as_str().into();
            let mut old_cnot: OldLossy = OldLossy::from(s.as_str());
            new_cnot.cnot(0, 1);
            old_cnot.cnot(0, 1);
            assert_eq!(new_cnot.to_string(), old_cnot.to_string(), "CNOT {s}");

            let mut new_cz: NewLossy = s.as_str().into();
            let mut old_cz: OldLossy = OldLossy::from(s.as_str());
            new_cz.cz(0, 1);
            old_cz.cz(0, 1);
            assert_eq!(new_cz.to_string(), old_cz.to_string(), "CZ {s}");
        }
    }
}

/// The extended Clifford set (`S†`, `√X`, `√X†`, `√Y`, `√Y†`, `CY`) exhaustively
/// over the lossy alphabet.
///
/// The new blanket `CliffordExtensions` (`ppvm-traits-2/src/pauli.rs`) expresses
/// each gate as a *product* of audited generators, while the old blanket wrote
/// one fused bit rule per gate. That makes the loss guard the sharpest case:
/// `cy` decomposes into `s(t)`, `cnot(c,t)`, `s_dag(t)`, so with a **lost control
/// and a present target** the two `S`-family steps still run where the old
/// whole-gate skip ran nothing. They must cancel (`S` and `S†` share the `z ⊕= x`
/// bit map), leaving the word untouched — that is what the `L`-containing rows
/// below pin against the old reference.
#[test]
fn clifford_extension_gates_match_old_exhaustive() {
    use ppvm_traits::traits::CliffordExtensions as OldCliffordExtensions;
    use ppvm_traits_2::CliffordExtensions as NewCliffordExtensions;

    const A: [char; 5] = ['I', 'X', 'Y', 'Z', 'L'];

    for c in A {
        let s = c.to_string();
        for (name, new_gate, old_gate) in [
            (
                "S†",
                &(|w: &mut NewLossy| w.s_dag(0)) as &dyn Fn(&mut NewLossy),
                &(|w: &mut OldLossy| w.s_dag(0)) as &dyn Fn(&mut OldLossy),
            ),
            ("√X", &|w| w.sqrt_x(0), &|w| w.sqrt_x(0)),
            ("√X†", &|w| w.sqrt_x_dag(0), &|w| w.sqrt_x_dag(0)),
            ("√Y", &|w| w.sqrt_y(0), &|w| w.sqrt_y(0)),
            ("√Y†", &|w| w.sqrt_y_dag(0), &|w| w.sqrt_y_dag(0)),
        ] {
            let mut new: NewLossy = s.as_str().into();
            let mut old: OldLossy = OldLossy::from(s.as_str());
            new_gate(&mut new);
            old_gate(&mut old);
            assert_eq!(new.to_string(), old.to_string(), "{name} {s}");
        }
    }

    for a in A {
        for b in A {
            let s: String = [a, b].iter().collect();
            let mut new: NewLossy = s.as_str().into();
            let mut old: OldLossy = OldLossy::from(s.as_str());
            new.cy(0, 1);
            old.cy(0, 1);
            assert_eq!(new.to_string(), old.to_string(), "CY {s}");

            // The stim alias must be the same gate.
            let mut new_alias: NewLossy = s.as_str().into();
            new_alias.zcy(0, 1);
            assert_eq!(new_alias.to_string(), old.to_string(), "ZCY {s}");
        }
    }
}

/// H and S exhaustively over the single-qubit lossy alphabet.
#[test]
fn single_qubit_gates_match_old_exhaustive() {
    for c in ['I', 'X', 'Y', 'Z', 'L'] {
        let s = c.to_string();

        let mut new_h: NewLossy = s.as_str().into();
        let mut old_h: OldLossy = OldLossy::from(s.as_str());
        new_h.h(0);
        old_h.h(0);
        assert_eq!(new_h.to_string(), old_h.to_string(), "H {s}");

        let mut new_s: NewLossy = s.as_str().into();
        let mut old_s: OldLossy = OldLossy::from(s.as_str());
        new_s.s(0);
        old_s.s(0);
        assert_eq!(new_s.to_string(), old_s.to_string(), "S {s}");
    }
}

// ---------------------------------------------------------------------------
// 4. Hash CONTRACT (not raw-digest equality)
// ---------------------------------------------------------------------------

#[test]
fn hash_writes_exactly_key_hash() {
    // Design: `Hash for Self` is exactly `state.write_u64(self.key_hash())`, so
    // the identity build-hasher reproduces `key_hash()` bit for bit.
    let bh = IdentityBuildHasher;
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let s = random_lossy_pauli_string(&mut rng, n);
            let w: NewLossy = s.as_str().into();
            assert_eq!(bh.hash_one(&w), w.key_hash(), "hash==key_hash for {s}");
        }
    }
}

#[test]
fn structurally_equal_keys_hash_equal() {
    // Equal lossy keys *including the loss plane* => equal digest.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let s = random_lossy_pauli_string(&mut rng, n);
            let a: NewLossy = s.as_str().into();
            let b: NewLossy = s.as_str().into();
            assert_eq!(a, b, "structural eq {s}");
            assert_eq!(a.key_hash(), b.key_hash(), "equal keys => equal digest {s}");
        }
    }
}

#[test]
fn loss_plane_participates_in_digest() {
    // Two words with the SAME present X/Z pattern but a DIFFERENT loss plane must
    // (with overwhelming probability) hash differently: loss is part of the
    // structural identity, and `combine_components` is domain-separated, so a
    // loss-only difference is not allowed to collapse. Sweep random pairs and
    // require the digest to change whenever the loss mask changes.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[2usize, 3, 5, 16, 60] {
            // Build a present word, then a variant with one extra site lost.
            let base: NewLossy = random_lossy_pauli_string(&mut rng, n).as_str().into();
            // Find a currently-present site to flip to Lost (skip if none).
            let present = (0..n).find(|&i| !NewLossyWord::is_lost(&base, i));
            if let Some(i) = present {
                let mut variant = base.clone();
                variant.set(i, LossySite::Lost);
                assert_ne!(base, variant, "flipping site {i} to Lost changed the word");
                assert_ne!(
                    base.key_hash(),
                    variant.key_hash(),
                    "loss-plane difference must change the digest (seed {seed} n {n})"
                );
            }
        }
    }
}

#[test]
fn avalanche_low_collision_distribution() {
    // Distribution property (the design's stated hash contract, not a type-level
    // guarantee): over many distinct random 60-qubit *lossy* words, the low 12
    // bits hashbrown uses for a bucket must spread out, and full-digest
    // collisions must be rare.
    use std::collections::HashSet;

    let mut low_bits: HashSet<u16> = HashSet::new();
    let mut digests: HashSet<u64> = HashSet::new();
    let mut seen_words: HashSet<String> = HashSet::new();
    let mut rng = seeded_rng(0x105537);
    let mut distinct = 0usize;

    while distinct < 4096 {
        let s = random_lossy_pauli_string(&mut rng, 60);
        if !seen_words.insert(s.clone()) {
            continue; // count only structurally-distinct keys
        }
        distinct += 1;
        let w: NewLossy = s.as_str().into();
        let d = w.key_hash();
        low_bits.insert((d & 0xfff) as u16);
        digests.insert(d);
    }

    // Full-digest collisions on 4096 random 60-qubit lossy words: expect ~none.
    assert!(
        digests.len() >= distinct - 1,
        "too many full-digest collisions: {} distinct digests for {distinct} keys",
        digests.len(),
    );
    // Low 12 bits (4096 buckets): expected occupancy ≈63% (~2590). Require > 2000.
    assert!(
        low_bits.len() > 2000,
        "low 12 bits collapsed into {} buckets (of 4096)",
        low_bits.len(),
    );
}
