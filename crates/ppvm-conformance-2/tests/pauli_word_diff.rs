// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential correctness: the new `ppvm-pauli-word-2::PauliWord` must agree
//! with the old `ppvm-pauli-word` reference on every *observable* Pauli-algebra
//! operation, driven by the shared seeded generators of `ppvm-conformance-2`.
//!
//! What is compared (design: `traits-2-implementation-plan.md` Phase 2):
//! construction from a shared Pauli string, per-site `get(i)`, `weight()`,
//! `iter()`, `x_bit`/`z_bit`, the single- and n-qubit twisted product (result
//! bits *and* the emitted `iᵏ` phase), and Clifford conjugation for H/S/CNOT/CZ
//! (resulting word bits).
//!
//! We compare **observable algebra, not raw hash digests** — the two crates
//! share the `phase/mul.rs` kernel but finalize the structural hash with
//! different folds by design, so digest equality is neither expected nor
//! required (the hash *contract* is checked separately below).
//!
//! ## Clifford phase deltas: a design-accepted asymmetry
//!
//! The old bare `PauliWord` and the new bare `PauliWord` both realize the
//! *bit-only* (`Sp(2n,2)`) Clifford map and drop the conjugation sign — the new
//! crate's `PhaseTrack` is a documented no-op on a phaseless word
//! (`ppvm-pauli-word-2/src/clifford.rs`). So the differential here can only
//! compare the **resulting word**. The old *phased* word tracks a real phase
//! delta; we assert the new word's bits equal that phased word's bits (the `Sp`
//! part the bare word is responsible for), and the sign semantics are validated
//! against the ℤ[i] matrix oracle in `pauli_word_lean.rs`.

use std::hash::BuildHasher;

use ppvm_conformance_2::{GateOp, random_circuit, random_pauli_string, seeded_rng};

// New crate under test.
use ppvm_pauli_word_2::PauliWord as NewWord;
use ppvm_traits_2::{
    Clifford as NewClifford, IdentityBuildHasher, Indexable, KeyProduct, Pauli as NewPauli,
    PauliBits, Phase, Word as NewWordTrait,
};

/// The new `Pauli` carries no `Display`; render it the same way the old one does.
fn new_pauli_char(p: NewPauli) -> char {
    match p {
        NewPauli::I => 'I',
        NewPauli::X => 'X',
        NewPauli::Y => 'Y',
        NewPauli::Z => 'Z',
    }
}

// Old reference crate.
use ppvm_pauli_word::phase::PhasedPauliWord;
use ppvm_pauli_word::word::PauliWord as OldBareWord;
use ppvm_traits::traits::{Clifford as OldClifford, PauliIter, PauliWordTrait};

/// A 64-qubit-capacity new word with the default fxhash finalizer.
type New = NewWord<u64>;
/// A 64-qubit-capacity old bare word (phaseless bit map).
type OldBare = OldBareWord<u64>;
/// A 64-qubit-capacity old phased word (the `Sp` bits + a `ℤ/4` phase).
type OldPhased = PhasedPauliWord<u64>;

/// Seeds every property sweep replays under.
const SEEDS: [u64; 12] = [1, 2, 3, 7, 11, 42, 99, 123, 777, 2024, 31337, 88888];
/// Qubit widths exercised (≤ 64, the `u64` backing capacity).
const WIDTHS: [usize; 6] = [1, 2, 3, 5, 16, 60];

// ---------------------------------------------------------------------------
// 1. Construction, inspection: get / weight / iter / x_bit / z_bit
// ---------------------------------------------------------------------------

#[test]
fn construction_and_inspection_match_old() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let s = random_pauli_string(&mut rng, n);
            let new: New = s.as_str().into();
            let old: OldBare = OldBare::from(s.as_str());

            // n_sites / n_qubits.
            assert_eq!(new.n_sites(), old.n_qubits(), "width {n} seed {seed}: {s}");

            // Whole-word rendering (per-site `get` composed).
            assert_eq!(new.to_string(), old.to_string(), "display {s}");

            // Per-site get(i), and x_bit/z_bit reads.
            for i in 0..n {
                assert_eq!(
                    new_pauli_char(new.get(i)).to_string(),
                    old.get(i).to_string(),
                    "get({i}) of {s}"
                );
                assert_eq!(new.x_bit(i), old.get_xbit(i), "x_bit({i}) of {s}");
                assert_eq!(new.z_bit(i), old.get_zbit(i), "z_bit({i}) of {s}");
            }

            // weight().
            assert_eq!(new.weight(), old.weight(), "weight of {s}");

            // iter() sequence.
            let new_iter: Vec<String> = new.iter().map(|p| new_pauli_char(p).to_string()).collect();
            let old_iter: Vec<String> = PauliIter::iter(&old).map(|p| p.to_string()).collect();
            assert_eq!(new_iter, old_iter, "iter of {s}");
        }
    }
}

// ---------------------------------------------------------------------------
// 2. Twisted product: result bits AND emitted phase iᵏ
// ---------------------------------------------------------------------------

/// The old phased product of two phaseless (`+`) words yields, in its result
/// `phase` field, exactly the phase exponent `k` that the new `key_mul` emits —
/// because both start at phase `0`, so `result.phase == phaseExp(v, w)`.
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
fn single_qubit_product_matches_old_exhaustive() {
    for &v in &["I", "X", "Y", "Z"] {
        for &w in &["I", "X", "Y", "Z"] {
            assert_eq!(new_product(v, w), old_product(v, w), "{v} * {w}");
        }
    }
}

#[test]
fn n_qubit_product_matches_old() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let v = random_pauli_string(&mut rng, n);
            let w = random_pauli_string(&mut rng, n);
            let (nw, np) = new_product(&v, &w);
            let (ow, op) = old_product(&v, &w);
            assert_eq!(nw, ow, "product bits {v} * {w}");
            assert_eq!(np, op, "product phase exponent {v} * {w}");
        }
    }
}

// ---------------------------------------------------------------------------
// 3. Clifford conjugation: H / S / CNOT / CZ (resulting word bits)
// ---------------------------------------------------------------------------

/// Replay only the Clifford generators (`H`/`S`/`CNOT`) of a circuit; the
/// rotations `Rx`/`Rz` are not Clifford and are skipped by both backends.
fn replay_clifford_new(word: &mut New, circuit: &[GateOp]) {
    for &op in circuit {
        match op {
            GateOp::H(q) => word.h(q),
            GateOp::S(q) => word.s(q),
            GateOp::Cnot(c, t) => word.cnot(c, t),
            GateOp::Rx(..) | GateOp::Rz(..) => {}
        }
    }
}

fn replay_clifford_old_bare(word: &mut OldBare, circuit: &[GateOp]) {
    for &op in circuit {
        match op {
            GateOp::H(q) => word.h(q),
            GateOp::S(q) => word.s(q),
            GateOp::Cnot(c, t) => word.cnot(c, t),
            GateOp::Rx(..) | GateOp::Rz(..) => {}
        }
    }
}

fn replay_clifford_old_phased(word: &mut OldPhased, circuit: &[GateOp]) {
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
            let s = random_pauli_string(&mut rng, n);
            let circuit = random_circuit(&mut rng, n, 200);

            let mut new: New = s.as_str().into();
            let mut old_bare: OldBare = OldBare::from(s.as_str());
            let mut old_phased: OldPhased = format!("+{s}").as_str().into();

            replay_clifford_new(&mut new, &circuit);
            replay_clifford_old_bare(&mut old_bare, &circuit);
            replay_clifford_old_phased(&mut old_phased, &circuit);

            // New bare word == old bare word (both are the bit-only Sp map).
            assert_eq!(
                new.to_string(),
                old_bare.to_string(),
                "bare {s} seed {seed}"
            );
            // New bare word's bits == the phased word's bits (the Sp part the
            // bare word owns; the phase delta the phased word additionally
            // tracks is dropped by design — see module docs).
            assert_eq!(
                new.to_string(),
                old_phased.word.to_string(),
                "phased-bits {s} seed {seed}"
            );
        }
    }
}

/// The extended Clifford set (`S†`, `√X`, `√X†`, `√Y`, `√Y†`, `CY`) — the new
/// blanket `CliffordExtensions` of `ppvm-traits-2/src/pauli.rs` against the old
/// blanket `impl<T: PauliWordTrait> CliffordExtensions for T`.
///
/// The new blanket derives each gate from audited generators (`√X ≃ H·S·H`, …)
/// rather than from a hand-written bit rule per gate, so the `Sp(2n,2)` map it
/// produces is exactly what must be diffed here; the accompanying signs (dropped
/// by a phaseless word) are diffed in `phased_pauli_word_diff.rs`.
#[test]
fn clifford_extension_bit_maps_match_old_exhaustive() {
    use ppvm_traits::traits::CliffordExtensions as OldCliffordExtensions;
    use ppvm_traits_2::CliffordExtensions as NewCliffordExtensions;

    for a in ["I", "X", "Y", "Z"] {
        let mut new_s: New = a.into();
        let mut old_s: OldBare = OldBare::from(a);
        new_s.s_dag(0);
        old_s.s_dag(0);
        assert_eq!(new_s.to_string(), old_s.to_string(), "S† {a}");

        let mut new_sx: New = a.into();
        let mut old_sx: OldBare = OldBare::from(a);
        new_sx.sqrt_x(0);
        old_sx.sqrt_x(0);
        assert_eq!(new_sx.to_string(), old_sx.to_string(), "√X {a}");

        let mut new_sxd: New = a.into();
        let mut old_sxd: OldBare = OldBare::from(a);
        new_sxd.sqrt_x_dag(0);
        old_sxd.sqrt_x_dag(0);
        assert_eq!(new_sxd.to_string(), old_sxd.to_string(), "√X† {a}");

        let mut new_sy: New = a.into();
        let mut old_sy: OldBare = OldBare::from(a);
        new_sy.sqrt_y(0);
        old_sy.sqrt_y(0);
        assert_eq!(new_sy.to_string(), old_sy.to_string(), "√Y {a}");

        let mut new_syd: New = a.into();
        let mut old_syd: OldBare = OldBare::from(a);
        new_syd.sqrt_y_dag(0);
        old_syd.sqrt_y_dag(0);
        assert_eq!(new_syd.to_string(), old_syd.to_string(), "√Y† {a}");

        for b in ["I", "X", "Y", "Z"] {
            let s = format!("{a}{b}");
            let mut new: New = s.as_str().into();
            let mut old: OldBare = OldBare::from(s.as_str());
            new.cy(0, 1);
            old.cy(0, 1);
            assert_eq!(new.to_string(), old.to_string(), "CY {s}");
        }
    }
}

/// CZ is not emitted by `random_circuit`, so cover the full two-qubit table
/// directly, comparing the new bit map against the old bare word's.
#[test]
fn cz_bit_map_matches_old_exhaustive() {
    for a in ["I", "X", "Y", "Z"] {
        for b in ["I", "X", "Y", "Z"] {
            let s = format!("{a}{b}");
            let mut new: New = s.as_str().into();
            let mut old: OldBare = OldBare::from(s.as_str());
            new.cz(0, 1);
            old.cz(0, 1);
            assert_eq!(new.to_string(), old.to_string(), "CZ {s}");
        }
    }
}

// ---------------------------------------------------------------------------
// 4. Hash CONTRACT (not raw-digest equality)
// ---------------------------------------------------------------------------

#[test]
fn hash_writes_exactly_key_hash() {
    // Design: `Hash for Self` is exactly `state.write_u64(self.key_hash())`,
    // so the identity build-hasher reproduces `key_hash()` bit for bit.
    let bh = IdentityBuildHasher;
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let s = random_pauli_string(&mut rng, n);
            let w: New = s.as_str().into();
            assert_eq!(bh.hash_one(&w), w.key_hash(), "hash==key_hash for {s}");
        }
    }
}

#[test]
fn structurally_equal_keys_hash_equal() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let s = random_pauli_string(&mut rng, n);
            let a: New = s.as_str().into();
            let b: New = s.as_str().into();
            assert_eq!(a, b, "structural eq {s}");
            assert_eq!(a.key_hash(), b.key_hash(), "equal keys => equal digest {s}");
        }
    }
}

#[test]
fn avalanche_low_collision_distribution() {
    // A distribution property test (design's stated hash contract, not a
    // type-level guarantee): over many distinct random 60-qubit words, the low
    // 12 bits hashbrown would use for a bucket must spread out, and full-digest
    // collisions must be rare.
    use std::collections::HashSet;

    let mut low_bits: HashSet<u16> = HashSet::new();
    let mut digests: HashSet<u64> = HashSet::new();
    let mut seen_words: HashSet<String> = HashSet::new();
    let mut rng = seeded_rng(0xC0FFEE);
    let mut distinct = 0usize;

    while distinct < 4096 {
        let s = random_pauli_string(&mut rng, 60);
        if !seen_words.insert(s.clone()) {
            continue; // ensure we count only structurally-distinct keys
        }
        distinct += 1;
        let w: New = s.as_str().into();
        let d = w.key_hash();
        low_bits.insert((d & 0xfff) as u16);
        digests.insert(d);
    }

    // Full-digest collisions on 4096 random 60-qubit words: expect ~none.
    assert!(
        digests.len() >= distinct - 1,
        "too many full-digest collisions: {} distinct digests for {distinct} keys",
        digests.len(),
    );
    // Low 12 bits (4096 buckets) should be well spread: with 4096 keys into
    // 4096 buckets the expected occupancy is ~63% (≈2590). Require > 2000.
    assert!(
        low_bits.len() > 2000,
        "low 12 bits collapsed into {} buckets (of 4096)",
        low_bits.len(),
    );
}
