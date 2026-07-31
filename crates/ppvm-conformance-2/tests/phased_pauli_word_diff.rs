// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential correctness: the new `ppvm-phased-pauli-word-2::PhasedPauliWord`
//! (`Phased<PauliWord>`) must agree with the old `ppvm-pauli-word::PhasedPauliWord`
//! reference on every *observable* phased-Pauli operation, driven by the shared
//! seeded generators of `ppvm-conformance-2`.
//!
//! Unlike the bare `PauliWord` diff (`pauli_word_diff.rs`), which can only compare
//! the `Sp(2n,2)` *bits* because the bare word's `PhaseTrack` is a no-op, the
//! phased word carries an explicit `ℤ₄` phase. So here we diff the **phase too**:
//!
//! * construction from a shared signed string (`"+iXYZ"` …): phase + bits;
//! * `Word` delegation (`n_sites`/`get`/`weight`/`iter`);
//! * the phased product — result X/Z bits **and** the accumulated `ℤ₄` phase;
//! * Clifford conjugation for `H`/`S`/`CNOT`/`CZ` — bits **and** the phase delta
//!   the bare word dropped (a seeded `random_circuit`'s Clifford gates replayed on
//!   matched words, plus the full two-qubit `CZ` table `random_circuit` omits);
//! * the phase accessor (`phase()` / `is_positive()`).
//!
//! We compare **observable algebra, not raw hash digests**: `Phased<W>` is
//! deliberately non-`Indexable`, so there is no digest to diff.

use ppvm_conformance_2::{GateOp, random_circuit, random_pauli_string, seeded_rng};

// New crate under test.
use ppvm_phased_pauli_word_2::PhasedPauliWord as NewPhased;
use ppvm_traits_2::{Clifford as NewClifford, Word as NewWordTrait};

// Old reference crate.
use ppvm_pauli_word::phase::PhasedPauliWord as OldPhasedTy;
use ppvm_traits::traits::{Clifford as OldClifford, PauliIter, PauliWordTrait};

/// A 64-qubit-capacity old phased word (the `Sp` bits + a `ℤ/4` phase `u8`).
type OldPhased = OldPhasedTy<u64>;

/// Seeds every property sweep replays under.
const SEEDS: [u64; 12] = [1, 2, 3, 7, 11, 42, 99, 123, 777, 2024, 31337, 88888];
/// Qubit widths exercised (≤ 64, the `u64` backing capacity).
const WIDTHS: [usize; 6] = [1, 2, 3, 5, 16, 60];
/// The four `ℤ₄` sign prefixes, in exponent order (`+ → 0`, `+i → 1`, …).
const PREFIXES: [&str; 4] = ["+", "+i", "-", "-i"];

/// Build a matched (new, old) pair from the *same* signed Pauli string.
fn pair(signed: &str) -> (NewPhased, OldPhased) {
    (NewPhased::from(signed), OldPhased::from(signed))
}

/// The new phase as the `u8` exponent the old crate stores (`+ → 0`, `+i → 1`,
/// `- → 2`, `-i → 3`), so the two phase encodings can be diffed directly.
fn new_phase_u8(w: &NewPhased) -> u8 {
    w.phase().exponent()
}

// ---------------------------------------------------------------------------
// 1. Construction + phase accessor + Word delegation
// ---------------------------------------------------------------------------

#[test]
fn construction_phase_and_delegation_match_old() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            for prefix in PREFIXES {
                let body = random_pauli_string(&mut rng, n);
                let signed = format!("{prefix}{body}");
                let (new, old) = pair(&signed);

                // Construction: phase exponent + whole-word rendering agree.
                assert_eq!(new_phase_u8(&new), old.phase, "phase of {signed}");
                assert_eq!(new.to_string(), old.to_string(), "display {signed}");

                // Phase accessor.
                assert_eq!(new.is_positive(), old.is_positive(), "is_positive {signed}");

                // Word delegation: n_sites / get / weight / iter, all against the
                // old word (the phase must not perturb site inspection).
                assert_eq!(new.n_sites(), old.n_qubits(), "n_sites {signed}");
                for i in 0..n {
                    assert_eq!(
                        format!("{:?}", new.get(i)).chars().next().unwrap(),
                        old.get(i).to_string().chars().next().unwrap(),
                        "get({i}) of {signed}"
                    );
                }
                assert_eq!(new.weight(), old.word.weight(), "weight {signed}");
                let new_iter: Vec<char> = new
                    .iter()
                    .map(|p| format!("{p:?}").chars().next().unwrap())
                    .collect();
                let old_iter: Vec<char> = PauliIter::iter(&old.word)
                    .map(|p| p.to_string().chars().next().unwrap())
                    .collect();
                assert_eq!(new_iter, old_iter, "iter of {signed}");
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 2. Phased product: result bits AND accumulated ℤ₄ phase
// ---------------------------------------------------------------------------

/// `(word_string, phase_u8)` of the new phased product `a · b`.
fn new_product(a: &NewPhased, b: &NewPhased) -> (String, u8) {
    let prod = a * b;
    (prod.word().to_string(), new_phase_u8(&prod))
}

/// `(word_string, phase_u8)` of the old phased product `a · b`.
fn old_product(a: OldPhased, b: OldPhased) -> (String, u8) {
    let prod = a * b;
    (prod.word.to_string(), prod.phase)
}

#[test]
fn single_qubit_product_matches_old_exhaustive() {
    // Every (phase, Pauli) × (phase, Pauli) single-qubit pair: 256 cases.
    for pl in PREFIXES {
        for l in ["I", "X", "Y", "Z"] {
            for pr in PREFIXES {
                for r in ["I", "X", "Y", "Z"] {
                    let (na, oa) = pair(&format!("{pl}{l}"));
                    let (nb, ob) = pair(&format!("{pr}{r}"));
                    assert_eq!(
                        new_product(&na, &nb),
                        old_product(oa, ob),
                        "{pl}{l} * {pr}{r}"
                    );
                }
            }
        }
    }
}

#[test]
fn n_qubit_product_matches_old_bits_and_phase() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            for pl in PREFIXES {
                for pr in PREFIXES {
                    let v = format!("{pl}{}", random_pauli_string(&mut rng, n));
                    let w = format!("{pr}{}", random_pauli_string(&mut rng, n));
                    let (na, oa) = pair(&v);
                    let (nb, ob) = pair(&w);
                    assert_eq!(
                        new_product(&na, &nb),
                        old_product(oa, ob),
                        "seed {seed}: {v} * {w}"
                    );
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 3. Clifford conjugation: H / S / CNOT / CZ — bits AND phase delta
// ---------------------------------------------------------------------------

fn replay_clifford_new(word: &mut NewPhased, circuit: &[GateOp]) {
    for &op in circuit {
        match op {
            GateOp::H(q) => word.h(q),
            GateOp::S(q) => word.s(q),
            GateOp::Cnot(c, t) => word.cnot(c, t),
            GateOp::Rx(..) | GateOp::Rz(..) => {}
        }
    }
}

fn replay_clifford_old(word: &mut OldPhased, circuit: &[GateOp]) {
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
fn clifford_replay_matches_old_bits_and_phase() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[2usize, 3, 5, 16] {
            // A random starting phase exercises phase carry-through under signed
            // conjugation (e.g. H(−iY) = +iY), not just deltas onto `+`.
            let prefix = PREFIXES[rng_prefix(&mut rng)];
            let signed = format!("{prefix}{}", random_pauli_string(&mut rng, n));
            let circuit = random_circuit(&mut rng, n, 200);

            let (mut new, mut old) = pair(&signed);
            replay_clifford_new(&mut new, &circuit);
            replay_clifford_old(&mut old, &circuit);

            // Bits AND phase must agree after the full Clifford replay.
            assert_eq!(
                new.word().to_string(),
                old.word.to_string(),
                "bits {signed} seed {seed}"
            );
            assert_eq!(new_phase_u8(&new), old.phase, "phase {signed} seed {seed}");
            // Whole-word display folds both together.
            assert_eq!(
                new.to_string(),
                old.to_string(),
                "display {signed} seed {seed}"
            );
        }
    }
}

/// A uniform index into [`PREFIXES`].
fn rng_prefix(rng: &mut rand::rngs::StdRng) -> usize {
    use rand::RngExt;
    rng.random_range(0..4usize)
}

/// `CZ` is not emitted by `random_circuit`, so cover the full two-qubit table
/// (every phase × every two-qubit Pauli), comparing new bits + phase vs old.
#[test]
fn cz_matches_old_exhaustive_bits_and_phase() {
    for prefix in PREFIXES {
        for a in ["I", "X", "Y", "Z"] {
            for b in ["I", "X", "Y", "Z"] {
                let signed = format!("{prefix}{a}{b}");
                let (mut new, mut old) = pair(&signed);
                new.cz(0, 1);
                old.cz(0, 1);
                assert_eq!(
                    new.word().to_string(),
                    old.word.to_string(),
                    "CZ bits {signed}"
                );
                assert_eq!(new_phase_u8(&new), old.phase, "CZ phase {signed}");
            }
        }
    }
}

/// A named single-qubit generator paired for the new and old phased words.
type GatePair = (&'static str, fn(&mut NewPhased), fn(&mut OldPhased));

/// The single-qubit generators (including the pure-sign `X`/`Y`/`Z`) exhaustively
/// vs old, over every starting phase — pins each per-gate phase delta.
#[test]
fn single_qubit_gates_match_old_exhaustive() {
    let gates: [GatePair; 5] = [
        ("H", |w| w.h(0), |w| w.h(0)),
        ("S", |w| w.s(0), |w| w.s(0)),
        ("X", |w| w.x(0), |w| w.x(0)),
        ("Y", |w| w.y(0), |w| w.y(0)),
        ("Z", |w| w.z(0), |w| w.z(0)),
    ];
    for (name, gnew, gold) in gates {
        for prefix in PREFIXES {
            for p in ["I", "X", "Y", "Z"] {
                let signed = format!("{prefix}{p}");
                let (mut new, mut old) = pair(&signed);
                gnew(&mut new);
                gold(&mut old);
                assert_eq!(
                    (new.word().to_string(), new_phase_u8(&new)),
                    (old.word.to_string(), old.phase),
                    "{name} {signed}"
                );
            }
        }
    }
}
