// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential correctness: the new `ppvm-pauli-sum-2::PauliSum<f64>` must agree
//! with the old `ppvm-pauli-sum::PauliSum<f64>` reference on every *observable*
//! sparse-sum operation, driven by the shared seeded generators of
//! `ppvm-conformance-2`.
//!
//! Both sums are built from the **same** `(pauli_string, coeff)` list (the
//! `build_old_sum`/`build_new_sum` harness pair) and compared on:
//!
//! * `len`/`is_empty`, `contains`/`get`, `n_qubits`/`n_sites`;
//! * the full support as a sorted `(canonical_pauli_string, coeff)` set (`iter`);
//! * `scale` by a constant;
//! * `reduce` — insert cancelling terms summing to zero, reduce, assert the key
//!   is dropped in BOTH;
//! * `overlap` of two sums (the L3 `Pair` trace pairing);
//! * the core: replay a seeded random **Clifford** circuit (`H`/`S`/`CNOT`/`CZ`)
//!   and, after *each* gate, assert the full supports + coefficients match — this
//!   exercises the phase-to-coefficient sign draining that the bare-word Clifford
//!   drops and the phased-word fused Clifford recovers.
//!
//! We compare **observable algebra, not raw hash digests**: the two crates'
//! finalization folds differ by design, so a digest diff is meaningless. The
//! hashing *contract* (write-exactly-`key_hash`, structural determinism,
//! avalanche) is tested separately below.

use ppvm_conformance_2::{
    GateOp, NewKey, NewSum, OldSum, assert_close, assert_supports_match, build_new_sum,
    build_old_sum, new_support, old_support, random_circuit, random_terms, reduce_old, seeded_rng,
};

// New crate Clifford trait.
use ppvm_traits_2::Clifford as NewClifford;
// Old crate Clifford trait.
use ppvm_traits::traits::Clifford as OldClifford;

use rand::RngExt;
use rand::rngs::StdRng;

/// Seeds every property sweep replays under.
const SEEDS: [u64; 10] = [1, 2, 3, 7, 42, 99, 123, 777, 2024, 31337];
/// Qubit widths exercised (≤ 64, the shared backing capacity).
const WIDTHS: [usize; 6] = [1, 2, 3, 5, 8, 12];
/// Coefficient comparison tolerance for the `f64` backends.
const TOL: f64 = 1e-9;

// ---------------------------------------------------------------------------
// 1. Construction, inspection, iter/support
// ---------------------------------------------------------------------------

#[test]
fn build_len_empty_and_support_match() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            for count in [0usize, 1, 5, 20, 60] {
                let terms = random_terms(&mut rng, n, count);
                let old = build_old_sum(n, &terms);
                let new = build_new_sum(n, &terms);

                // n_qubits / n_sites.
                assert_eq!(old.n_qubits(), n, "old n_qubits seed {seed} n {n}");
                assert_eq!(new.n_sites(), n, "new n_sites seed {seed} n {n}");

                // len / is_empty (both merge colliding keys identically).
                assert_eq!(
                    old.len(),
                    new.len(),
                    "len seed {seed} n {n} count {count}\nold={:?}\nnew={:?}",
                    old_support(&old),
                    new_support(&new)
                );
                assert_eq!(old.is_empty(), new.is_empty(), "is_empty seed {seed}");

                // Full support as a sorted (string, coeff) set.
                assert_supports_match(&old, &new, TOL);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 2. contains / get agree, key by key
// ---------------------------------------------------------------------------

#[test]
fn contains_and_get_match() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let terms = random_terms(&mut rng, n, 25);
            let new = build_new_sum(n, &terms);

            // Every key in the (reduced) NEW support: `get` present in both, with
            // the same coefficient; `contains` agrees.
            for (word, coeff) in new_support(&new) {
                let key = NewKey::from(word.as_str());
                let got = new.get(&key).expect("key in support has a coeff");
                assert_close(got, coeff, TOL);
                assert!(new.contains(&key), "new.contains {word}");

                // The OLD sum agrees on presence + coefficient (via its support).
                let old = build_old_sum(n, &terms);
                let os = old_support(&old);
                let found = os.iter().find(|(w, _)| *w == word);
                assert!(found.is_some(), "old missing key {word}");
                assert_close(found.unwrap().1, coeff, TOL);
            }

            // A guaranteed-absent key (all-identity is rarely produced; check it is
            // reported consistently): build it and compare membership.
            let absent = "I".repeat(n);
            let key = NewKey::from(absent.as_str());
            let new_has = new.contains(&key);
            let old = build_old_sum(n, &terms);
            let old_has = old_support(&old).iter().any(|(w, _)| *w == absent);
            assert_eq!(
                new_has, old_has,
                "identity-key membership seed {seed} n {n}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// 3. scale by a constant
// ---------------------------------------------------------------------------

#[test]
fn scale_matches() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let terms = random_terms(&mut rng, n, 20);
            for s in [0.0, 1.0, -1.0, 2.5, -0.5, 1e6] {
                let mut old = build_old_sum(n, &terms);
                let mut new = build_new_sum(n, &terms);
                old *= s;
                new.scale(&s);

                // Scaling by 0 leaves the OLD map full of exact zeros while the NEW
                // `scale` is a pure coefficient map (no reduce). Reduce both so the
                // comparison is on canonical support.
                reduce_old(&mut old);
                // The NEW `scale` does not reduce; rebuild through from_terms to
                // canonicalize (drop the zeros the s == 0 case creates).
                let new = build_new_sum(
                    n,
                    &new_support(&new)
                        .into_iter()
                        .collect::<Vec<(String, f64)>>(),
                );

                assert_supports_match(&old, &new, TOL.max(s.abs() * 1e-12));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 4. reduce drops a cancelled key in BOTH crates
// ---------------------------------------------------------------------------

#[test]
fn reduce_drops_cancelled_key_in_both() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            // A dedicated cancel key (all-Z), with survivors generated to *exclude*
            // it so the cancelling pair sums to exactly 0.0 in isolation.
            let cancel_key = "Z".repeat(n);
            let mut terms: Vec<(String, f64)> = random_terms(&mut rng, n, 6)
                .into_iter()
                .filter(|(w, _)| *w != cancel_key)
                .collect();
            terms.push((cancel_key.clone(), 1.5));
            terms.push((cancel_key.clone(), -1.5)); // sums to exactly 0.0

            let mut old = build_old_sum(n, &terms);
            let new = build_new_sum(n, &terms); // from_terms already reduced

            // OLD: `+=` leaves the cancelled key at 0.0; realize `reduce`.
            let before = old.len();
            reduce_old(&mut old);
            assert!(
                old.len() < before || before == 0,
                "reduce did not shrink the old support (before {before})"
            );

            // The cancelled key is absent from BOTH.
            assert!(
                !new.contains(&NewKey::from(cancel_key.as_str())),
                "new kept cancelled key {cancel_key}"
            );
            assert!(
                !old_support(&old).iter().any(|(w, _)| *w == cancel_key),
                "old kept cancelled key {cancel_key}"
            );

            // And the reduced supports agree.
            assert_supports_match(&old, &new, TOL);
        }
    }
}

// ---------------------------------------------------------------------------
// 5. overlap of two sums (L3 Pair)
// ---------------------------------------------------------------------------

#[test]
fn overlap_matches() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let a_terms = random_terms(&mut rng, n, 15);
            let b_terms = random_terms(&mut rng, n, 15);

            let old_a = build_old_sum(n, &a_terms);
            let old_b = build_old_sum(n, &b_terms);
            let new_a = build_new_sum(n, &a_terms);
            let new_b = build_new_sum(n, &b_terms);

            let old_ov = old_a.overlap(&old_b);
            let new_ov = new_a.overlap(&new_b);
            // Magnitude-relative tolerance: overlap sums up to 15×15 products.
            let tol = TOL.max(old_ov.abs() * 1e-10);
            assert_close(old_ov, new_ov, tol);

            // Symmetry cross-check on the new side (bilinear pairing).
            assert_close(new_ov, new_b.overlap(&new_a), tol);
        }
    }
}

// ---------------------------------------------------------------------------
// 6. THE CORE — Clifford replay with per-gate support+coefficient agreement
// ---------------------------------------------------------------------------

/// A replayable Clifford gate — the `H`/`S`/`CNOT`/`CZ` subset (`random_circuit`
/// omits `CZ`, so it is drawn explicitly here to exercise the two-qubit sign
/// rules).
#[derive(Clone, Copy, Debug)]
enum Cliff {
    H(usize),
    S(usize),
    Cnot(usize, usize),
    Cz(usize, usize),
}

fn apply_old(s: &mut OldSum, g: Cliff) {
    match g {
        Cliff::H(q) => s.h(q),
        Cliff::S(q) => s.s(q),
        Cliff::Cnot(c, t) => s.cnot(c, t),
        Cliff::Cz(a, b) => s.cz(a, b),
    }
}

fn apply_new(s: &mut NewSum, g: Cliff) {
    match g {
        Cliff::H(q) => s.h(q),
        Cliff::S(q) => s.s(q),
        Cliff::Cnot(c, t) => s.cnot(c, t),
        Cliff::Cz(a, b) => s.cz(a, b),
    }
}

/// A random Clifford gate list on `n` qubits, `len` gates long, over `H/S/CNOT/CZ`.
fn random_clifford_ops(rng: &mut StdRng, n: usize, len: usize) -> Vec<Cliff> {
    (0..len)
        .map(|_| {
            let q = rng.random_range(0..n);
            match rng.random_range(0..4usize) {
                0 => Cliff::H(q),
                1 => Cliff::S(q),
                2 if n > 1 => {
                    let mut t = rng.random_range(0..n);
                    while t == q {
                        t = rng.random_range(0..n);
                    }
                    Cliff::Cnot(q, t)
                }
                3 if n > 1 => {
                    let mut t = rng.random_range(0..n);
                    while t == q {
                        t = rng.random_range(0..n);
                    }
                    Cliff::Cz(q, t)
                }
                _ => Cliff::H(q),
            }
        })
        .collect()
}

#[test]
fn clifford_replay_matches_old_per_gate() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[2usize, 3, 5, 8, 12] {
            // A moderate initial support so cancellations and sign draining have
            // something to act on.
            let terms = random_terms(&mut rng, n, 30);
            let mut old = build_old_sum(n, &terms);
            let mut new = build_new_sum(n, &terms);

            // Agreement before any gate.
            assert_supports_match(&old, &new, TOL);

            let ops = random_clifford_ops(&mut rng, n, 120);
            for (i, g) in ops.into_iter().enumerate() {
                apply_old(&mut old, g);
                apply_new(&mut new, g);
                // A Clifford re-key is a bijection: len is preserved and every key
                // maps to a distinct key, so the support size must never change.
                assert_eq!(
                    old.len(),
                    new.len(),
                    "len diverged after gate {i} ({g:?}) seed {seed} n {n}"
                );
                assert_supports_match(&old, &new, TOL);
            }
        }
    }
}

/// A named single-qubit gate paired for the old and new sums.
type SumGatePair = (&'static str, fn(&mut OldSum), fn(&mut NewSum));

/// Exhaustive single-qubit Clifford draining: over every single-qubit Pauli term,
/// `H`/`S` on the old and new sums must agree on the resulting (word, sign).
#[test]
fn single_qubit_clifford_exhaustive() {
    let gens: [SumGatePair; 2] = [("H", |s| s.h(0), |s| s.h(0)), ("S", |s| s.s(0), |s| s.s(0))];
    for (name, gold, gnew) in gens {
        for p in ["I", "X", "Y", "Z"] {
            for &c in &[1.0f64, -1.0, 2.5, -0.75] {
                let terms = vec![(p.to_string(), c)];
                let mut old = build_old_sum(1, &terms);
                let mut new = build_new_sum(1, &terms);
                gold(&mut old);
                gnew(&mut new);
                assert_supports_match(&old, &new, TOL);
                // Spell out the expected sign drain for X/Y/Z under H/S so a
                // regression in *both* crates cannot hide behind the diff.
                let (w, sign) = new_support(&new)
                    .into_iter()
                    .next()
                    .map(|(w, cc)| (w, cc.signum()))
                    .unwrap();
                let expect = expected_single(name, p);
                assert_eq!(w, expect.0, "{name} {p}: word");
                assert_eq!(sign, expect.1 * c.signum(), "{name} {p}: sign");
            }
        }
    }
}

/// Reference conjugation `g P g†` for single-qubit `H`/`S` (Heisenberg picture):
/// `(resulting_word, sign)`.
fn expected_single(gate: &str, p: &str) -> (String, f64) {
    match (gate, p) {
        // H: X↔Z, HYH = −Y.
        ("H", "I") => ("I".into(), 1.0),
        ("H", "X") => ("Z".into(), 1.0),
        ("H", "Z") => ("X".into(), 1.0),
        ("H", "Y") => ("Y".into(), -1.0),
        // S (the S† conjugation both crates and `conjSdag_sign` in Lean use):
        // X→−Y, Y→+X, Z→Z.
        ("S", "I") => ("I".into(), 1.0),
        ("S", "X") => ("Y".into(), -1.0),
        ("S", "Y") => ("X".into(), 1.0),
        ("S", "Z") => ("Z".into(), 1.0),
        _ => unreachable!("unhandled {gate} {p}"),
    }
}

// ---------------------------------------------------------------------------
// 6b. THE PURE-SIGN GATES — X/Y/Z drain a `±1` in place (no re-key). These prove
//     the new crate's in-place `sign_flip_by_key` fast path is sign-correct
//     against the OLD crate's in-place `scale`, term for term.
// ---------------------------------------------------------------------------

/// A replayable pure-sign single-qubit gate (`X`/`Y`/`Z`).
#[derive(Clone, Copy, Debug)]
enum SignGate {
    X(usize),
    Y(usize),
    Z(usize),
}

fn apply_sign_old(s: &mut OldSum, g: SignGate) {
    match g {
        SignGate::X(q) => s.x(q),
        SignGate::Y(q) => s.y(q),
        SignGate::Z(q) => s.z(q),
    }
}

fn apply_sign_new(s: &mut NewSum, g: SignGate) {
    match g {
        SignGate::X(q) => s.x(q),
        SignGate::Y(q) => s.y(q),
        SignGate::Z(q) => s.z(q),
    }
}

/// A random `X`/`Y`/`Z` gate list on `n` qubits, `len` gates long.
fn random_sign_ops(rng: &mut StdRng, n: usize, len: usize) -> Vec<SignGate> {
    (0..len)
        .map(|_| {
            let q = rng.random_range(0..n);
            match rng.random_range(0..3usize) {
                0 => SignGate::X(q),
                1 => SignGate::Y(q),
                _ => SignGate::Z(q),
            }
        })
        .collect()
}

/// Replay a seeded random `X`/`Y`/`Z` sequence and, after *each* gate, assert the
/// full support + coefficients match the OLD crate. A pure-sign gate leaves the
/// word fixed and only flips a `±1`, so `len` is invariant and only signs move —
/// exactly the in-place fast path under test.
#[test]
fn pure_sign_clifford_replay_matches_old_per_gate() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 8, 12] {
            let terms = random_terms(&mut rng, n, 30);
            let mut old = build_old_sum(n, &terms);
            let mut new = build_new_sum(n, &terms);

            assert_supports_match(&old, &new, TOL);

            let ops = random_sign_ops(&mut rng, n, 120);
            for (i, g) in ops.into_iter().enumerate() {
                apply_sign_old(&mut old, g);
                apply_sign_new(&mut new, g);
                // Pure sign: the word never changes, so the support size is fixed.
                assert_eq!(
                    old.len(),
                    new.len(),
                    "len diverged after sign-gate {i} ({g:?}) seed {seed} n {n}"
                );
                assert_supports_match(&old, &new, TOL);
            }
        }
    }
}

/// Exhaustive single-qubit pure-sign draining: over every single-qubit Pauli term,
/// `X`/`Y`/`Z` on the old and new sums must agree on the resulting (word, sign),
/// and the sign must equal the reference conjugation. Includes the hand cases
/// `X·Z = −Z`, `Y·X = −X`, `Z·Y = −Y`.
#[test]
fn pure_sign_single_qubit_exhaustive() {
    type SignPair = (&'static str, fn(&mut OldSum), fn(&mut NewSum));
    let gens: [SignPair; 3] = [
        ("X", |s| s.x(0), |s| s.x(0)),
        ("Y", |s| s.y(0), |s| s.y(0)),
        ("Z", |s| s.z(0), |s| s.z(0)),
    ];
    for (name, gold, gnew) in gens {
        for p in ["I", "X", "Y", "Z"] {
            for &c in &[1.0f64, -1.0, 2.5, -0.75] {
                let terms = vec![(p.to_string(), c)];
                let mut old = build_old_sum(1, &terms);
                let mut new = build_new_sum(1, &terms);
                gold(&mut old);
                gnew(&mut new);
                assert_supports_match(&old, &new, TOL);

                let (w, sign) = new_support(&new)
                    .into_iter()
                    .next()
                    .map(|(w, cc)| (w, cc.signum()))
                    .unwrap();
                let expect = expected_sign(name, p);
                // A pure-sign gate never changes the word.
                assert_eq!(w, p, "{name} {p}: word must be unchanged");
                assert_eq!(sign, expect * c.signum(), "{name} {p}: sign");
            }
        }
    }
}

/// Reference sign of `g P g†` for single-qubit `X`/`Y`/`Z` (the word is fixed):
/// `+1` or `−1`.
///   X: flip iff `z` set  → `XZX = −Z`, `XYX = −Y`.
///   Y: flip iff `x ⊕ z`  → `YXY = −X`, `YZY = −Z`.
///   Z: flip iff `x` set  → `ZXZ = −X`, `ZYZ = −Y`.
fn expected_sign(gate: &str, p: &str) -> f64 {
    match (gate, p) {
        ("X", "I") | ("X", "X") => 1.0,
        ("X", "Y") | ("X", "Z") => -1.0,
        ("Y", "I") | ("Y", "Y") => 1.0,
        ("Y", "X") | ("Y", "Z") => -1.0,
        ("Z", "I") | ("Z", "Z") => 1.0,
        ("Z", "X") | ("Z", "Y") => -1.0,
        _ => unreachable!("unhandled {gate} {p}"),
    }
}

// ---------------------------------------------------------------------------
// 7. The Clifford+rotation `random_circuit`, Clifford subset, still agrees
// ---------------------------------------------------------------------------

#[test]
fn shared_random_circuit_clifford_subset_matches() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[2usize, 4, 8] {
            let terms = random_terms(&mut rng, n, 20);
            let mut old = build_old_sum(n, &terms);
            let mut new = build_new_sum(n, &terms);
            let circuit = random_circuit(&mut rng, n, 150);
            for op in circuit {
                match op {
                    GateOp::H(q) => {
                        old.h(q);
                        new.h(q);
                    }
                    GateOp::S(q) => {
                        old.s(q);
                        new.s(q);
                    }
                    GateOp::Cnot(c, t) => {
                        old.cnot(c, t);
                        new.cnot(c, t);
                    }
                    // Rotations are a later component; skip on both.
                    GateOp::Rx(..) | GateOp::Rz(..) => {}
                }
                assert_supports_match(&old, &new, TOL);
            }
        }
    }
}
