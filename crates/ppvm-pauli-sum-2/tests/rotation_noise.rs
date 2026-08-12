// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the non-Clifford `RotationOne` branch and the diagonal
//! `PauliError` channel on `PauliSum<f64>`, ported from
//! `ppvm-pauli-sum/src/sum/rot1.rs` and `sum/noise.rs`.

use ppvm_pauli_sum_2::{PauliSum, PauliWord};
use ppvm_traits_2::{PauliError, RotationOne};
use rand::SeedableRng;
use rand::rngs::SmallRng;

/// The `-2` sum backends are density-matrix-like: the noise channels are
/// analytic coefficient scalings and never draw. Call sites thread this
/// fixed-seed RNG only to satisfy the injected-RNG trait surface.
fn rng() -> SmallRng {
    SmallRng::seed_from_u64(0)
}

fn pw(s: &str) -> PauliWord {
    PauliWord::from(s)
}

fn approx(sum: &PauliSum, key: &str, want: f64) {
    let got = sum.get(&pw(key)).unwrap_or(0.0);
    assert!(
        (got - want).abs() < 1e-12,
        "coeff at {key}: got {got}, want {want}"
    );
}

// --- RotationOne: rx/ry/rz on each Pauli, against the old rot1 test tables. ----

#[test]
fn rx_branches_match_old_tables() {
    let theta = 2.1_f64;
    let (s, c) = (theta.sin(), theta.cos());

    // X commutes with axis X → unchanged.
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("X"), 1.0)]);
    sum.rx(0, theta);
    assert_eq!(sum.len(), 1);
    approx(&sum, "X", 1.0);

    // Y anticommutes: Y → cos·Y + (−sin)·Z.
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("Y"), 1.0)]);
    sum.rx(0, theta);
    approx(&sum, "Y", c);
    approx(&sum, "Z", -s);

    // Z anticommutes: Z → cos·Z + sin·Y.
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("Z"), 1.0)]);
    sum.rx(0, theta);
    approx(&sum, "Z", c);
    approx(&sum, "Y", s);

    // I commutes → unchanged.
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("I"), 1.0)]);
    sum.rx(0, theta);
    approx(&sum, "I", 1.0);
}

#[test]
fn ry_branches_match_old_tables() {
    let theta = 2.1_f64;
    let (s, c) = (theta.sin(), theta.cos());

    // X anticommutes: X → cos·X + sin·Z.
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("X"), 1.0)]);
    sum.ry(0, theta);
    approx(&sum, "X", c);
    approx(&sum, "Z", s);

    // Y commutes → unchanged.
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("Y"), 1.0)]);
    sum.ry(0, theta);
    assert_eq!(sum.len(), 1);
    approx(&sum, "Y", 1.0);

    // Z anticommutes: Z → cos·Z + (−sin)·X.
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("Z"), 1.0)]);
    sum.ry(0, theta);
    approx(&sum, "Z", c);
    approx(&sum, "X", -s);
}

#[test]
fn rz_branches_match_old_tables() {
    let theta = 2.1_f64;
    let (s, c) = (theta.sin(), theta.cos());

    // X anticommutes: X → cos·X + (−sin)·Y.
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("X"), 1.0)]);
    sum.rz(0, theta);
    approx(&sum, "X", c);
    approx(&sum, "Y", -s);

    // Y anticommutes: Y → cos·Y + sin·X.
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("Y"), 1.0)]);
    sum.rz(0, theta);
    approx(&sum, "Y", c);
    approx(&sum, "X", s);

    // Z commutes → unchanged.
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("Z"), 1.0)]);
    sum.rz(0, theta);
    assert_eq!(sum.len(), 1);
    approx(&sum, "Z", 1.0);
}

#[test]
fn rotation_branch_merges_colliding_key() {
    // A branch key that collides with an existing term must merge (unlike a
    // bijective Clifford re-key). Start with cos-heavy Z and a Y that rx will
    // branch onto: rx(Z) contributes sin·Y, which merges with the stored Y.
    let theta = 0.5_f64;
    let (s, c) = (theta.sin(), theta.cos());
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("Z"), 1.0), (pw("Y"), 1.0)]);
    sum.rx(0, theta);
    // Z → cos·Z + sin·Y ; Y → cos·Y + (−sin)·Z.
    approx(&sum, "Z", c - s);
    approx(&sum, "Y", c + s);
    assert_eq!(sum.len(), 2);
}

// --- PauliError: diagonal per-term eigenvalue scale. ---------------------------

#[test]
fn pauli_error_scales_each_pauli_by_its_eigenvalue() {
    // p = [pX, pY, pZ] = [0.0, 0.1, 0.2].
    // λ_X = 1 − 2(pY+pZ) = 0.4 ; λ_Y = 1 − 2(pX+pZ) = 0.6 ; λ_Z = 1 − 2(pX+pY) = 0.8.
    let p = [0.0_f64, 0.1, 0.2];

    let mut sum: PauliSum = PauliSum::from_terms(
        1,
        [
            (pw("I"), 1.0),
            (pw("X"), 1.0),
            (pw("Y"), 1.0),
            (pw("Z"), 1.0),
        ],
    );
    sum.pauli_error(0, p, &mut rng());
    approx(&sum, "I", 1.0);
    approx(&sum, "X", 0.4);
    approx(&sum, "Y", 0.6);
    approx(&sum, "Z", 0.8);
    // Diagonal: no key moves, no new terms.
    assert_eq!(sum.len(), 4);
}

#[test]
fn pauli_error_zero_eigenvalue_keeps_the_term_matching_old() {
    // λ_X = 1 − 2(pY + pZ) = 1 − 2(0.25 + 0.25) = 0 exactly, so the X-term is
    // scaled to 0 — and must STAY in the support with coefficient 0.
    //
    // This pins behaviour parity with the old crate, which has no `reduce` and
    // no drop-zero path anywhere: `PauliSum::scale` can only mutate, never
    // remove, and old's exact-map `PartialEq` depends on zero terms surviving
    // (`ppvm-pauli-sum/tests/loss.rs::test_reset_channel` asserts equality with
    // a clone after a channel multiplies coefficients by 0).
    //
    // An earlier revision of this crate dropped the term here, justified by a
    // "reduced-canonical-form invariant" that is this design's own invention
    // rather than old's behaviour — a divergence under the behaviour-preserving
    // prime directive (gap `ps2.zero.behaviour`). Restored to match old.
    let p = [0.0_f64, 0.25, 0.25];
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("X"), 1.0), (pw("Z"), 1.0)]);
    sum.pauli_error(0, p, &mut rng());
    assert_eq!(
        sum.get(&pw("X")),
        Some(0.0),
        "zeroed X-term must REMAIN in the support (old never removes it)"
    );
    assert_eq!(sum.len(), 2, "both terms remain; nothing is removed");
    // λ_Z = 1 − 2(pX + pY) = 1 − 2(0.25) = 0.5 for the Z-term.
    approx(&sum, "Z", 0.5);
}

#[test]
fn pauli_error_only_touches_target_qubit() {
    // On a 2-qubit term, only the Pauli at the target qubit selects the factor.
    // Symmetric depolarizing-ish: λ for any non-I Pauli = 1 − 2(0.1) = 0.8.
    let p = [0.05_f64, 0.05, 0.05];
    let mut sum: PauliSum = PauliSum::from_terms(2, [(pw("XZ"), 1.0), (pw("IZ"), 1.0)]);
    sum.pauli_error(0, p, &mut rng());
    // Qubit 0 is X on the first term (scaled by 0.8) and I on the second (kept).
    approx(&sum, "XZ", 0.8);
    approx(&sum, "IZ", 1.0);
}
