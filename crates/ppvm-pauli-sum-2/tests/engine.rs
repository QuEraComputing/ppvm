// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the `ppvm-pauli-sum-2` graded engine + Clifford
//! propagation, driven through the public API (`PauliSum` alias).

use ppvm_pauli_sum_2::{
    CoefficientThreshold, CombinedPolicy, MaxPauliWeight, PauliSum, PauliWord, Sum,
};
use ppvm_traits_2::{Clifford, Coefficient};

fn pw(s: &str) -> PauliWord {
    PauliWord::from(s)
}

#[test]
fn from_terms_combines_and_reduces() {
    let sum: PauliSum = PauliSum::from_terms(
        2,
        [
            (pw("XI"), 1.0),
            (pw("IZ"), 2.0),
            (pw("XI"), 3.0),  // collides with the first → 4.0
            (pw("IZ"), -2.0), // cancels the second → dropped by reduce
        ],
    );
    assert_eq!(sum.n_sites(), 2);
    assert_eq!(sum.len(), 1);
    assert_eq!(sum.get(&pw("XI")), Some(4.0));
    assert!(!sum.contains(&pw("IZ")));
}

#[test]
fn ghz_propagation_zz_to_iz() {
    // The old ppvm-pauli-sum doctest: ZZ propagated (Heisenberg) through
    // CNOT(0,1) then H(0) collapses to IZ with coefficient 1.0.
    let mut sum: PauliSum = PauliSum::from_terms(2, [(pw("ZZ"), 1.0)]);
    sum.cnot(0, 1);
    sum.h(0);
    assert_eq!(sum.len(), 1);
    assert_eq!(sum.get(&pw("IZ")), Some(1.0));
}

#[test]
fn clifford_rekey_replaces_not_merges() {
    // X on q0 conjugated by H(0) becomes Z on q0; the old X key must be gone.
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("X"), 1.0)]);
    sum.h(0);
    assert_eq!(sum.len(), 1, "re-key must replace, not leave a stale key");
    assert_eq!(sum.get(&pw("Z")), Some(1.0));
    assert!(!sum.contains(&pw("X")));
}

#[test]
fn single_qubit_signs_match_conjugation_tables() {
    // H Y H = -Y (word unchanged, sign flips).
    let mut hy: PauliSum = PauliSum::from_terms(1, [(pw("Y"), 1.0)]);
    hy.h(0);
    assert_eq!(hy.get(&pw("Y")), Some(-1.0));

    // Z X Z = -X.
    let mut zx: PauliSum = PauliSum::from_terms(1, [(pw("X"), 1.0)]);
    zx.z(0);
    assert_eq!(zx.get(&pw("X")), Some(-1.0));

    // S X S(dagger convention) = -Y in this crate's backward Heisenberg run.
    let mut sx: PauliSum = PauliSum::from_terms(1, [(pw("X"), 1.0)]);
    sx.s(0);
    assert_eq!(sx.get(&pw("Y")), Some(-1.0));
    assert!(!sx.contains(&pw("X")));
}

#[test]
fn cnot_sign_matches_table() {
    // +XZ under CNOT(0,1) becomes -YY (from the phased-word conjugation table).
    let mut sum: PauliSum = PauliSum::from_terms(2, [(pw("XZ"), 1.0)]);
    sum.cnot(0, 1);
    assert_eq!(sum.len(), 1);
    assert_eq!(sum.get(&pw("YY")), Some(-1.0));
}

#[test]
fn cz_sign_matches_table() {
    // +XY under CZ(0,1) becomes -YX.
    let mut sum: PauliSum = PauliSum::from_terms(2, [(pw("XY"), 1.0)]);
    sum.cz(0, 1);
    assert_eq!(sum.get(&pw("YX")), Some(-1.0));
}

#[test]
fn h_is_involutive_on_coefficients() {
    let mut sum: PauliSum = PauliSum::from_terms(1, [(pw("Y"), 3.0)]);
    sum.h(0);
    sum.h(0);
    assert_eq!(sum.get(&pw("Y")), Some(3.0));
}

#[test]
fn scale_multiplies_all_coefficients() {
    let mut sum: PauliSum = PauliSum::from_terms(2, [(pw("XI"), 2.0), (pw("IZ"), -1.5)]);
    sum.scale(&2.0);
    assert_eq!(sum.get(&pw("XI")), Some(4.0));
    assert_eq!(sum.get(&pw("IZ")), Some(-3.0));
}

#[test]
fn max_pauli_weight_policy_truncates() {
    let mut sum: Sum<_, MaxPauliWeight> = PauliSum::<f64, MaxPauliWeight>::from_terms_with_policy(
        3,
        MaxPauliWeight(1),
        [(pw("XII"), 1.0), (pw("XYI"), 1.0), (pw("XYZ"), 1.0)],
    );
    sum.truncate();
    assert_eq!(sum.len(), 1);
    assert_eq!(sum.get(&pw("XII")), Some(1.0));
}

#[test]
fn coefficient_threshold_keeps_boundary() {
    let mut sum: Sum<_, CoefficientThreshold> =
        PauliSum::<f64, CoefficientThreshold>::from_terms_with_policy(
            2,
            CoefficientThreshold { threshold: 1.0 },
            [(pw("XI"), 1.0), (pw("IZ"), 0.5), (pw("XY"), 2.0)],
        );
    sum.truncate();
    // 1.0 is exactly the threshold → kept (>=). 0.5 dropped.
    assert!(sum.contains(&pw("XI")));
    assert!(sum.contains(&pw("XY")));
    assert!(!sum.contains(&pw("IZ")));
}

#[test]
fn combined_policy_applies_both() {
    type P = CombinedPolicy<MaxPauliWeight, CoefficientThreshold>;
    let policy: P = CombinedPolicy(MaxPauliWeight(1), CoefficientThreshold { threshold: 1.0 });
    let mut sum: Sum<_, P> = PauliSum::<f64, P>::from_terms_with_policy(
        2,
        policy,
        [
            (pw("XI"), 2.0),  // weight 1, |c|>=1 → kept
            (pw("XY"), 5.0),  // weight 2 → dropped by MaxPauliWeight
            (pw("IZ"), 0.25), // weight 1 but |c|<1 → dropped by threshold
        ],
    );
    sum.truncate();
    assert_eq!(sum.len(), 1);
    assert!(sum.contains(&pw("XI")));
}

#[test]
fn coefficient_mul_sign_is_correct() {
    // Guard the sign convention the Clifford drain relies on.
    assert_eq!(3.0_f64.mul_sign(-1), -3.0);
    assert_eq!(3.0_f64.mul_sign(1), 3.0);
}
