// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the `ppvm-pauli-sum-2` graded engine + Clifford
//! propagation, driven through the public API (`PauliSum` alias).

use ppvm_pauli_sum_2::{
    CoefficientThreshold, CombinedPolicy, MaxPauliWeight, PauliSum, PauliWord, RekeyProducer, Sum,
};
use ppvm_traits_2::{Clifford, Coefficient};

fn pw(s: &str) -> PauliWord {
    PauliWord::from(s)
}

#[test]
fn from_terms_combines_and_keeps_zeros() {
    let mut sum: PauliSum = PauliSum::from_terms(
        2,
        [
            (pw("XI"), 1.0),
            (pw("IZ"), 2.0),
            (pw("XI"), 3.0),  // collides with the first → 4.0
            (pw("IZ"), -2.0), // cancels the second → stays, at exactly 0.0
        ],
    );
    assert_eq!(sum.n_sites(), 2);
    // Construction accumulates but does NOT reduce: the cancelled key survives at
    // 0.0, exactly as `n` applications of old's `sum += (word, coeff)` do.
    assert_eq!(sum.len(), 2);
    assert_eq!(sum.get(&pw("XI")), Some(4.0));
    assert_eq!(sum.get(&pw("IZ")), Some(0.0));
    assert!(sum.contains(&pw("IZ"), &0.0));

    // `reduce` is the caller-driven canonicalizer.
    sum.reduce();
    assert_eq!(sum.len(), 1);
    assert!(!sum.contains_key(&pw("IZ")));
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
    assert!(!sum.contains_key(&pw("X")));
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
    assert!(!sx.contains_key(&pw("X")));
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
    assert!(sum.contains_key(&pw("XI")));
    assert!(sum.contains_key(&pw("XY")));
    assert!(!sum.contains_key(&pw("IZ")));
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
    assert!(sum.contains_key(&pw("XI")));
}

#[test]
fn coefficient_mul_sign_is_correct() {
    // Guard the sign convention the Clifford drain relies on.
    assert_eq!(3.0_f64.mul_sign(-1), -3.0);
    assert_eq!(3.0_f64.mul_sign(1), 3.0);
}

// ---------------------------------------------------------------------------
// `apply` + `RekeyProducer` — the generic `TermProducer` path.
//
// Every hot gate family bypasses `apply` for a storage fast path
// (`rekey_bijective`, `rotate_in_place`, `scale_by_key`, `sign_flip_by_key`), so
// without these the producer path — the one a future branching channel, tableau
// gate or multiply producer will take — has no coverage at all.
// ---------------------------------------------------------------------------

/// String reversal is a bijection on Pauli words: `apply` must **replace** the
/// support with the image, not merge onto it (merging would leave both `k` and
/// `φ(k)` and double the support).
#[test]
fn apply_replaces_the_support_through_a_bijective_producer() {
    let mut sum: PauliSum = PauliSum::from_terms(3, [(pw("XIZ"), 2.0), (pw("IYI"), -1.5)]);
    sum.apply(RekeyProducer::new(|k: &PauliWord, c: &f64| {
        let reversed: String = k.to_string().chars().rev().collect();
        (pw(&reversed), *c)
    }));

    assert_eq!(sum.len(), 2, "a bijection preserves the support size");
    assert_eq!(sum.get(&pw("ZIX")), Some(2.0));
    assert_eq!(sum.get(&pw("IYI")), Some(-1.5));
    assert!(!sum.contains_key(&pw("XIZ")), "the pre-image must be gone");

    // Applying it twice is the identity — and exercises the store's reused batch
    // on a second call.
    sum.apply(RekeyProducer::new(|k: &PauliWord, c: &f64| {
        let reversed: String = k.to_string().chars().rev().collect();
        (pw(&reversed), *c)
    }));
    assert_eq!(sum.get(&pw("XIZ")), Some(2.0));
    assert_eq!(sum.get(&pw("IYI")), Some(-1.5));
}

/// Colliding produced keys are **accumulated** (the free-module merge), and an
/// exact cancellation stays in the support at `0.0` — `apply` runs no `reduce`.
#[test]
fn apply_accumulates_collisions_and_keeps_exact_zeros() {
    let mut sum: PauliSum = PauliSum::from_terms(2, [(pw("XI"), 1.5), (pw("IZ"), -1.5)]);
    // Collapse everything onto one key: not a bijection, so the merge must sum.
    sum.apply(RekeyProducer::new(|_: &PauliWord, c: &f64| (pw("II"), *c)));
    assert_eq!(sum.len(), 1);
    assert_eq!(sum.get(&pw("II")), Some(0.0), "the cancellation is kept");
}

/// `apply` never truncates: a produced term below the policy's floor survives
/// until the caller asks (behavioural contract 1).
#[test]
fn apply_does_not_truncate() {
    let mut sum: PauliSum<f64, CoefficientThreshold> =
        PauliSum::<f64, CoefficientThreshold>::from_terms_with_policy(
            2,
            CoefficientThreshold { threshold: 0.5 },
            [(pw("XI"), 1.0)],
        );
    sum.apply(RekeyProducer::new(|k: &PauliWord, c: &f64| {
        (k.clone(), *c * 0.1)
    }));
    assert_eq!(sum.get(&pw("XI")), Some(0.1), "sub-threshold term survives");
    sum.truncate();
    assert!(!sum.contains_key(&pw("XI")), "…until the caller truncates");
}

// ---------------------------------------------------------------------------
// Sum + sum: free-module addition (pointwise coefficient ADDITION).
//
// The Lean-adjudicated correction of old's `extend`-based right-biased
// overwrite — see `Sum::add_sum`.
// ---------------------------------------------------------------------------

#[test]
fn sum_addition_accumulates_shared_keys() {
    let a: PauliSum = PauliSum::from_terms(2, [(pw("ZZ"), 1.0), (pw("XI"), 4.0)]);
    let b: PauliSum = PauliSum::from_terms(2, [(pw("ZZ"), 2.0), (pw("IY"), -1.0)]);

    let mut got = a.clone();
    got += &b;
    // Old's `extend` would leave `ZZ` at 2.0 (a right-biased overwrite); free
    // module addition is pointwise addition.
    assert_eq!(got.get(&pw("ZZ")), Some(3.0));
    assert_eq!(got.get(&pw("XI")), Some(4.0));
    assert_eq!(got.get(&pw("IY")), Some(-1.0));
    assert_eq!(got.len(), 3);

    // The by-value operator and the `&a + &b` form agree.
    let mut by_value = a.clone();
    by_value += b.clone();
    assert!(by_value == got);
    assert!(&a + &b == got);
}

#[test]
fn sum_addition_keeps_exact_cancellations_and_does_not_truncate() {
    let mut a: PauliSum<f64, CoefficientThreshold> =
        PauliSum::<f64, CoefficientThreshold>::from_terms_with_policy(
            1,
            CoefficientThreshold { threshold: 0.5 },
            [(pw("Z"), 1.0)],
        );
    let b: PauliSum<f64, CoefficientThreshold> =
        PauliSum::<f64, CoefficientThreshold>::from_terms_with_policy(
            1,
            CoefficientThreshold { threshold: 0.5 },
            [(pw("Z"), -1.0), (pw("X"), 0.1)],
        );
    a += &b;
    assert_eq!(
        a.get(&pw("Z")),
        Some(0.0),
        "cancellation stays in the support"
    );
    assert_eq!(a.get(&pw("X")), Some(0.1), "sub-threshold term stays too");
    assert_eq!(a.len(), 2);
}

#[test]
fn sum_addition_is_commutative_and_has_the_empty_sum_as_identity() {
    let a: PauliSum = PauliSum::from_terms(2, [(pw("ZZ"), 1.0), (pw("XI"), 4.0)]);
    let b: PauliSum = PauliSum::from_terms(2, [(pw("ZZ"), 2.0), (pw("IY"), -1.0)]);
    assert!(&a + &b == &b + &a);

    let zero: PauliSum = PauliSum::new(2);
    assert!(&a + &zero == a);
}
