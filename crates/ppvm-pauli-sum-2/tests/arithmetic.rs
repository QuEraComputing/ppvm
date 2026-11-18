// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

// Borrowed and owned operator forms are both part of the API under test.
#![allow(clippy::op_ref)]

use approx::{AbsDiffEq, RelativeEq};
use num::Complex;
use ppvm_pauli_sum_2::{CoefficientThreshold, LossyPauliSum, NoPolicy, PauliSum};
use ppvm_traits_2::Clifford;
use std::panic::{AssertUnwindSafe, catch_unwind};

#[test]
fn addition_and_subtraction_accumulate_all_ownership_forms() {
    let left = PauliSum::<f64>::from_terms(1, [("X".into(), 2.0), ("Z".into(), 3.0)]);
    let right = PauliSum::<f64>::from_terms(1, [("X".into(), -2.0), ("Y".into(), 4.0)]);
    let added =
        PauliSum::<f64>::from_terms(1, [("X".into(), 0.0), ("Y".into(), 4.0), ("Z".into(), 3.0)]);
    let subtracted = PauliSum::<f64>::from_terms(
        1,
        [("X".into(), 4.0), ("Y".into(), -4.0), ("Z".into(), 3.0)],
    );
    assert_eq!(&left + &right, added);
    assert_eq!(left.clone() + &right, added);
    assert_eq!(&left + right.clone(), added);
    assert_eq!(left.clone() + right.clone(), added);
    assert_eq!(&left - &right, subtracted);
    assert_eq!(left.clone() - &right, subtracted);
    assert_eq!(&left - right.clone(), subtracted);
    assert_eq!(left.clone() - right.clone(), subtracted);
    let mut actual = left.clone();
    actual += &right;
    assert_eq!(actual, added);
    actual = left.clone();
    actual += right.clone();
    assert_eq!(actual, added);
    actual = left.clone();
    actual -= &right;
    assert_eq!(actual, subtracted);
    actual = left.clone();
    actual -= right;
    assert_eq!(actual, subtracted);
    actual += ("X", -4.0);
    assert_eq!(actual.get(&"X".into()), Some(0.0));
    assert_eq!(left.get(&"X".into()), Some(2.0));
}

#[test]
fn arithmetic_keeps_left_configuration_and_explicit_zeros() {
    let left = PauliSum::from_terms_with_policy(
        1,
        CoefficientThreshold { threshold: 10.0 },
        [("X".into(), 0.25)],
    )
    .preserving(["X".into()]);
    let right = PauliSum::from_terms_with_policy(
        1,
        CoefficientThreshold { threshold: 20.0 },
        [("Z".into(), 0.5)],
    )
    .preserving(["Z".into()]);
    for sum in [&left + &right, &left - &right] {
        assert_eq!(sum.len(), 2);
        assert_eq!(sum.policy(), left.policy());
        assert_eq!(sum.preserved_keys(), left.preserved_keys());
    }
    let mut cancelled = &left - &left;
    assert_eq!(cancelled.get(&"X".into()), Some(0.0));
    cancelled.reduce();
    assert!(cancelled.is_empty());
}

#[test]
fn width_errors_happen_before_mutation_including_empty_sums() {
    let original = PauliSum::<f64>::from_terms(1, [("X".into(), 2.0)]);
    let mut left = original.clone();
    for right in [
        PauliSum::<f64>::new(2),
        PauliSum::from_terms(2, [("XI".into(), 1.0)]),
    ] {
        assert!(catch_unwind(AssertUnwindSafe(|| left += &right)).is_err());
        assert_eq!(left, original);
        assert!(catch_unwind(AssertUnwindSafe(|| left -= &right)).is_err());
        assert_eq!(left, original);
    }
    assert!(catch_unwind(AssertUnwindSafe(|| left += ("XX", 1.0))).is_err());
    assert_eq!(left, original);
    assert!(catch_unwind(|| PauliSum::<f64>::new(1) + PauliSum::<f64>::new(2)).is_err());
}

#[test]
fn scalars_support_owned_and_borrowed_operands_without_reducing() {
    let sum = PauliSum::<f64>::from_terms(1, [("X".into(), 2.0), ("Z".into(), -3.0)]);
    let expected = PauliSum::<f64>::from_terms(1, [("X".into(), 4.0), ("Z".into(), -6.0)]);
    let scalar = 2.0;
    assert_eq!(&sum * scalar, expected);
    assert_eq!(&sum * &scalar, expected);
    assert_eq!(sum.clone() * scalar, expected);
    assert_eq!(sum.clone() * &scalar, expected);
    let mut actual = sum.clone();
    actual *= &scalar;
    assert_eq!(actual, expected);
    actual = sum.clone();
    actual *= scalar;
    assert_eq!(actual, expected);
    actual *= 0.0;
    assert_eq!(actual.len(), 2);
    assert_eq!(actual.get(&"X".into()), Some(0.0));
    assert_eq!(sum.get(&"X".into()), Some(2.0));
}

#[test]
fn complex_and_lossy_sums_use_the_same_arithmetic() {
    let a = LossyPauliSum::<Complex<f64>>::from_terms(1, [("L".into(), Complex::new(1.0, 2.0))]);
    let b = LossyPauliSum::<Complex<f64>>::from_terms(1, [("L".into(), Complex::new(2.0, -1.0))]);
    assert_eq!((&a + &b).get(&"L".into()), Some(Complex::new(3.0, 1.0)));
    assert_eq!((&a - &b).get(&"L".into()), Some(Complex::new(-1.0, 3.0)));
    let scalar = Complex::new(2.0, -1.0);
    let expected =
        LossyPauliSum::<Complex<f64>>::from_terms(1, [("L".into(), Complex::new(4.0, 3.0))]);
    assert_eq!(&a * scalar, expected);
    assert_eq!(&a * &scalar, expected);
    assert_eq!(a.clone() * scalar, expected);
    assert_eq!(a.clone() * &scalar, expected);
    let mut actual = a.clone();
    actual *= &scalar;
    assert_eq!(actual, expected);
    actual = a.clone();
    actual *= scalar;
    assert_eq!(actual, expected);
    assert!(actual.abs_diff_eq(&expected, 1e-12));
    assert!(actual.relative_eq(&expected, 1e-12, 1e-12));
}

#[test]
fn equality_compares_values_not_policy_or_scratch() {
    let a = PauliSum::from_terms_with_policy(
        1,
        CoefficientThreshold { threshold: 1.0 },
        [("X".into(), 2.0), ("Z".into(), 3.0)],
    );
    let mut b = PauliSum::with_capacity(1, CoefficientThreshold { threshold: 2.0 }, 64)
        .preserving(["X".into()]);
    b += ("Z", 3.0);
    b += ("X", 2.0);
    b.h(0);
    b.h(0);
    assert_eq!(a, b);
    assert!(a.relative_eq(&b, 0.0, 0.0));
    b += ("I", 0.0);
    assert_ne!(a, b);
    assert!(!a.abs_diff_eq(&b, 100.0));
    b.reduce();
    assert_eq!(a, b);
    assert_ne!(PauliSum::<f64>::new(1), PauliSum::<f64>::new(2));
}

#[test]
fn approximate_comparisons_require_matching_keys_and_width() {
    let a = PauliSum::<f64>::from_terms(1, [("X".into(), 1000.0)]);
    let b = PauliSum::<f64>::from_terms(1, [("X".into(), 1000.01)]);
    assert_ne!(a, b);
    assert!(a.abs_diff_eq(&b, 0.02));
    assert!(!a.abs_diff_eq(&b, 0.001));
    assert!(a.relative_eq(&b, 0.0, 1e-4));
    assert!(!a.relative_eq(&b, 0.0, 1e-6));
    let different_key = PauliSum::<f64>::from_terms(1, [("Z".into(), 1000.0)]);
    assert!(!a.abs_diff_eq(&different_key, 1e9));
    assert!(!a.relative_eq(&different_key, 1e9, 1e9));
    assert!(!PauliSum::<f64>::new(1).abs_diff_eq(&PauliSum::new(2), 1e9));
    let nan = PauliSum::<f64>::from_terms(1, [("X".into(), f64::NAN)]);
    assert_ne!(nan, nan);
    assert!(!nan.abs_diff_eq(&nan, 1e-6));
}

#[test]
fn formatting_sorts_by_weight_and_includes_zero_terms() {
    let sum = PauliSum::<f64>::from_terms(
        2,
        [
            ("XX".into(), -2.0),
            ("II".into(), 1.23456789),
            ("IX".into(), 0.0),
        ],
    );
    assert_eq!(format!("{sum}"), "1.235 * II + 0.000 * IX + -2.000 * XX");
    assert_eq!(
        format!("{sum:?}"),
        "1.23456789 * II + 0.00000000 * IX + -2.00000000 * XX"
    );
    let empty = PauliSum::<f64>::with_capacity(2, NoPolicy, 10);
    assert_eq!(format!("{empty}"), "");
    assert_eq!(format!("{empty:?}"), "");
    let lossy = LossyPauliSum::<f64>::from_terms(2, [("LX".into(), -0.5)]);
    assert_eq!(format!("{lossy}"), "-0.500 * LX");
    let complex = PauliSum::<Complex<f64>>::from_terms(1, [("X".into(), Complex::new(1.0, -2.0))]);
    assert_eq!(format!("{complex}"), "1.000-2.000i * X");
}
