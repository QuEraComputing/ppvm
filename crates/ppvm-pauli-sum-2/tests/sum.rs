// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::Complex;
use ppvm_pauli_sum_2::{
    CoefficientThreshold, LossyPauliSum, NoPolicy, PauliSum, PauliWord, Policy,
};
use ppvm_traits_2::Retain;
use std::panic::{AssertUnwindSafe, catch_unwind};

#[test]
fn empty_sums_retain_their_width() {
    let mut sum = PauliSum::<f64>::with_capacity(3, NoPolicy, 8);
    assert_eq!(sum.n_sites(), 3);
    assert!(sum.is_empty());
    sum.add_term("XYZ".into(), 2.0);
    let mut cloned = sum.clone();
    sum.clear();
    assert_eq!(sum.n_sites(), 3);
    assert!(sum.is_empty());
    assert_eq!(cloned.get(&"XYZ".into()), Some(2.0));
    cloned.add_term("XYZ".into(), -2.0);
    cloned.reduce();
    assert!(cloned.is_empty());
    assert_eq!(cloned.n_sites(), 3);

    let mut zero_width = PauliSum::<f64>::new(0);
    zero_width.add_term("".into(), 1.0);
    assert_eq!(zero_width.get(&"".into()), Some(1.0));
}

#[test]
fn construction_and_insertion_accumulate_without_reducing() {
    let x = PauliWord::from("X");
    let z = PauliWord::from("Z");
    let mut sum = PauliSum::<f64>::from_terms(1, [(x, 2.0), (x, -2.0), (z, 3.0)]);
    assert_eq!(sum.len(), 2);
    assert!(sum.contains(&x, &0.0));
    assert!(!sum.contains(&x, &2.0));
    assert!(sum.contains_key(&x));
    assert!(!sum.contains_key(&"Y".into()));
    assert_eq!(sum.policy(), &NoPolicy);
    sum.truncate();
    assert_eq!(sum.len(), 2);
    sum.add_term(z, 1.0);
    sum.scale(&2.0);
    assert_eq!(sum.get(&z), Some(8.0));
    sum.reduce();
    assert_eq!(sum.iter().collect::<Vec<_>>(), vec![(z, 8.0)]);
    let mut visited = Vec::new();
    sum.for_each_ref(|key, coeff| visited.push((*key, *coeff)));
    assert_eq!(visited, vec![(z, 8.0)]);
}

#[test]
fn lossy_alias_uses_the_same_core_operations() {
    let mut sum = LossyPauliSum::<f64>::from_terms_with_policy(
        2,
        NoPolicy,
        [("XL".into(), 2.0), ("XI".into(), 3.0), ("XL".into(), -1.0)],
    );
    sum.scale(&2.0);
    assert_eq!(sum.len(), 2);
    assert_eq!(sum.get(&"XL".into()), Some(2.0));
    assert_eq!(sum.get(&"XI".into()), Some(6.0));
    let other = LossyPauliSum::<f64>::from_terms(2, [("XL".into(), 3.0)]);
    assert_eq!(sum.overlap(&other), 6.0);
}

#[test]
fn complex_pairing_accepts_different_policies_without_truncating() {
    let x = PauliWord::from("X");
    let mut left = PauliSum::<Complex<f64>>::from_terms(1, [(x, Complex::new(1.0, 2.0))]);
    left.scale(&Complex::new(0.0, 1.0));
    assert_eq!(left.get(&x), Some(Complex::new(-2.0, 1.0)));
    let right = PauliSum::from_terms_with_policy(
        1,
        CoefficientThreshold { threshold: 10.0 },
        [(x, Complex::new(3.0, 4.0))],
    );
    assert_eq!(left.overlap(&right), Complex::new(-10.0, -5.0));
    assert_eq!(right.overlap(&left), Complex::new(-10.0, -5.0));
    assert_eq!(left.hermitian_overlap(&right), Complex::new(-2.0, -11.0));
    assert_eq!(right.hermitian_overlap(&left), Complex::new(-2.0, 11.0));
}

#[test]
fn invalid_insertions_leave_the_sum_unchanged() {
    let mut sum = PauliSum::<f64>::from_terms(2, [("XI".into(), 3.0)]);
    let result = catch_unwind(AssertUnwindSafe(|| sum.add_term("X".into(), 2.0)));
    assert!(result.is_err());
    assert_eq!(sum.len(), 1);
    assert_eq!(sum.get(&"XI".into()), Some(3.0));

    assert!(catch_unwind(|| PauliSum::<f64>::from_terms(2, [("X".into(), 1.0)])).is_err());
    assert!(catch_unwind(|| LossyPauliSum::<f64>::from_terms(2, [("L".into(), 1.0)])).is_err());
}

#[test]
fn overlap_rejects_different_widths_even_for_empty_sums() {
    let left = PauliSum::<f64>::new(1);
    let right = PauliSum::<f64, _>::with_policy(2, CoefficientThreshold::default());
    assert!(catch_unwind(|| left.overlap(&right)).is_err());
    assert!(catch_unwind(|| left.hermitian_overlap(&right)).is_err());
}

#[test]
fn explicit_policies_need_neither_default_nor_clone() {
    struct KeepAbove(f64);

    impl Policy<PauliWord, f64> for KeepAbove {
        fn truncate<M: Retain<PauliWord, f64>>(&self, storage: &mut M) {
            storage.retain(|_, coeff| coeff.abs() >= self.0);
        }
    }

    let mut sums = [
        PauliSum::with_policy(1, KeepAbove(2.0)),
        PauliSum::with_capacity(1, KeepAbove(2.0), 4),
        PauliSum::from_terms_with_policy(1, KeepAbove(2.0), []),
    ];
    for sum in &mut sums {
        sum.add_term("X".into(), 1.0);
        sum.add_term("Z".into(), 3.0);
        assert_eq!(sum.len(), 2);
        sum.truncate();
        assert_eq!(sum.len(), 1);
        assert_eq!(sum.get(&"Z".into()), Some(3.0));
        sum.clear();
        assert!(sum.is_empty());
        assert_eq!(sum.policy().0, 2.0);
    }
}
