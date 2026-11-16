// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::Complex;
use ppvm_pauli_sum_2::{
    CoefficientThreshold, CombinedPolicy, LossyPauliSum, MaxLossWeight, MaxPauliWeight, PauliSum,
};
use std::panic::{AssertUnwindSafe, catch_unwind};

#[test]
fn threshold_is_inclusive_and_truncation_is_explicit() {
    let mut sum = PauliSum::from_terms_with_policy(
        1,
        CoefficientThreshold { threshold: 0.5 },
        [
            ("I".into(), 0.25),
            ("I".into(), 0.25),
            ("X".into(), -0.5),
            ("Y".into(), 0.25),
            ("Z".into(), 1.0),
            ("Z".into(), -1.0),
        ],
    );
    assert_eq!(sum.len(), 4);
    sum.truncate();
    assert_eq!(sum.len(), 2);
    assert_eq!(sum.get(&"I".into()), Some(0.5));
    assert_eq!(sum.get(&"X".into()), Some(-0.5));
    sum.scale(&0.5);
    assert_eq!(sum.len(), 2);
    sum.truncate();
    assert!(sum.is_empty());
}

#[test]
fn threshold_uses_complex_magnitude() {
    let mut sum = PauliSum::from_terms_with_policy(
        1,
        CoefficientThreshold { threshold: 5.0 },
        [
            ("X".into(), Complex::new(3.0, 4.0)),
            ("Y".into(), Complex::new(0.0, -4.0)),
            ("Z".into(), Complex::new(-6.0, 0.0)),
        ],
    );
    sum.truncate();
    assert_eq!(sum.len(), 2);
    assert!(sum.contains_key(&"X".into()));
    assert!(sum.contains_key(&"Z".into()));
}

#[test]
fn pauli_weight_counts_nonidentity_and_lost_sites() {
    let terms = [("II".into(), 1.0), ("XI".into(), 2.0), ("XY".into(), 3.0)];
    for (limit, expected) in [(0, 1), (1, 2), (usize::MAX, 3)] {
        let mut sum = PauliSum::from_terms_with_policy(2, MaxPauliWeight(limit), terms);
        sum.truncate();
        assert_eq!(sum.len(), expected);
        assert_eq!(sum.get(&"II".into()), Some(1.0));
    }

    let mut sum = LossyPauliSum::from_terms_with_policy(
        2,
        MaxPauliWeight(1),
        [("LI".into(), 1.0), ("LX".into(), 2.0)],
    );
    sum.truncate();
    assert_eq!(sum.len(), 1);
    assert_eq!(sum.get(&"LI".into()), Some(1.0));
}

#[test]
fn loss_weight_only_counts_lost_sites() {
    let terms = [("XY".into(), 1.0), ("XL".into(), 2.0), ("LL".into(), 3.0)];
    for (limit, expected) in [(0, 1), (1, 2), (usize::MAX, 3)] {
        let mut sum = LossyPauliSum::from_terms_with_policy(2, MaxLossWeight(limit), terms.clone());
        sum.truncate();
        assert_eq!(sum.len(), expected);
        assert_eq!(sum.get(&"XY".into()), Some(1.0));
    }

    let mut ordinary = PauliSum::from_terms_with_policy(2, MaxLossWeight(0), [("XY".into(), 1.0)]);
    ordinary.truncate();
    assert_eq!(ordinary.len(), 1);
}

#[test]
fn default_limits_match_the_source_policies() {
    assert_eq!(CoefficientThreshold::default().threshold, 1e-12);
    assert_eq!(MaxPauliWeight::default().0, usize::MAX);
    assert_eq!(MaxLossWeight::default().0, 10);
    let mut sum = LossyPauliSum::<f64, MaxLossWeight>::from_terms(
        11,
        [("LLLLLLLLLLX".into(), 1.0), ("LLLLLLLLLLL".into(), 2.0)],
    );
    sum.truncate();
    assert_eq!(sum.len(), 1);
    assert_eq!(sum.get(&"LLLLLLLLLLX".into()), Some(1.0));
}

#[test]
fn nested_policies_filter_independently_and_preserve_exact_coefficients() {
    let mut sum = LossyPauliSum::from_terms_with_policy(
        3,
        CombinedPolicy(
            CoefficientThreshold { threshold: 0.5 },
            CombinedPolicy(MaxPauliWeight(2), MaxLossWeight(1)),
        ),
        [
            ("XLI".into(), 0.5),
            ("XII".into(), 0.25),
            ("XYZ".into(), 2.0),
            ("LLI".into(), 3.0),
            ("LLL".into(), 0.0),
        ],
    )
    .preserving(["LLL".into(), "LLL".into(), "YII".into()]);
    assert_eq!(sum.preserved_keys().len(), 2);
    for _ in 0..2 {
        sum.truncate();
        assert_eq!(sum.len(), 2);
        assert_eq!(sum.get(&"XLI".into()), Some(0.5));
        assert_eq!(sum.get(&"LLL".into()), Some(0.0));
        assert!(!sum.contains_key(&"YII".into()));
    }
    sum.add_term("YII".into(), 0.125);
    sum.truncate();
    assert_eq!(sum.get(&"YII".into()), Some(0.125));
    sum.reduce();
    assert!(!sum.contains_key(&"LLL".into()));
    sum.clear_preserved_keys();
    sum.truncate();
    assert_eq!(sum.len(), 1);
    assert_eq!(sum.get(&"XLI".into()), Some(0.5));
}

#[test]
fn clearing_terms_keeps_exemptions_and_cloning_is_independent() {
    let mut sum = PauliSum::from_terms_with_policy(1, MaxPauliWeight(0), [("X".into(), 2.0)])
        .preserving(["X".into()]);
    let mut cloned = sum.clone();
    cloned.clear_preserved_keys();
    cloned.truncate();
    assert!(cloned.is_empty());
    sum.clear();
    assert!(sum.is_empty());
    assert_eq!(sum.preserved_keys().len(), 1);
    sum.add_term("X".into(), 3.0);
    sum.truncate();
    assert_eq!(sum.get(&"X".into()), Some(3.0));
}

#[test]
fn invalid_preserved_width_leaves_sum_unchanged() {
    let mut sum = PauliSum::from_terms_with_policy(2, MaxPauliWeight(0), [("XI".into(), 1.0)])
        .preserving(["XI".into()]);
    assert!(catch_unwind(AssertUnwindSafe(|| sum.preserve_key("X".into()))).is_err());
    assert_eq!(sum.preserved_keys().len(), 1);
    sum.truncate();
    assert_eq!(sum.get(&"XI".into()), Some(1.0));
}
