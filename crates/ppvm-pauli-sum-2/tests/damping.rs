// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use approx::assert_abs_diff_eq;
use ppvm_pauli_sum_2::{CoefficientThreshold, LossyPauliSum, PauliSum, PauliWord};
use ppvm_traits_2::AmplitudeDamping;

#[test]
fn amplitude_damping_merges_into_identity_after_original_updates() {
    let mut sum = LossyPauliSum::<f64>::new(1);
    for (word, c) in [("I", 2.0), ("X", 3.0), ("Y", 4.0), ("Z", 5.0), ("L", 6.0)] {
        sum += (word, c);
    }
    sum.amplitude_damping(0, 0.36);
    let mut expected = LossyPauliSum::<f64>::new(1);
    for (word, c) in [("I", 3.8), ("X", 2.4), ("Y", 3.2), ("Z", 3.2), ("L", 6.0)] {
        expected += (word, c);
    }
    assert_abs_diff_eq!(sum, expected, epsilon = 1e-12);
}

#[test]
fn full_damping_keeps_zero_terms() {
    let mut sum = PauliSum::<f64>::new(1);
    for word in ["I", "X", "Y", "Z"] {
        sum += (word, 1.0);
    }
    sum.amplitude_damping(0, 1.0);
    assert_eq!(sum.len(), 4);
    assert_eq!(sum.get(&PauliWord::from("I")), Some(2.0));
    for word in ["X", "Y", "Z"] {
        assert_eq!(sum.get(&PauliWord::from(word)), Some(0.0));
    }
}

#[test]
fn damping_does_not_reduce_or_truncate() {
    let mut sum = PauliSum::<f64, _>::with_policy(1, CoefficientThreshold { threshold: 10.0 });
    sum += ("Z", 1.0);
    sum.amplitude_damping(0, 0.0);
    assert_eq!(sum.get(&PauliWord::from("I")), Some(0.0));
    assert_eq!(sum.len(), 2);
    sum.reduce();
    assert_eq!(sum.len(), 1);
    sum.truncate();
    assert!(sum.is_empty());
}

#[test]
fn invalid_damping_target_panics_before_mutation() {
    use std::panic::{AssertUnwindSafe, catch_unwind};
    let mut sum = PauliSum::<f64>::new(1);
    sum += ("Z", 1.0);
    let original = sum.clone();
    assert!(catch_unwind(AssertUnwindSafe(|| sum.amplitude_damping(1, 0.1))).is_err());
    assert_eq!(sum, original);
}
