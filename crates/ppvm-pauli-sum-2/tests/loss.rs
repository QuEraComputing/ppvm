// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use approx::assert_abs_diff_eq;
use num::Complex;
use ppvm_pauli_sum_2::{CoefficientThreshold, LossyPauliSum, LossyPauliWord};
use ppvm_traits_2::{CorrelatedLossChannel, LossChannel, ResetLossChannel};

#[test]
fn single_loss_merges_branches_after_scaling_original_terms() {
    let mut sum = LossyPauliSum::<f64>::new(1);
    sum += ("L", 2.0);
    sum += ("I", 3.0);
    sum += ("X", 4.0);
    sum.loss_channel(0, 0.2);
    let mut expected = LossyPauliSum::<f64>::new(1);
    expected += ("L", 2.0);
    expected += ("I", 2.8);
    expected += ("X", 3.2);
    assert_abs_diff_eq!(sum, expected, epsilon = 1e-12);
}

#[test]
fn correlated_loss_preserves_distinct_branch_weights_and_collisions() {
    let mut sum = LossyPauliSum::<Complex<f64>>::new(2);
    let c = Complex::new(1.0, 2.0);
    for (word, v) in [("LL", 2.0), ("IL", 3.0), ("LI", 5.0), ("II", 7.0)] {
        sum += (word, c * v);
    }
    sum.correlated_loss_channel(0, 1, [0.1.into(), 0.2.into(), 0.3.into()]);
    let mut expected = LossyPauliSum::<Complex<f64>>::new(2);
    for (word, v) in [("LL", 2.0), ("IL", 2.7), ("LI", 4.1), ("II", 5.3)] {
        expected += (word, c * v);
    }
    assert_abs_diff_eq!(sum, expected, epsilon = 1e-12);
}

#[test]
fn loss_maps_preserve_the_full_identity_including_lost_sectors() {
    // I projects onto the live subspace; I+L is the full one-site identity.
    let mut identity = LossyPauliSum::<f64>::new(2);
    for word in ["II", "IL", "LI", "LL"] {
        identity += (word, 1.0);
    }
    let mut sum = identity.clone();
    sum.loss_channel(0, 0.2);
    sum.correlated_loss_channel(0, 1, [0.1, 0.2, 0.3]);
    sum.reset_loss_channel(1);
    assert_abs_diff_eq!(sum, identity, epsilon = 1e-12);
}

#[test]
fn reset_loss_accumulates_i_and_z_into_l_without_reprocessing_them() {
    let mut sum = LossyPauliSum::<f64>::new(1);
    for (word, c) in [("I", 2.0), ("Z", 3.0), ("L", 11.0), ("X", 7.0)] {
        sum += (word, c);
    }
    sum.reset_loss_channel(0);
    let mut expected = LossyPauliSum::<f64>::new(1);
    for (word, c) in [("I", 2.0), ("Z", 3.0), ("L", 5.0), ("X", 7.0)] {
        expected += (word, c);
    }
    assert_eq!(sum, expected);
}

#[test]
fn loss_keeps_zero_branches_and_does_not_apply_policy() {
    let mut sum = LossyPauliSum::<f64, _>::with_policy(2, CoefficientThreshold { threshold: 10.0 });
    sum += ("LL", 1.0);
    sum.correlated_loss_channel(0, 1, [0.0; 3]);
    assert_eq!(sum.len(), 4);
    assert_eq!(sum.get(&LossyPauliWord::from("II")), Some(0.0));
    sum.reduce();
    assert_eq!(sum.len(), 1);
    sum.reset_loss_channel(0);
    assert_eq!(sum.len(), 1);
    assert_eq!(sum.get(&LossyPauliWord::from("LL")), Some(0.0));
    sum.truncate();
    assert!(sum.is_empty());
}

#[test]
fn invalid_loss_targets_panic_before_mutation() {
    use std::panic::{AssertUnwindSafe, catch_unwind};
    let mut sum = LossyPauliSum::<f64>::new(1);
    sum += ("L", 1.0);
    let original = sum.clone();
    assert!(catch_unwind(AssertUnwindSafe(|| sum.loss_channel(1, 0.1))).is_err());
    assert!(
        catch_unwind(AssertUnwindSafe(
            || sum.correlated_loss_channel(0, 0, [0.1; 3])
        ))
        .is_err()
    );
    assert!(catch_unwind(AssertUnwindSafe(|| sum.reset_loss_channel(1))).is_err());
    assert_eq!(sum, original);
}
