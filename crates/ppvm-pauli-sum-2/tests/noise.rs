// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use approx::assert_abs_diff_eq;
use num::Complex;
use ppvm_pauli_sum_2::{CoefficientThreshold, LossyPauliSum, LossyPauliWord, PauliSum, PauliWord};
use ppvm_traits_2::{Depolarizing, Depolarizing2, PauliError, TwoQubitPauliError};

#[test]
fn pauli_channel_scales_complex_coefficients_and_preserves_loss() {
    let mut sum = LossyPauliSum::<Complex<f64>>::new(1);
    let c = Complex::new(2.0, -3.0);
    for word in ["I", "X", "Y", "Z", "L"] {
        sum += (word, c);
    }
    sum.pauli_error(0, [0.1.into(), 0.2.into(), 0.3.into()]);
    for (word, factor) in [("I", 1.0), ("X", 0.0), ("Y", 0.2), ("Z", 0.4), ("L", 1.0)] {
        assert_abs_diff_eq!(
            sum.get(&LossyPauliWord::from(word)).unwrap(),
            c * factor,
            epsilon = 1e-12
        );
    }
    assert_eq!(sum.len(), 5);
}

#[test]
fn all_two_qubit_pauli_errors_have_the_declared_order() {
    let mut original = LossyPauliSum::<f64>::new(2);
    for a in "IXYZL".chars() {
        for b in "IXYZL".chars() {
            original += (format!("{a}{b}").as_str(), 1.0);
        }
    }
    let errors: Vec<_> = "IXYZ"
        .chars()
        .flat_map(|a| "IXYZ".chars().map(move |b| (a, b)))
        .skip(1)
        .collect();
    for (i, (a, b)) in errors.into_iter().enumerate() {
        let mut sum = original.clone();
        let mut p = [0.0; 15];
        p[i] = 1.0;
        sum.two_qubit_pauli_error(0, 1, p);
        for x in "IXYZL".chars() {
            for y in "IXYZL".chars() {
                let anti =
                    |error, word| error != 'I' && word != 'I' && word != 'L' && error != word;
                let sign = if anti(a, x) ^ anti(b, y) { -1.0 } else { 1.0 };
                assert_eq!(
                    sum.get(&LossyPauliWord::from(format!("{x}{y}").as_str())),
                    Some(sign)
                );
            }
        }
    }
}

#[test]
fn depolarization_matches_uniform_pauli_channels_including_loss() {
    let mut original = LossyPauliSum::<f64>::new(2);
    for a in "IXYZL".chars() {
        for b in "IXYZL".chars() {
            original += (format!("{a}{b}").as_str(), 1.0);
        }
    }
    for p in [0.0, 0.3, 1.0] {
        let mut sum = original.clone();
        let mut expected = original.clone();
        sum.depolarize1(0, p);
        expected.pauli_error(0, [p / 3.0; 3]);
        assert_abs_diff_eq!(sum, expected, epsilon = 1e-12);
        let mut sum = original.clone();
        let mut expected = original.clone();
        sum.depolarize2(0, 1, p);
        expected.two_qubit_pauli_error(0, 1, [p / 15.0; 15]);
        assert_abs_diff_eq!(sum, expected, epsilon = 1e-12);
    }
}

#[test]
fn noise_does_not_reduce_or_truncate() {
    let mut sum = PauliSum::<f64, _>::with_policy(1, CoefficientThreshold { threshold: 10.0 });
    sum += ("X", 1.0);
    sum.x_error(0, 0.5);
    assert_eq!(sum.get(&PauliWord::from("X")), Some(1.0));
    sum.z_error(0, 0.5);
    assert_eq!(sum.get(&PauliWord::from("X")), Some(0.0));
    sum.depolarize1(0, 0.75);
    assert_eq!(sum.len(), 1);
    sum.reduce();
    assert!(sum.is_empty());
    sum += ("Z", 1.0);
    sum.depolarize1(0, 0.3);
    assert_eq!(sum.len(), 1);
    sum.truncate();
    assert!(sum.is_empty());
}

#[test]
fn invalid_noise_targets_panic_before_mutation() {
    use std::panic::{AssertUnwindSafe, catch_unwind};
    let mut sum = PauliSum::<f64>::new(1);
    sum += ("Z", 1.0);
    let original = sum.clone();
    assert!(catch_unwind(AssertUnwindSafe(|| sum.pauli_error(1, [0.1; 3]))).is_err());
    assert!(
        catch_unwind(AssertUnwindSafe(
            || sum.two_qubit_pauli_error(0, 0, [0.01; 15])
        ))
        .is_err()
    );
    assert!(catch_unwind(AssertUnwindSafe(|| sum.depolarize1(1, 0.1))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| sum.depolarize2(0, 1, 0.1))).is_err());
    assert_eq!(sum, original);
    let mut empty = PauliSum::<f64>::new(0);
    assert!(catch_unwind(AssertUnwindSafe(|| empty.depolarize1(0, 0.1))).is_err());
}
