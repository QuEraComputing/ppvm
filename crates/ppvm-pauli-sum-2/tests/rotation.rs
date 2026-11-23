// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use approx::assert_abs_diff_eq;
use num::Complex;
use ppvm_pauli_sum_2::{CoefficientThreshold, LossyPauliSum, LossyPauliWord, PauliSum, PauliWord};
use ppvm_traits_2::{Pauli, RotationOne};

#[test]
fn colliding_branches_use_original_coefficients() {
    for (axis, positive, negative) in [
        (Pauli::X, "Z", "Y"),
        (Pauli::Y, "X", "Z"),
        (Pauli::Z, "Y", "X"),
    ] {
        let mut sum = PauliSum::<f64>::new(1);
        sum += (positive, 2.0);
        sum += (negative, 3.0);
        sum += ("I", 4.0);
        sum.rotate_1(axis, 0, 0.37);
        let (sin, cos) = 0.37_f64.sin_cos();
        let mut expected = PauliSum::<f64>::new(1);
        expected += (positive, 2.0 * cos - 3.0 * sin);
        expected += (negative, 3.0 * cos + 2.0 * sin);
        expected += ("I", 4.0);
        assert_abs_diff_eq!(sum, expected, epsilon = 1e-12);
    }
}

#[test]
fn complex_coefficients_with_real_angles() {
    let mut sum = PauliSum::<Complex<f64>>::new(1);
    let coeff = Complex::new(2.0, -3.0);
    sum += ("Z", coeff);
    sum.rx(0, 0.37_f64);
    let mut expected = PauliSum::<Complex<f64>>::new(1);
    expected += ("Z", coeff * 0.37_f64.cos());
    expected += ("Y", coeff * 0.37_f64.sin());
    assert_abs_diff_eq!(sum, expected, epsilon = 1e-12);
}

#[test]
fn zeros_and_truncation_remain_explicit() {
    let mut sum = PauliSum::<f64, _>::with_policy(1, CoefficientThreshold { threshold: 10.0 });
    sum += ("Z", 1.0);
    sum.rx(0, 0.0);
    assert_eq!(sum.get(&PauliWord::from("Z")), Some(1.0));
    assert_eq!(sum.get(&PauliWord::from("Y")), Some(0.0));
    sum.reduce();
    assert_eq!(sum.len(), 1);
    sum.truncate();
    assert!(sum.is_empty());
}

#[test]
fn exact_cancellation_keeps_zero_until_reduced() {
    // Z cancels exactly: sin(theta) * cos(theta) - cos(theta) * sin(theta).
    let (sin, cos) = 0.37_f64.sin_cos();
    let mut sum = PauliSum::<f64>::new(1);
    sum += ("Y", cos);
    sum += ("Z", sin);
    sum.rx(0, 0.37);
    assert_eq!(sum.get(&PauliWord::from("Z")), Some(0.0));
    assert_eq!(sum.len(), 2);
    sum.reduce();
    assert_eq!(sum.len(), 1);
}

#[test]
fn lost_and_commuting_sites_are_unchanged() {
    for (axis, commuting) in [(Pauli::X, "X"), (Pauli::Y, "Y"), (Pauli::Z, "Z")] {
        let mut sum = LossyPauliSum::<f64>::new(1);
        sum += ("L", 2.0);
        sum += ("I", 3.0);
        sum += (commuting, 4.0);
        let original = sum.clone();
        sum.rotate_1(axis, 0, 0.37);
        assert_eq!(sum, original);
    }
    let mut sum = LossyPauliSum::<f64>::new(2);
    sum += ("LZ", 2.0);
    sum.rx(1, 0.37);
    assert_abs_diff_eq!(
        sum.get(&LossyPauliWord::from("LY")).unwrap(),
        2.0 * 0.37_f64.sin(),
        epsilon = 1e-12
    );
}

#[test]
fn identity_and_invalid_indices() {
    use std::panic::{AssertUnwindSafe, catch_unwind};
    let mut sum = PauliSum::<f64>::new(1);
    sum += ("Z", 1.0);
    let original = sum.clone();
    sum.rotate_1(Pauli::I, 0, 0.37);
    assert_eq!(sum, original);
    for axis in [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z] {
        assert!(catch_unwind(AssertUnwindSafe(|| sum.rotate_1(axis, 1, 0.37))).is_err());
        assert_eq!(sum, original);
    }
    let mut empty = PauliSum::<f64>::new(0);
    assert!(catch_unwind(AssertUnwindSafe(|| empty.rx(0, 0.37))).is_err());
}

#[test]
fn sequence_matches_existing_engine() {
    use old::RotationOne as _;
    use ppvm_pauli_sum::prelude as old;
    let mut legacy: old::PauliSum<old::config::fxhash::ByteF64<1>> =
        old::PauliSum::builder().n_qubits(2).build();
    let mut sum = PauliSum::<f64>::new(2);
    for (i, word) in ["II", "XY", "YZ", "ZX", "IZ", "IY"].iter().enumerate() {
        let coeff = (i as f64 - 2.0) / 3.0;
        legacy += (old::PauliWord::<[u8; 1]>::from(*word), coeff);
        sum += (*word, coeff);
    }
    for step in 0..24 {
        let theta = (step as f64 - 8.0) * 0.17;
        let q = step % 2;
        match step % 3 {
            0 => {
                legacy.rx(q, theta);
                sum.rx(q, theta);
            }
            1 => {
                legacy.ry(q, theta);
                sum.ry(q, theta);
            }
            _ => {
                legacy.rz(q, theta);
                sum.rz(q, theta);
            }
        }
        for (key, coeff) in legacy.iter() {
            let key = PauliWord::from(key.to_string().as_str());
            assert_abs_diff_eq!(sum.get(&key).unwrap(), *coeff, epsilon = 1e-12);
        }
        assert_eq!(sum.len(), legacy.len());
    }
}
