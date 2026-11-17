// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::Complex;
use ppvm_pauli_sum_2::{CoefficientThreshold, LossyPauliSum, LossyPauliWord, PauliSum};
use ppvm_phased_pauli_word_2::Phased;
use ppvm_traits_2::{Clifford, CliffordBatch};
use std::panic::{AssertUnwindSafe, catch_unwind};

type Gate<T> = fn(&mut T, usize, usize);

// Single-qubit gates ignore the second index.
fn gates<T: Clifford>() -> [Gate<T>; 13] {
    [
        |s, q, _| s.x(q),
        |s, q, _| s.y(q),
        |s, q, _| s.z(q),
        |s, q, _| s.h(q),
        |s, q, _| s.s(q),
        |s, q, _| s.s_dag(q),
        |s, q, _| s.sqrt_x(q),
        |s, q, _| s.sqrt_x_dag(q),
        |s, q, _| s.sqrt_y(q),
        |s, q, _| s.sqrt_y_dag(q),
        T::cnot,
        T::cz,
        T::cy,
    ]
}

fn assert_same(actual: &PauliSum, expected: &PauliSum) {
    assert_eq!(actual.n_sites(), expected.n_sites());
    assert_eq!(actual.len(), expected.len());
    for (key, coeff) in expected.iter() {
        assert_eq!(actual.get(&key), Some(coeff));
    }
}

#[test]
fn propagation_matches_phased_words_including_loss() {
    let alphabet = ['I', 'X', 'Y', 'Z', 'L'];
    let mut original = LossyPauliSum::<f64>::new(2);
    for first in alphabet {
        for second in alphabet {
            let word = LossyPauliWord::from(format!("{first}{second}").as_str());
            original.add_term(word, original.len() as f64);
        }
    }
    for (gate_index, (apply_sum, apply_word)) in gates::<LossyPauliSum>()
        .into_iter()
        .zip(gates::<Phased<LossyPauliWord>>())
        .enumerate()
    {
        for (a, b) in [(0, 1), (1, 0)] {
            let mut sum = original.clone();
            apply_sum(&mut sum, a, b);
            assert_eq!(sum.len(), original.len());
            for (key, coeff) in original.iter() {
                let mut expected = Phased::new(key);
                apply_word(&mut expected, a, b);
                let coeff = if expected.is_positive() {
                    coeff
                } else {
                    -coeff
                };
                assert_eq!(
                    sum.get(expected.word()),
                    Some(coeff),
                    "gate {gate_index}, ({a}, {b})"
                );
            }
        }
    }
}

#[test]
fn complex_coefficients_absorb_signs_and_preserve_overlap() {
    let mut left = PauliSum::<Complex<f64>>::from_terms(
        2,
        [
            ("XY".into(), Complex::new(2.0, 3.0)),
            ("IZ".into(), Complex::new(-1.0, 4.0)),
        ],
    );
    let mut right = PauliSum::from_terms(2, [("XY".into(), Complex::new(1.0, -2.0))]);
    let overlap = left.hermitian_overlap(&right);
    left.z(0);
    assert_eq!(left.get(&"XY".into()), Some(Complex::new(-2.0, -3.0)));
    right.z(0);
    for sum in [&mut left, &mut right] {
        sum.h(1);
        sum.cnot(0, 1);
    }
    assert_eq!(left.get(&"YZ".into()), Some(Complex::new(2.0, 3.0)));
    assert_eq!(left.hermitian_overlap(&right), overlap);
}

#[test]
fn gates_keep_zero_terms_and_leave_truncation_explicit() {
    let mut sum = PauliSum::from_terms_with_policy(
        1,
        CoefficientThreshold { threshold: 1.0 },
        [("X".into(), 0.0), ("Z".into(), 0.25)],
    )
    .preserving(["X".into()]);
    sum.h(0);
    assert_eq!(sum.len(), 2);
    assert_eq!(sum.get(&"Z".into()), Some(0.0));
    assert_eq!(sum.get(&"X".into()), Some(0.25));
    assert!(sum.preserved_keys().contains(&"X".into()));
    sum.truncate();
    assert_eq!(sum.len(), 1);
    assert_eq!(sum.get(&"X".into()), Some(0.25));
}

#[test]
fn invalid_gate_arguments_do_not_change_terms() {
    for original in [
        PauliSum::new(2),
        PauliSum::from_terms(2, [("XY".into(), 2.0)]),
    ] {
        for gate in gates::<PauliSum>() {
            let mut sum = original.clone();
            assert!(catch_unwind(AssertUnwindSafe(|| gate(&mut sum, 2, 0))).is_err());
            assert_same(&sum, &original);
        }
        for gate in [PauliSum::cnot, PauliSum::cz, PauliSum::cy] {
            for (a, b) in [(0, 0), (0, 2)] {
                let mut sum = original.clone();
                assert!(catch_unwind(AssertUnwindSafe(|| gate(&mut sum, a, b))).is_err());
                assert_same(&sum, &original);
            }
        }
    }
}

#[test]
fn batches_and_aliases_follow_gate_order() {
    let mut actual = PauliSum::from_terms(
        3,
        [
            ("XYZ".into(), 2.0),
            ("IXZ".into(), -3.0),
            ("ZIX".into(), 0.0),
        ],
    );
    let mut expected = actual.clone();
    actual.s_many(&[0, 2, 0]);
    actual.cnot_many(&[(0, 1), (1, 2), (0, 1)]);
    for q in [0, 2, 0] {
        expected.s(q);
    }
    for (a, b) in [(0, 1), (1, 2), (0, 1)] {
        expected.cnot(a, b);
    }
    assert_same(&actual, &expected);

    let aliases: [(Gate<PauliSum>, Gate<PauliSum>); 4] = [
        (PauliSum::cx, PauliSum::cnot),
        (PauliSum::zcx, PauliSum::cnot),
        (PauliSum::zcz, PauliSum::cz),
        (PauliSum::zcy, PauliSum::cy),
    ];
    for (alias, gate) in aliases {
        alias(&mut actual, 0, 1);
        gate(&mut expected, 0, 1);
        assert_same(&actual, &expected);
    }
}
