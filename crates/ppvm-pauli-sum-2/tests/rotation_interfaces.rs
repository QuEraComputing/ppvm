// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use approx::assert_abs_diff_eq;
use num::Complex;
use ppvm_pauli_sum_2::{LossyPauliSum, LossyPauliWord, PauliSum};
use ppvm_traits_2::{
    CRx, Pauli, RotXY, RotationOne, RotationOneBatch, RotationTwo, RotationTwoBatch, TGate, U3Gate,
};

const AXES: [(Pauli, char); 4] = [
    (Pauli::I, 'I'),
    (Pauli::X, 'X'),
    (Pauli::Y, 'Y'),
    (Pauli::Z, 'Z'),
];
type C = Complex<f64>;
type Matrix = [[C; 4]; 4];

fn multiply(a: Matrix, b: Matrix) -> Matrix {
    std::array::from_fn(|i| std::array::from_fn(|j| (0..4).map(|k| a[i][k] * b[k][j]).sum()))
}

fn pauli(a: Pauli, b: Pauli) -> Matrix {
    fn entry(p: Pauli, i: usize, j: usize) -> C {
        match p {
            Pauli::I => C::new((i == j) as u8 as f64, 0.0),
            Pauli::X => C::new((i != j) as u8 as f64, 0.0),
            Pauli::Y if i != j => C::new(0.0, if i == 0 { -1.0 } else { 1.0 }),
            Pauli::Z if i == j => C::new(if i == 0 { 1.0 } else { -1.0 }, 0.0),
            _ => C::new(0.0, 0.0),
        }
    }
    std::array::from_fn(|i| {
        std::array::from_fn(|j| entry(a, i / 2, j / 2) * entry(b, i % 2, j % 2))
    })
}

fn rotation(a: Pauli, b: Pauli, theta: f64) -> Matrix {
    let generator = pauli(a, b);
    let identity = pauli(Pauli::I, Pauli::I);
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            identity[i][j] * (theta / 2.0).cos() - C::i() * generator[i][j] * (theta / 2.0).sin()
        })
    })
}

fn observable() -> PauliSum<C> {
    let mut sum = PauliSum::new(2);
    for (i, (_, a)) in AXES.iter().enumerate() {
        for (j, (_, b)) in AXES.iter().enumerate() {
            sum += (
                format!("{a}{b}").as_str(),
                C::new((i * 4 + j + 1) as f64 / 17.0, (i as f64 - j as f64) / 9.0),
            );
        }
    }
    sum
}

fn matrix(sum: &PauliSum<C>) -> Matrix {
    let mut out = [[C::new(0.0, 0.0); 4]; 4];
    for (a, ca) in AXES {
        for (b, cb) in AXES {
            let key = ppvm_pauli_sum_2::PauliWord::from(format!("{ca}{cb}").as_str());
            let coeff = sum.get(&key).unwrap_or_default();
            let term = pauli(a, b);
            for i in 0..4 {
                for j in 0..4 {
                    out[i][j] += coeff * term[i][j];
                }
            }
        }
    }
    out
}

fn check_conjugation(sum: &PauliSum<C>, unitary: Matrix) {
    let adjoint = std::array::from_fn(|i| std::array::from_fn(|j| unitary[j][i].conj()));
    let expected = multiply(multiply(adjoint, matrix(&observable())), unitary);
    let actual = matrix(sum);
    for i in 0..4 {
        for j in 0..4 {
            assert_abs_diff_eq!(actual[i][j], expected[i][j], epsilon = 1e-12);
        }
    }
}

#[test]
fn all_two_qubit_axes_match_dense_conjugation() {
    for (a, _) in AXES {
        for (b, _) in AXES {
            for theta in [0.0, 0.37, -0.81] {
                let mut sum = observable();
                sum.rotate_2(a, b, 0, 1, theta);
                check_conjugation(&sum, rotation(a, b, theta));
            }
        }
    }
}

#[test]
fn lost_sites_reduce_to_surviving_rotation() {
    for (a, _) in AXES {
        for (b, _) in AXES {
            for word in ["LX", "LY", "LZ", "XL", "YL", "ZL", "LL"] {
                let mut sum = LossyPauliSum::<f64>::new(2);
                sum += (word, 2.0);
                let mut expected = sum.clone();
                if word.starts_with('L') {
                    expected.rotate_1(b, 1, 0.37);
                } else {
                    expected.rotate_1(a, 0, 0.37);
                }
                sum.rotate_2(a, b, 0, 1, 0.37);
                assert_abs_diff_eq!(sum, expected, epsilon = 1e-12);
            }
        }
    }
}

#[test]
fn composite_interfaces_match_dense_unitaries() {
    let mut sum = observable();
    sum.u3(1, 0.37, -0.51, 0.82);
    check_conjugation(
        &sum,
        multiply(
            multiply(
                rotation(Pauli::I, Pauli::Z, -0.51),
                rotation(Pauli::I, Pauli::Y, 0.37),
            ),
            rotation(Pauli::I, Pauli::Z, 0.82),
        ),
    );

    let mut sum = observable();
    sum.r(1, 0.51, 0.37);
    check_conjugation(
        &sum,
        multiply(
            multiply(
                rotation(Pauli::I, Pauli::Z, 0.51),
                rotation(Pauli::I, Pauli::X, 0.37),
            ),
            rotation(Pauli::I, Pauli::Z, -0.51),
        ),
    );

    let mut sum = observable();
    sum.t(1);
    check_conjugation(
        &sum,
        rotation(Pauli::I, Pauli::Z, std::f64::consts::FRAC_PI_4),
    );
    sum.t_dag(1);
    assert_abs_diff_eq!(sum, observable(), epsilon = 1e-12);

    for theta in [0.0, 0.37, -0.81] {
        let mut sum = observable();
        sum.crx(0, 1, theta);
        // Build a controlled Rx directly: identity in the control=0 block.
        let rx = rotation(Pauli::I, Pauli::X, theta);
        let mut unitary = pauli(Pauli::I, Pauli::I);
        for i in 2..4 {
            for j in 2..4 {
                unitary[i][j] = rx[i][j];
            }
        }
        check_conjugation(&sum, unitary);
    }
}

#[test]
fn batches_match_ordered_scalar_calls() {
    type Single = fn(&mut PauliSum<C>, usize, f64);
    type Singles = fn(&mut PauliSum<C>, &[usize], f64);
    let singles: [(Single, Singles); 3] = [
        (PauliSum::rx, PauliSum::rx_many),
        (PauliSum::ry, PauliSum::ry_many),
        (PauliSum::rz, PauliSum::rz_many),
    ];
    for (scalar, batch) in singles {
        let mut sum = observable();
        let mut expected = sum.clone();
        batch(&mut sum, &[1, 0, 1], 0.37);
        for q in [1, 0, 1] {
            scalar(&mut expected, q, 0.37);
        }
        assert_abs_diff_eq!(sum, expected, epsilon = 1e-12);
    }
    type Pair = fn(&mut PauliSum<C>, usize, usize, f64);
    type Pairs = fn(&mut PauliSum<C>, &[(usize, usize)], f64);
    let pairs: [(Pair, Pairs); 9] = [
        (PauliSum::rxx, PauliSum::rxx_many),
        (PauliSum::rxy, PauliSum::rxy_many),
        (PauliSum::rxz, PauliSum::rxz_many),
        (PauliSum::ryx, PauliSum::ryx_many),
        (PauliSum::ryy, PauliSum::ryy_many),
        (PauliSum::ryz, PauliSum::ryz_many),
        (PauliSum::rzx, PauliSum::rzx_many),
        (PauliSum::rzy, PauliSum::rzy_many),
        (PauliSum::rzz, PauliSum::rzz_many),
    ];
    for (scalar, batch) in pairs {
        let mut sum = observable();
        let mut expected = sum.clone();
        batch(&mut sum, &[(0, 1), (1, 0)], 0.37);
        scalar(&mut expected, 0, 1, 0.37);
        scalar(&mut expected, 1, 0, 0.37);
        assert_abs_diff_eq!(sum, expected, epsilon = 1e-12);
    }
    let mut sum = observable();
    sum.t_many(&[0, 1]);
    sum.t_dag_many(&[0, 1]);
    assert_abs_diff_eq!(sum, observable(), epsilon = 1e-12);
}

#[test]
fn invalid_pairs_panic_before_mutation() {
    use std::panic::{AssertUnwindSafe, catch_unwind};
    for (a, b) in [(0, 0), (0, 2), (2, 0)] {
        let mut sum = observable();
        assert!(catch_unwind(AssertUnwindSafe(|| sum.rxx(a, b, 0.37))).is_err());
        assert_eq!(sum, observable());
        assert!(catch_unwind(AssertUnwindSafe(|| sum.crx(a, b, 0.37))).is_err());
        assert_eq!(sum, observable());
    }
    let mut empty = PauliSum::<f64>::new(0);
    assert!(catch_unwind(AssertUnwindSafe(|| empty.rzz(0, 1, 0.37))).is_err());
}

#[test]
fn two_qubit_zero_branches_remain_until_reduced() {
    let mut sum = PauliSum::<f64>::new(2);
    sum += ("ZI", 2.0);
    sum.rxx(0, 1, 0.0);
    assert_eq!(sum.len(), 2);
    assert_eq!(sum.get(&ppvm_pauli_sum_2::PauliWord::from("YX")), Some(0.0));
    sum.reduce();
    assert_eq!(sum.len(), 1);

    let mut sum = observable();
    sum.crx(0, 1, C::new(0.37, 0.0));
    let mut expected = observable();
    expected.crx(0, 1, 0.37);
    assert_abs_diff_eq!(sum, expected, epsilon = 1e-12);
}

#[test]
fn controlled_rotation_loss_follows_generator_decomposition() {
    let mut sum = LossyPauliSum::<f64>::new(2);
    sum += ("LZ", 1.0);
    sum.crx(0, 1, 0.37);
    assert_abs_diff_eq!(
        sum.get(&LossyPauliWord::from("LZ")).unwrap(),
        1.0,
        epsilon = 1e-12
    );
    assert_abs_diff_eq!(
        sum.get(&LossyPauliWord::from("LY")).unwrap(),
        0.0,
        epsilon = 1e-12
    );
    let mut sum = LossyPauliSum::<f64>::new(2);
    sum += ("XL", 1.0);
    let mut expected = sum.clone();
    expected.rz(0, -0.37 / 2.0);
    sum.crx(0, 1, 0.37);
    assert_abs_diff_eq!(sum, expected, epsilon = 1e-12);
}

#[test]
fn specialized_pairs_match_generic_with_loss_and_complex_coefficients() {
    type Gate = fn(&mut LossyPauliSum<C>, usize, usize, f64);
    let gates: [(Pauli, Gate); 3] = [
        (Pauli::X, LossyPauliSum::rxx),
        (Pauli::Y, LossyPauliSum::ryy),
        (Pauli::Z, LossyPauliSum::rzz),
    ];
    for (axis, gate) in gates {
        for theta in [0.0, 0.37, -0.81] {
            for (a, b) in [(0, 1), (1, 0)] {
                let mut sum = LossyPauliSum::<C>::new(2);
                for (i, x) in "IXYZL".chars().enumerate() {
                    for (j, y) in "IXYZL".chars().enumerate() {
                        sum += (
                            format!("{x}{y}").as_str(),
                            C::new(i as f64 - 2.0, j as f64 - 1.0),
                        );
                    }
                }
                let mut expected = sum.clone();
                expected.rotate_2(axis, axis, a, b, theta);
                gate(&mut sum, a, b, theta);
                assert_abs_diff_eq!(sum, expected, epsilon = 1e-12);
            }
        }
    }
}
