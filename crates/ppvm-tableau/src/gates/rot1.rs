// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bitvec::view::BitView;
use num::PrimInt;
use num::{
    Complex, One, Zero,
    complex::{Complex64, ComplexFloat},
};

use crate::prelude::*;

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> RotationOne<T>
    for GeneralizedTableau<T, I, C>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    T::Coeff: Zero + One + Send + Sync + num::Num + PartialOrd,
    I: TableauIndex + Send + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
{
    fn rotate_1(&mut self, axis: Pauli, addr0: usize, theta: <T as Config>::Coeff) {
        if self.is_lost_or_leaked(addr0) {
            return;
        }
        let (sin, cos) = (theta * 0.5.into()).sin_cos();

        let complex_cos: Complex<T::Coeff> = Complex {
            re: cos,
            im: T::Coeff::zero(),
        };

        let i_complex_sin: Complex<T::Coeff> = Complex {
            re: T::Coeff::zero(),
            im: -sin,
        };

        self.branch_with_coefficients(addr0, axis, complex_cos, i_complex_sin);
    }
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> RotXY<T> for GeneralizedTableau<T, I, C>
where
    GeneralizedTableau<T, I, C>: RotationOne<T>,
{
    fn r(&mut self, addr0: usize, axis_angle: T::Coeff, theta: T::Coeff) {
        // R(axis_angle, θ) = RZ(axis_angle)·RX(θ)·RZ(−axis_angle). The tableau
        // runs in the Schrödinger picture, so the sub-rotations are applied in
        // forward order: RZ(−axis_angle) first, then RX(θ), then RZ(axis_angle).
        self.rz(addr0, -axis_angle.clone());
        self.rx(addr0, theta);
        self.rz(addr0, axis_angle);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_pauli_sum::config::fxhash::ByteF64;
    use std::f64::consts::{FRAC_PI_2, PI};

    type TestConfig = ByteF64<1>;
    type TestTableau = GeneralizedTableau<TestConfig>;

    /// rx(π) = exp(-i·π/2·X) acts as -iX on |0⟩, giving -i|1⟩.
    #[test]
    fn test_rx_pi_flips_zero_to_one() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.rx(0, PI);

        // make sure we don't branch
        assert_eq!(tab.coefficients.len(), 1);

        let result = tab.measure(0);
        assert!(result.unwrap(), "Expected to measure 1, got {:?}", result);
    }

    #[test]
    fn test_ry_pi_flips_zero_to_one() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.ry(0, PI);
        // make sure we don't branch
        assert_eq!(tab.coefficients.len(), 1);
        let result = tab.measure(0);
        assert!(result.unwrap(), "Expected to measure 1, got {:?}", result);
    }

    #[test]
    fn test_rz_leaves_zero_invariant() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.rz(0, 0.123);
        // make sure we don't branch
        assert_eq!(tab.coefficients.len(), 1);
        let result = tab.measure(0);
        assert!(!result.unwrap(), "Expected to measure 0, got {:?}", result);
    }

    #[test]
    fn test_rx_pi_2_statistics() {
        let mut tab: TestTableau = GeneralizedTableau::new_with_seed(1, 1e-12, 0);
        tab.rx(0, FRAC_PI_2);
        // make sure we branch
        assert_eq!(tab.coefficients.len(), 2);

        let trials = 100;
        let mut count_one = 0;
        for i in 0..trials {
            let mut tmp_tab = tab.fork(Some(i));
            let result = tmp_tab.measure(0);
            count_one += result.unwrap() as u8;
        }
        println!("{}", count_one);
        assert!(35 < count_one && count_one < 65);
    }

    #[test]
    fn test_ry_round_trip() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-10);
        tab.ry(0, 2.0 * PI);

        assert_eq!(tab.coefficients.len(), 1);

        assert!(!tab.measure(0).unwrap());
    }

    /// R(axis_angle=0, θ=0) = identity: |0⟩ stays |0⟩, no branching.
    #[test]
    fn test_r_identity() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.r(0, 0.0, 0.0);
        assert_eq!(tab.coefficients.len(), 1);
        assert!(!tab.measure(0).unwrap());
    }

    /// R(axis_angle=0, θ=π) = RX(π): flips |0⟩ → |1⟩ with no branching.
    #[test]
    fn test_r_axis_zero_is_rx() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.r(0, 0.0, PI);
        assert_eq!(tab.coefficients.len(), 1);
        assert!(tab.measure(0).unwrap());
    }

    /// R(axis_angle=π/2, θ=π) = RY(π): flips |0⟩ → |1⟩ with no branching.
    #[test]
    fn test_r_axis_half_pi_is_ry() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.r(0, FRAC_PI_2, PI);
        assert_eq!(tab.coefficients.len(), 1);
        assert!(tab.measure(0).unwrap());
    }

    /// A partial rotation about an in-plane axis branches into two terms.
    #[test]
    fn test_r_branches() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.r(0, 0.37 * PI, FRAC_PI_2);
        assert_eq!(tab.coefficients.len(), 2);
    }

    /// R(axis_angle, θ) must give identical per-seed measurement statistics
    /// to the manual decomposition RZ(axis_angle)·RX(θ)·RZ(−axis_angle).
    #[test]
    fn test_r_matches_rz_rx_rz() {
        let (axis_angle, theta) = (0.21 * PI, 0.34 * PI);

        let mut tab_r: TestTableau = GeneralizedTableau::new_with_seed(1, 1e-12, 0);
        tab_r.r(0, axis_angle, theta);

        let mut tab_manual: TestTableau = GeneralizedTableau::new_with_seed(1, 1e-12, 0);
        tab_manual.rz(0, -axis_angle);
        tab_manual.rx(0, theta);
        tab_manual.rz(0, axis_angle);

        for seed in 0..200 {
            let result_r = tab_r.fork(Some(seed)).measure(0).unwrap();
            let result_manual = tab_manual.fork(Some(seed)).measure(0).unwrap();
            assert_eq!(result_r, result_manual, "mismatch at seed {}", seed);
        }
    }

    #[test]
    fn test_two_qubit_case() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-10);
        tab.rx(0, PI);

        assert_eq!(tab.coefficients.len(), 1);

        assert!(!tab.measure(1).unwrap());
        assert!(tab.measure(0).unwrap());

        tab.rx(1, PI);

        assert_eq!(tab.coefficients.len(), 1);

        assert!(tab.measure(1).unwrap());
        assert!(tab.measure(0).unwrap());

        tab.rx(0, PI);

        assert_eq!(tab.coefficients.len(), 1);

        assert!(tab.measure(1).unwrap());
        assert!(!tab.measure(0).unwrap());

        tab.rx(1, PI);

        assert_eq!(tab.coefficients.len(), 1);

        assert!(!tab.measure(1).unwrap());
        assert!(!tab.measure(0).unwrap());
    }
}
