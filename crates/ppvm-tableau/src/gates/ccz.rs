// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bitvec::view::BitView;
use fxhash::FxHashMap;
use num::{
    Complex, One, PrimInt, Zero,
    complex::{Complex64, ComplexFloat},
};

use crate::prelude::*;

impl<T, I, C> CCZGate<T> for GeneralizedTableau<T, I, C>
where
    T: Config,
    <T::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I>,
    T::Coeff: One + Zero + Clone + Send + Sync + num::Num + PartialOrd,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
    I: TableauIndex + Send + Sync,
{
    fn ccz(&mut self, a: usize, b: usize, c: usize) {
        assert!(a != b && a != c && b != c, "CCZ requires distinct qubits");
        assert!(a < self.n_qubits() && b < self.n_qubits() && c < self.n_qubits());
        if self.any_inactive(&[a, b, c]) {
            return;
        }

        // CCZ = (3I + Za + Zb + Zc - ZaZb - ZaZc - ZbZc + ZaZbZc) / 4.
        let mut merged = FxHashMap::<I, Complex<T::Coeff>>::default();
        for (mask, weight) in [0.75, 0.25, 0.25, -0.25, 0.25, -0.25, -0.25, 0.25]
            .into_iter()
            .enumerate()
        {
            let mut branch = self.coefficients.clone();
            for (bit, qubit) in [a, b, c].into_iter().enumerate() {
                if mask & (1 << bit) != 0 {
                    self.compute_coefficients_after_pauli_apply_untruncated(
                        &mut branch,
                        qubit,
                        Pauli::Z,
                    );
                }
            }
            let factor: Complex<T::Coeff> = Complex64::new(weight, 0.0).into();
            for (coeff, idx) in branch {
                *merged.entry(idx).or_insert(Complex::zero()) += factor * coeff;
            }
        }

        // Coherent branches can cancel or reinforce: truncate only their sum.
        let cutoff_sq = self.coefficient_threshold.clone() * self.coefficient_threshold.clone();
        self.coefficients = C::new();
        self.coefficients.reserve(merged.len());
        for (idx, coeff) in merged {
            if coeff.norm_sqr() > cutoff_sq {
                self.coefficients.unsafe_insert(idx, coeff);
            }
        }
    }
}
