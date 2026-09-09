// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::arithmetic::Coefficient;

// A rotation angle that yields `(sin, cos)` already in coefficient domain `C`.
pub trait Angle<C: Coefficient> {
    /// Return `(sin θ, cos θ)` in the coefficient domain `C`.
    fn sin_cos(&self) -> (C, C);
}

impl Angle<f64> for f64 {
    #[inline]
    fn sin_cos(&self) -> (f64, f64) {
        num::traits::Float::sin_cos(*self)
    }
}

/// The complex-coefficient angle domain, i.e. the defaulted `A = C` case of
/// [`crate::gates::RotationOne`] at `C = Complex<f64>`.
impl Angle<num::Complex<f64>> for num::Complex<f64> {
    #[inline]
    fn sin_cos(&self) -> (num::Complex<f64>, num::Complex<f64>) {
        let (s, c) = num::traits::Float::sin_cos(self.re);
        (num::Complex::new(s, 0.0), num::Complex::new(c, 0.0))
    }
}

// A **real** angle driving a complex-coefficient sum.
impl Angle<num::Complex<f64>> for f64 {
    #[inline]
    fn sin_cos(&self) -> (num::Complex<f64>, num::Complex<f64>) {
        let (s, c) = num::traits::Float::sin_cos(*self);
        (num::Complex::new(s, 0.0), num::Complex::new(c, 0.0))
    }
}
