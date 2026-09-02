// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The coefficient scalar of a Pauli-sum basis: real (`f64`) on the
//! plain adaptive path, complex on the momentum-sector orbit-rep path.
//!
//! Every truncation and 1-norm decision in this crate only ever needs a
//! magnitude and a zero, so the two paths share one implementation
//! parameterised by [`Coeff`] instead of a real and a complex copy.

use num::Complex;

/// A Pauli-sum coefficient: `f64` or `Complex<f64>`.
pub(crate) trait Coeff: Copy + Send + Sync {
    /// Absolute value (`f64::abs`) / modulus (`Complex::norm`).
    fn mag(self) -> f64;

    fn zero() -> Self;
}

impl Coeff for f64 {
    #[inline]
    fn mag(self) -> f64 {
        self.abs()
    }

    #[inline]
    fn zero() -> Self {
        0.0
    }
}

impl Coeff for Complex<f64> {
    #[inline]
    fn mag(self) -> f64 {
        self.norm()
    }

    #[inline]
    fn zero() -> Self {
        Complex::new(0.0, 0.0)
    }
}
