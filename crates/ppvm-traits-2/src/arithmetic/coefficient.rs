// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub};

pub trait Coefficient:
    PartialEq
    + Clone
    + num::Zero
    + Neg<Output = Self>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + AddAssign<Self>
    + MulAssign<Self>
    + std::iter::Sum
    + Send
    + Sync
{
    /// Multiply by `sign ∈ {-1, +1}` (encoded as `i8`).
    fn mul_sign(&self, sign: i8) -> Self;

    /// Multiply this coefficient in place by `sign ∈ {-1, +1}`.
    #[inline]
    fn mul_sign_assign(&mut self, sign: i8) {
        *self = self.mul_sign(sign)
    }

    ///  Accumulates a borrowed coefficient.
    #[inline]
    fn add_assign_ref(&mut self, rhs: &Self) {
        *self += rhs.clone();
    }

    /// Add this coefficient to itself. Numeric implementations may use their
    /// native multiply-by-two operation; exact rings retain the additive default.
    #[inline(always)]
    fn doubled(&self) -> Self {
        self.clone() + self.clone()
    }

    /// Nonnegative magnitude. Exposes a property of the value for a `Policy` to
    /// threshold; it does not itself decide any cutoff. Replaces the old
    /// `Coefficient::cutoff`.
    fn magnitude(&self) -> f64;
}

impl Coefficient for f64 {
    #[inline]
    fn mul_sign(&self, sign: i8) -> Self {
        (sign as f64) * (*self)
    }

    #[inline(always)]
    fn doubled(&self) -> Self {
        *self * 2.0
    }

    #[inline]
    fn magnitude(&self) -> f64 {
        self.abs()
    }
}

impl Coefficient for num::Complex<f64> {
    #[inline]
    fn mul_sign(&self, sign: i8) -> Self {
        (sign as f64) * (*self)
    }

    #[inline(always)]
    fn doubled(&self) -> Self {
        *self * 2.0
    }

    #[inline]
    fn magnitude(&self) -> f64 {
        self.norm()
    }
}
