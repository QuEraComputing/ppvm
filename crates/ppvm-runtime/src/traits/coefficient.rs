// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

use num::traits::Float;

/// Numeric coefficient type usable inside a
/// [`PauliSum`](crate::sum::PauliSum).
///
/// `Coefficient` bundles every arithmetic operation a Pauli-propagation
/// step needs — addition, multiplication, signed multiplication, half,
/// sin/cos, and a `cutoff` predicate used by truncation strategies. The
/// built-in `f64` impl covers the common case; `Complex<f64>` is
/// available when phase tracking is required.
pub trait Coefficient:
    PartialEq
    + Clone
    + num::Zero
    + From<f64>
    + std::ops::Neg<Output = Self>
    + std::ops::Add<f64, Output = Self>
    + std::ops::Add<Self, Output = Self>
    + std::ops::Sub<f64, Output = Self>
    + std::ops::Sub<Self, Output = Self>
    + std::ops::Mul<f64, Output = Self>
    + std::ops::Mul<Self, Output = Self>
    + std::ops::AddAssign<f64>
    + std::ops::AddAssign<Self>
    + std::ops::MulAssign<f64>
    + std::ops::MulAssign<Self>
    + std::iter::Sum
    + Sync
    + Send
{
    /// Multiply by `sign ∈ {-1, +1}` (encoded as `i8`).
    fn mul_sign(&self, sign: i8) -> Self;
    /// Divide by two.
    fn half(&self) -> Self;
    /// Return `(sin θ, cos θ)`.
    fn sin_cos(&self) -> (Self, Self);

    /// Determine whether this coefficient should be cutoff
    /// Returns `true`, if the coefficient should be cut, and `false` else.
    fn cutoff(&self, threshold: f64) -> bool;
}

impl Coefficient for f64 {
    fn mul_sign(&self, sign: i8) -> Self {
        (sign as f64) * (*self)
    }

    fn half(&self) -> Self {
        *self / 2.0
    }

    fn sin_cos(&self) -> (Self, Self) {
        Float::sin_cos(*self)
    }

    fn cutoff(&self, threshold: f64) -> bool {
        self.abs() < threshold
    }
}

/// A [`Coefficient`] extended with multiplication by a fourth-root-of-unity
/// phase. Used by [`PhasedPauliWord`](crate::phase::PhasedPauliWord).
pub trait ComplexCoefficient: Coefficient {
    /// multiply by phase encoded as:
    ///
    /// |  | sign | imag |
    /// |--|------|------|
    /// |+1|    0 |    0 |
    /// |+i|    0 |    1 |
    /// |-1|    1 |    0 |
    /// |-i|    1 |    1 |
    fn mul_phase(&self, phase: u8) -> Self;
}

impl Coefficient for num::complex::Complex<f64> {
    fn cutoff(&self, threshold: f64) -> bool {
        self.norm() < threshold
    }

    fn half(&self) -> Self {
        *self / 2.0
    }

    fn mul_sign(&self, sign: i8) -> Self {
        (sign as f64) * (*self)
    }

    fn sin_cos(&self) -> (Self, Self) {
        let (s, c) = Float::sin_cos(self.re);
        (
            num::complex::Complex::new(s, 0.0),
            num::complex::Complex::new(c, 0.0),
        )
    }
}

impl ComplexCoefficient for num::complex::Complex<f64> {
    fn mul_phase(&self, phase: u8) -> Self {
        match phase % 4 {
            0 => *self,
            1 => num::complex::Complex::new(-self.im, self.re),
            2 => -*self,
            3 => num::complex::Complex::new(self.im, -self.re),
            _ => unreachable!(),
        }
    }
}
