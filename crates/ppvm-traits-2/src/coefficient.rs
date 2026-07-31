// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The coefficient value domain and the rotation-angle domain.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Coefficient, angle, and
//! truncation". `Coefficient` is the old `ppvm_traits::Coefficient` bundle with
//! its non-value-domain responsibilities extracted: `sin_cos` moves to
//! [`Angle`], `cutoff` is replaced by [`Coefficient::magnitude`] (a property a
//! `Policy` thresholds), and the partial `0.5·x` operation moves to the
//! [`Halvable`] capability. Crucially it also drops the old `Mul<f64>` bound —
//! the one bound that excluded exact rings — so `GaussianInt` (`ℤ[i]`),
//! `Complex<Rational>`, and cyclotomic integers can be `Coefficient`s without a
//! lossy division.

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub};

/// Value-domain ring arithmetic only — no rotation, no truncation, no
/// bare-`f64` scaling (so exact rings like `GaussianInt` qualify).
///
/// Design: §"Coefficient, angle, and truncation". This is the value domain the
/// graded algebra of `docs/design` §"The map is a graded algebra over `C[K]`"
/// builds `C[K]` over; the module laws that consume it are machine-checked in
/// `lean/PPVM/Algebra/GradedMap.lean` (`accumulate_comm`, `accumulate_assoc`,
/// `scale_scale`).
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

    /// Nonnegative magnitude. Exposes a property of the value for a `Policy` to
    /// threshold; it does not itself decide any cutoff. Replaces the old
    /// `Coefficient::cutoff`.
    ///
    /// Design: §"Coefficient, angle, and truncation". The threshold comparison
    /// that consumes this lives in `Policy` (`ppvm-pauli-sum-2`), whose keep-rule
    /// boundary against the tableau path is machine-checked in
    /// `lean/PPVM/Algebra/Truncation.lean` (`cutoff_mismatch`).
    fn magnitude(&self) -> f64;
}

/// A coefficient ring in which halving (`0.5·x`) is total and exact: the
/// capability the projective computational-basis measurement kernel needs to
/// apply the `(I ± Z)/2` projectors.
///
/// This is split out of [`Coefficient`] for the same reason `sin_cos` became
/// [`Angle`] and `i`/conjugation became [`crate::algebra::ImaginaryUnit`] /
/// [`crate::algebra::Conjugate`]: it is **not** value-domain ring arithmetic and
/// **not** required by propagation (Clifford, rotation, and the twisted product
/// never halve). Keeping it a `Coefficient` bound is exactly the `Mul<f64>`
/// mistake in another guise — `0.5·(1+i)` leaves `ℤ[i]`, so an exact ring could
/// only satisfy it with a lossy integer `/2` for which `half(x) + half(x) != x`.
/// Making it a separate capability lets exact rings be `Coefficient`s (the L4
/// witness `GaussianInt` deliberately does **not** implement `Halvable`) while
/// the float and `Complex<f64>` domains that back Phase-1 measurement do.
///
/// Design: §"Coefficient, angle, and truncation".
pub trait Halvable: Coefficient {
    /// Divide by two. Impls must be exact: `x.half() + x.half() == x`.
    fn half(&self) -> Self;
}

/// A rotation angle that yields `(sin, cos)` already in coefficient domain `C`.
///
/// Design: §"Coefficient, angle, and truncation". Defaulting the angle to the
/// coefficient (see [`crate::gates::RotationOne`]) recovers today's `rx(theta: C)`
/// while permitting a symbolic/parametric angle over an `f64`-coefficient sum.
/// The 2-D rotation this drives on a branching term pair is norm-preserving and
/// angle-additive — machine-checked in `lean/PPVM/Instantiations/Rotation.lean`
/// (`rot_norm_sq`, `rot_rot`).
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

impl Coefficient for f64 {
    #[inline]
    fn mul_sign(&self, sign: i8) -> Self {
        (sign as f64) * (*self)
    }

    #[inline]
    fn magnitude(&self) -> f64 {
        self.abs()
    }
}

impl Halvable for f64 {
    #[inline]
    fn half(&self) -> Self {
        *self / 2.0
    }
}

impl Coefficient for num::Complex<f64> {
    #[inline]
    fn mul_sign(&self, sign: i8) -> Self {
        (sign as f64) * (*self)
    }

    #[inline]
    fn magnitude(&self) -> f64 {
        self.norm()
    }
}

impl Halvable for num::Complex<f64> {
    #[inline]
    fn half(&self) -> Self {
        *self / 2.0
    }
}
