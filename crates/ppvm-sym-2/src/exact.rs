// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`GaussianInt`] — the Gaussian integers `ℤ[i]`, the **exact** coefficient
//! ring that closes the L4 loop.
//!
//! # Why this type exists
//!
//! The `-2` design split `Halvable`, `Angle<C>`, `ImaginaryUnit` and `Conjugate`
//! off `Coefficient` *specifically* so that exact rings stay expressible: gap
//! `t2.coefficient.1` records that keeping `Coefficient::half` would "foreclose
//! exact rings (`0.5·(1+i) ∉ ℤ[i]`)", and dropping the old `Mul<f64>` bound was
//! the other half of the same argument. That is a *claim about the trait
//! tower*, and a claim is only worth what its witness is worth — so this type
//! validates it by construction:
//!
//! * its representation is **two `i64`s**; there is no `f64` anywhere in it;
//! * it implements [`Coefficient`], [`ImaginaryUnit`] and [`Conjugate`], which
//!   is exactly the capability set `ppvm-pauli-sum-2`'s
//!   [`Multiply`](ppvm_traits_2::Multiply) (L4) requires; and
//! * it deliberately does **not** implement
//!   [`Halvable`](ppvm_traits_2::Halvable) — `half(1 + i)` has no `ℤ[i]` value
//!   satisfying `x.half() + x.half() == x`, which is the whole point.
//!
//! If any bound in the tower still secretly forced a float, this type would not
//! compile; `tests/exact_ring.rs` drives `Sum<HashMapStore<PauliWord<…>,
//! GaussianInt>, …>` through `multiply_into` to show the L4 product genuinely
//! runs over it.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"The map is a graded algebra
//! over `C[K]`"; implementation plan §Phase 5. Lean: `lean/PPVM/Pauli/Matrix.lean`
//! (`iU_sq`, `star_iU`) and `lean/PPVM/Algebra/Twisted.lean` (`tmul_assoc`,
//! `twistedConv_assoc`, `iPow_add`) — all stated over an arbitrary commutative
//! ring with `i⁴ = 1`, which is precisely what `ℤ[i]` is.

use ppvm_traits_2::{Coefficient, Conjugate, ImaginaryUnit};

/// A Gaussian integer `re + im·i` with `re, im ∈ ℤ`.
///
/// Exact: every ring operation is integer arithmetic. Arithmetic overflow
/// panics in debug builds and wraps in release, the same contract `i64` has
/// everywhere else in Rust.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, PartialOrd, Ord)]
pub struct GaussianInt {
    /// The real part.
    pub re: i64,
    /// The imaginary part.
    pub im: i64,
}

impl GaussianInt {
    /// `re + im·i`.
    #[inline]
    pub const fn new(re: i64, im: i64) -> Self {
        Self { re, im }
    }

    /// The rational integer `n + 0i`.
    #[inline]
    pub const fn from_int(n: i64) -> Self {
        Self { re: n, im: 0 }
    }

    /// The **algebraic** norm `N(a + bi) = a² + b²`, exact in `ℤ` and strictly
    /// multiplicative.
    ///
    /// [`Coefficient::magnitude`] returns its square root (the modulus) because
    /// the trait's `N(xy) = N(x)N(y)` law is stated for the modulus; the field
    /// norm here is the exact integer companion, useful when a test wants no
    /// floating point at all.
    #[inline]
    pub const fn norm_sq(self) -> i64 {
        self.re * self.re + self.im * self.im
    }
}

impl std::fmt::Display for GaussianInt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.im {
            0 => write!(f, "{}", self.re),
            im if im < 0 => write!(f, "{}-{}i", self.re, -im),
            im => write!(f, "{}+{}i", self.re, im),
        }
    }
}

impl std::ops::Add for GaussianInt {
    type Output = GaussianInt;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl std::ops::AddAssign for GaussianInt {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.re += rhs.re;
        self.im += rhs.im;
    }
}

impl std::ops::Sub for GaussianInt {
    type Output = GaussianInt;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl std::ops::Neg for GaussianInt {
    type Output = GaussianInt;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}

impl std::ops::Mul for GaussianInt {
    type Output = GaussianInt;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl std::ops::MulAssign for GaussianInt {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl std::iter::Sum for GaussianInt {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::default(), |a, b| a + b)
    }
}

impl num::Zero for GaussianInt {
    #[inline]
    fn zero() -> Self {
        Self::new(0, 0)
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.re == 0 && self.im == 0
    }
}

impl num::One for GaussianInt {
    #[inline]
    fn one() -> Self {
        Self::new(1, 0)
    }
}

impl Coefficient for GaussianInt {
    /// Exact integer sign flip — no float, no allocation.
    #[inline]
    fn mul_sign(&self, sign: i8) -> Self {
        if sign < 0 { -*self } else { *self }
    }

    /// The complex modulus `√(a² + b²)`.
    ///
    /// A genuine **absolute value** on `ℤ[i]`, so unlike [`crate::Term`] this
    /// impl satisfies every clause of the documented law: nonnegative, zero only
    /// at `0`, subadditive (triangle inequality) and multiplicative
    /// (`|zw| = |z||w|`). The `f64` here is a *readout* for a `Policy` threshold,
    /// not part of the value representation — the ring itself stays exact. The
    /// exact integer companion is [`GaussianInt::norm_sq`].
    #[inline]
    fn magnitude(&self) -> f64 {
        ((self.norm_sq()) as f64).sqrt()
    }
}

impl ImaginaryUnit for GaussianInt {
    #[inline]
    fn imaginary_unit() -> Self {
        Self::new(0, 1)
    }

    /// The exact component swap `(re, im) ↦ (−im, re)`.
    #[inline]
    fn mul_i(&self) -> Self {
        Self::new(-self.im, self.re)
    }
}

impl Conjugate for GaussianInt {
    #[inline]
    fn conj(&self) -> Self {
        Self::new(self.re, -self.im)
    }
}

impl From<i64> for GaussianInt {
    #[inline]
    fn from(value: i64) -> Self {
        Self::from_int(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::{One, Zero};

    #[test]
    fn imaginary_unit_law_holds_exactly() {
        // `i·i == −one()` — `lean/PPVM/Pauli/Matrix.lean` `iU_sq`. Exact ring, so
        // this is `assert_eq!` with zero tolerance.
        let i = GaussianInt::imaginary_unit();
        assert_eq!(i * i, -GaussianInt::one());
        assert_eq!(i * i * i * i, GaussianInt::one());
    }

    #[test]
    fn conj_of_i_is_minus_i_exactly() {
        // `lean/PPVM/Pauli/Matrix.lean` `star_iU`.
        let i = GaussianInt::imaginary_unit();
        assert_eq!(i.conj(), -i);
    }

    #[test]
    fn star_ring_laws() {
        let a = GaussianInt::new(3, -5);
        let b = GaussianInt::new(-2, 7);
        assert_eq!(a.conj().conj(), a);
        assert_eq!((a + b).conj(), a.conj() + b.conj());
        assert_eq!((a * b).conj(), a.conj() * b.conj());
    }

    #[test]
    fn magnitude_is_multiplicative_and_zero_only_at_zero() {
        let a = GaussianInt::new(3, 4);
        let b = GaussianInt::new(1, 2);
        assert_eq!(a.magnitude(), 5.0);
        assert!(((a * b).magnitude() - a.magnitude() * b.magnitude()).abs() < 1e-12);
        assert_eq!(GaussianInt::zero().magnitude(), 0.0);
        assert!(GaussianInt::new(0, 1).magnitude() > 0.0);
        // Subadditivity (the clause the ℓ¹ truncation bound needs).
        assert!((a + b).magnitude() <= a.magnitude() + b.magnitude() + 1e-12);
    }

    #[test]
    fn ring_is_associative_and_distributive_exactly() {
        let a = GaussianInt::new(2, -1);
        let b = GaussianInt::new(-3, 4);
        let c = GaussianInt::new(5, 6);
        assert_eq!((a * b) * c, a * (b * c));
        assert_eq!(a * (b + c), a * b + a * c);
    }

    #[test]
    fn mul_i_agrees_with_the_ring_product() {
        for re in -3..=3i64 {
            for im in -3..=3i64 {
                let z = GaussianInt::new(re, im);
                assert_eq!(z.mul_i(), z * GaussianInt::imaginary_unit());
            }
        }
    }

    #[test]
    fn mul_sign_is_exact() {
        let z = GaussianInt::new(7, -9);
        assert_eq!(z.mul_sign(1), z);
        assert_eq!(z.mul_sign(-1), -z);
    }

    #[test]
    fn display_round_trips_the_sign() {
        assert_eq!(GaussianInt::new(1, 0).to_string(), "1");
        assert_eq!(GaussianInt::new(1, 2).to_string(), "1+2i");
        assert_eq!(GaussianInt::new(1, -2).to_string(), "1-2i");
    }
}
