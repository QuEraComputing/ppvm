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
    /// Prefer moving coefficients through a bijective map re-key instead of
    /// cloning them from a borrowed bucket traversal.
    ///
    /// The two traversals are algebraically identical. The borrowed default is
    /// faster for small numeric coefficients because it walks a dense, stable
    /// table; heap-backed coefficients whose clone walks owned state may opt in
    /// to draining traversal to transfer that state without allocation.
    const PREFER_MOVED_REKEY: bool = false;

    /// Multiply by `sign ∈ {-1, +1}` (encoded as `i8`).
    fn mul_sign(&self, sign: i8) -> Self;

    /// Multiply this coefficient in place by `sign ∈ {-1, +1}`.
    ///
    /// The default preserves `*self = self.mul_sign(sign)`. Heap-backed rings
    /// can override it to mutate their allocation without cloning it first.
    #[inline]
    fn mul_sign_assign(&mut self, sign: i8) {
        *self = self.mul_sign(sign);
    }

    /// Accumulate a borrowed coefficient.
    ///
    /// The default preserves the value semantics of `*self += rhs.clone()`.
    /// Heap-backed rings may override it to clone only the entries they actually
    /// insert instead of first cloning an entire temporary coefficient.
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

    /// Transfer eigenvalues `(λ_X, λ_Z, λ_Y)` for Pauli probabilities
    /// `(p_X, p_Y, p_Z)`.
    #[inline(always)]
    fn pauli_error_factors(probabilities: [Self; 3]) -> [Self; 3]
    where
        Self: Sized + num::One,
    {
        let [px, py, pz] = probabilities;
        let one = Self::one();
        [
            one.clone() - py.doubled() - pz.doubled(),
            one.clone() - px.doubled() - py.doubled(),
            one - px.doubled() - pz.doubled(),
        ]
    }

    /// Nonnegative magnitude. Exposes a property of the value for a `Policy` to
    /// threshold; it does not itself decide any cutoff. Replaces the old
    /// `Coefficient::cutoff`.
    ///
    /// # Law
    ///
    /// `magnitude` must be an **absolute value** `N` on the coefficient ring:
    ///
    /// * `N(x) >= 0`,
    /// * `N(x) == 0` iff `x == 0`,
    /// * `N(x + y) <= N(x) + N(y)` (subadditive),
    /// * `N(x * y) == N(x) * N(y)` (multiplicative).
    ///
    /// This is the assumption the entire truncation guarantee rests on, so it is
    /// an impl obligation rather than a hint. Nonnegativity alone does **not**
    /// suffice: `N(x) = x²` is nonnegative, vanishes only at `0`, and is even
    /// multiplicative, yet the `ℓ¹` bound fails for it — machine-checked in
    /// `lean/PPVM/Algebra/Truncation.lean` (`l1_bound_needs_subadditive`).
    ///
    /// Design: §"Coefficient, angle, and truncation". The threshold comparison
    /// that consumes this lives in `Policy` (`ppvm-pauli-sum-2`), whose keep-rule
    /// boundary against the tableau path is machine-checked in
    /// `lean/PPVM/Algebra/Truncation.lean` (`cutoff_mismatch`); the `ℓ¹`
    /// truncation-error bound this law buys is `l1_bound_abv` there (any
    /// coefficient ring with such an `N`), specialized to the shipped
    /// `Complex<f64>`/`norm()` configuration by `l1_bound_norm` /
    /// `l1_bound_complex`.
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

/// The complex-coefficient angle domain, i.e. the defaulted `A = C` case of
/// [`crate::gates::RotationOne`] at `C = Complex<f64>`.
///
/// # Behaviour parity
///
/// The old `ppvm_traits::Coefficient` carried `sin_cos` as a *method on the
/// coefficient*, and it was implemented for `Complex<f64>` as
/// `(Complex::new(sin(re), 0), Complex::new(cos(re), 0))`
/// (`crates/ppvm-traits/src/traits/coefficient.rs`): the real part is taken as
/// the angle and the amplitudes come back purely real. Splitting `sin_cos` out
/// into [`Angle`] must not silently drop that instantiation — without this impl
/// a `Sum` over `Complex<f64>` coefficients could not be rotated at all
/// (`RotationOne<Complex<f64>>` would have no angle type), a capability the old
/// crate had. Ported verbatim, imaginary part of `theta` discarded exactly as
/// before.
impl Angle<num::Complex<f64>> for num::Complex<f64> {
    #[inline]
    fn sin_cos(&self) -> (num::Complex<f64>, num::Complex<f64>) {
        let (s, c) = num::traits::Float::sin_cos(self.re);
        (num::Complex::new(s, 0.0), num::Complex::new(c, 0.0))
    }
}

/// A **real** angle driving a complex-coefficient sum.
///
/// # Behaviour parity
///
/// The old rotation methods took `theta: impl Into<T::Coeff>` and
/// `Coefficient: From<f64>`, so `sum.rx(0, 0.1)` compiled on a `Complex<f64>`
/// sum: the `f64` was widened to `Complex::new(0.1, 0.0)` and then `sin_cos`
/// dropped the (zero) imaginary part again. The redesigned
/// [`crate::gates::RotationOne`] cannot take `impl Into<A>` — `A` is a free trait
/// parameter, so the conversion would leave it uninferable at the call site — so
/// that caller-visible spelling is preserved *here* instead, by making `f64`
/// itself an angle domain over complex coefficients. The amplitudes it returns
/// are exactly what the old widening produced.
impl Angle<num::Complex<f64>> for f64 {
    #[inline]
    fn sin_cos(&self) -> (num::Complex<f64>, num::Complex<f64>) {
        let (s, c) = num::traits::Float::sin_cos(*self);
        (num::Complex::new(s, 0.0), num::Complex::new(c, 0.0))
    }
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

    #[inline(always)]
    fn pauli_error_factors([px, py, pz]: [Self; 3]) -> [Self; 3] {
        [
            1.0 - py * 2.0 - pz * 2.0,
            1.0 - px * 2.0 - py * 2.0,
            1.0 - px * 2.0 - pz * 2.0,
        ]
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

    #[inline(always)]
    fn doubled(&self) -> Self {
        *self * 2.0
    }

    #[inline(always)]
    fn pauli_error_factors([px, py, pz]: [Self; 3]) -> [Self; 3] {
        let one = Self::new(1.0, 0.0);
        [
            one - py * 2.0 - pz * 2.0,
            one - px * 2.0 - py * 2.0,
            one - px * 2.0 - pz * 2.0,
        ]
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
