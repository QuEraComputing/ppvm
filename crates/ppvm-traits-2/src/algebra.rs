// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Algebra capabilities: key products, coefficient conjugation, and imaginary units.
//! [`Phase`] carries the residual fourth root of unity emitted by a key product.

use crate::arithmetic::Coefficient;

/// A fourth root of unity `iᵏ`, with `k` modulo four.
/// [`KeyProduct::key_mul`] returns this residual phase for the coefficient to absorb.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Phase {
    /// `i⁰ = +1`.
    Pos1,
    /// `i¹ = +i`.
    PosI,
    /// `i² = −1`.
    Neg1,
    /// `i³ = −i`.
    NegI,
}

impl Phase {
    /// The exponent `k ∈ {0, 1, 2, 3}` such that this phase equals `iᵏ`.
    #[inline]
    pub fn exponent(self) -> u8 {
        self as u8
    }

    /// The phase `iᵏ` for exponent `k` (taken mod 4).
    #[inline]
    pub fn from_exponent(k: u8) -> Self {
        match k & 3 {
            0 => Phase::Pos1,
            1 => Phase::PosI,
            2 => Phase::Neg1,
            _ => Phase::NegI,
        }
    }

    /// The identity phase: `p.compose(Self::one()) == p`.
    #[inline]
    pub fn one() -> Self {
        Phase::Pos1
    }

    /// Compose phases by adding exponents modulo four: `iᵃ · iᵇ = i^{a+b}`.
    /// Composition is commutative and can accumulate key-product phases before application.
    #[inline]
    pub fn compose(self, other: Self) -> Self {
        Phase::from_exponent(self.exponent() + other.exponent())
    }

    /// The group inverse `i^{-k} = i^{4-k}`, i.e. the phase `q` with
    /// `self.compose(q) == Phase::one()`.
    #[inline]
    pub fn inverse(self) -> Self {
        // 4 - k is exact for k ∈ {0,1,2,3}; the mod-4 reduction in
        // `from_exponent` sends k = 0 back to 0.
        Phase::from_exponent(4 - self.exponent())
    }

    /// Apply `iᵏ` to a coefficient through [`ImaginaryUnit::mul_i_pow`].
    /// Overrides can preserve symbolic representations; complex multiplication by `i`
    /// uses component swaps to preserve signed zeros and avoid non-finite contamination.
    #[inline]
    pub fn apply<C: ImaginaryUnit>(self, c: &C) -> C {
        c.mul_i_pow(self.exponent())
    }
}

/// `iᵃ · iᵇ = i^{a+b}` — [`compose`](Phase::compose) as the `*` operator, so a
/// residual-phase accumulator reads `acc *= phase` / `acc = a * b`.
impl core::ops::Mul for Phase {
    type Output = Phase;

    #[inline]
    fn mul(self, rhs: Phase) -> Phase {
        self.compose(rhs)
    }
}

impl core::ops::MulAssign for Phase {
    #[inline]
    fn mul_assign(&mut self, rhs: Phase) {
        *self = self.compose(rhs);
    }
}

/// A key product returning a key and residual phase `i^{β(u,v)}`.
/// Key multiplication must be associative, and phases must satisfy the cocycle law:
/// `β(u,v) + β(u·v,w) == β(v,w) + β(u,v·w)` modulo four.
pub trait KeyProduct: Eq + Clone {
    /// Product of two keys, with the phase it produces (folded onto the coeff).
    fn key_mul(&self, other: &Self) -> (Self, Phase);
}

/// A commutative coefficient ring with an imaginary unit.
/// Implementations must satisfy `Self::imaginary_unit() * Self::imaginary_unit()
/// == -Self::one()`, hence `i⁴ = 1`.
pub trait ImaginaryUnit: Coefficient + num::One {
    /// The imaginary unit `i`; impls must satisfy
    /// `Self::imaginary_unit() * Self::imaginary_unit() == -Self::one()`.
    fn imaginary_unit() -> Self;

    /// Multiply by `i`; the default uses ring multiplication.
    /// Complex values override this with `(re, im) ↦ (-im, re)` to preserve signed
    /// zeros and avoid contaminating both components when one is non-finite.
    #[inline]
    fn mul_i(&self) -> Self {
        self.clone() * Self::imaginary_unit()
    }

    /// Multiply by `iᵏ`, with `k` modulo four; used by [`Phase::apply`].
    /// Overrides may fold phases into symbolic data, including when `k == 0`.
    /// The result must denote `iᵏ · self`, with `mul_i_pow(1) == mul_i()`.
    #[inline]
    fn mul_i_pow(&self, k: u8) -> Self {
        match k & 3 {
            0 => self.clone(),
            1 => self.mul_i(),
            2 => -(self.clone()),
            _ => -(self.mul_i()),
        }
    }
}

/// A commutative ring involution used by [`crate::containers::Pair::hermitian_overlap`].
/// Requires `conj(conj(a)) == a`, preservation of addition and multiplication,
/// and `conj(i) == -i` when the ring also implements [`ImaginaryUnit`].
pub trait Conjugate: Coefficient {
    /// The ring involution applied to this value.
    fn conj(&self) -> Self;
}

impl ImaginaryUnit for num::Complex<f64> {
    #[inline]
    fn imaginary_unit() -> Self {
        num::Complex::new(0.0, 1.0)
    }

    /// The old `ComplexCoefficient::mul_phase(1)` component swap, verbatim
    /// (`crates/ppvm-traits/src/traits/coefficient.rs`): total on non-finite
    /// components and sign-of-zero exact, unlike the generic `self * i`.
    #[inline]
    fn mul_i(&self) -> Self {
        num::Complex::new(-self.im, self.re)
    }
}

impl Conjugate for num::Complex<f64> {
    #[inline]
    fn conj(&self) -> Self {
        num::Complex::conj(self)
    }
}

impl Conjugate for f64 {
    /// Conjugation is the identity on a real ring.
    #[inline]
    fn conj(&self) -> Self {
        *self
    }
}
