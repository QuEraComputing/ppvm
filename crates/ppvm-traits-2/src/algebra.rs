// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Algebra capabilities that lift `C[K]` from a module to a (twisted) algebra:
//! the key product [`KeyProduct`] and the two coefficient-ring capabilities it
//! and the sesquilinear pairing need, [`ImaginaryUnit`] and [`Conjugate`], plus
//! the [`Phase`] the Pauli key product emits.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"The map is a graded
//! algebra over `C[K]`" (L4 and its coefficient capabilities).

use crate::coefficient::Coefficient;

/// A fourth root of unity `iᵏ`, `k ∈ ℤ/4`, i.e. an element of `{1, i, −1, −i}`.
///
/// The Pauli key product is not closed on keys: `v·w = iᵏ (v⊕w)`, so
/// [`KeyProduct::key_mul`] returns the residual `Phase` for the coefficient to
/// absorb. `iᵏ` already spans `{1, i, −1, −i}`, so no separate `±` is carried.
///
/// Design: §"The map is a graded algebra over `C[K]`" (`KeyProduct`). The
/// exponent `k` is the packed `phaseExp` of `lean/PPVM/Pauli/Phase.lean`
/// (`phaseExp_eq_ref`).
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
        match self {
            Phase::Pos1 => 0,
            Phase::PosI => 1,
            Phase::Neg1 => 2,
            Phase::NegI => 3,
        }
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

    /// The group identity `i⁰ = +1` of the `ℤ/4` phase group.
    ///
    /// The neutral element for [`compose`](Self::compose): `p.compose(one()) ==
    /// p` for every phase `p`.
    #[inline]
    pub fn one() -> Self {
        Phase::Pos1
    }

    /// The `ℤ/4` group product `iᵃ · iᵇ = i^{a+b}` (exponents added mod 4).
    ///
    /// This is the first-class group operation on phases: a `KeyProduct` chain
    /// can accumulate the residual phases each [`KeyProduct::key_mul`] emits with
    /// `compose` (or the equivalent [`Mul`] impl) without a coefficient in hand,
    /// deferring the [`apply`](Self::apply) fold onto the coefficient to the end.
    /// The group is abelian, so `a.compose(b) == b.compose(a)`.
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

    /// Fold this phase onto a coefficient: return `iᵏ · c`.
    ///
    /// This is the `iPow` fold of `lean/PPVM/Algebra/Twisted.lean` — the phase a
    /// [`KeyProduct::key_mul`] emits is absorbed by the coefficient here. Needs a
    /// primitive fourth root of unity, hence the [`ImaginaryUnit`] bound.
    #[inline]
    pub fn apply<C: ImaginaryUnit>(self, c: &C) -> C {
        match self {
            Phase::Pos1 => c.clone(),
            Phase::PosI => c.clone() * C::imaginary_unit(),
            Phase::Neg1 => -(c.clone()),
            Phase::NegI => -(c.clone() * C::imaginary_unit()),
        }
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

/// A key whose set carries a product — the (projective) group structure that
/// lifts `C[K]` from a module to an algebra.
///
/// The keys form a group only *up to phase*: the product is **not closed on
/// keys**, it emits an `iᵏ`, which is why `key_mul` returns `(Self, Phase)` and
/// why `C[PauliWord]` is a **2-cocycle-twisted** group algebra. That the phase
/// exponent is a genuine 2-cocycle (hence the twisted product is associative) is
/// machine-checked in `lean/PPVM/Pauli/Phase.lean` (`phaseExp_cocycle`) and
/// `lean/PPVM/Algebra/Twisted.lean` (`tmul_assoc`).
///
/// Design: §"The map is a graded algebra over `C[K]`" (`KeyProduct`).
pub trait KeyProduct: Eq + Clone {
    /// Product of two keys, with the phase it produces (folded onto the coeff).
    fn key_mul(&self, other: &Self) -> (Self, Phase);
}

/// The phase capability L4 needs, over a **commutative** coefficient ring: a
/// distinguished primitive fourth root of unity `i`.
///
/// Impls must satisfy `Self::imaginary_unit() * Self::imaginary_unit() ==
/// -Self::one()` (hence `i⁴ = 1`). This is strictly weaker than requiring
/// `Complex<f64>`: `GaussianInt` (`ℤ[i]`), `Complex<Rational>`, and cyclotomic
/// integers all satisfy it, so L4 does not foreclose exact Pauli multiplication.
///
/// Design: §"The map is a graded algebra over `C[K]`" (`ImaginaryUnit`). Law
/// machine-checked in `lean/PPVM/Pauli/Matrix.lean` (`iU_sq`: `iU * iU = -1`),
/// and the twisted product is associative over any commutative ring with
/// `i⁴ = 1` in `lean/PPVM/Algebra/Twisted.lean` (`tmul_assoc`).
pub trait ImaginaryUnit: Coefficient + num::One {
    /// The imaginary unit `i`; impls must satisfy
    /// `Self::imaginary_unit() * Self::imaginary_unit() == -Self::one()`.
    fn imaginary_unit() -> Self;
}

/// A coefficient ring carrying a ring involution (a commutative `*`-ring):
/// complex conjugation on `Complex<f64>` / `GaussianInt` / cyclotomic integers,
/// and the identity on real rings.
///
/// Supplies exactly the conjugation the sesquilinear
/// [`crate::graded::Pair::hermitian_overlap`] needs; nothing in propagation
/// requires it, so — like [`ImaginaryUnit`] — it is a separate capability, not a
/// `Coefficient` bound.
///
/// Laws (commutative `*`-ring): `conj(conj(a)) == a`, `conj(a + b) ==
/// conj(a) + conj(b)`, `conj(a · b) == conj(a) · conj(b)`; and when the ring is
/// also [`ImaginaryUnit`], `conj(i) == −i`.
///
/// Design: §"The map is a graded algebra over `C[K]`" (`Conjugate`). The
/// `conj(i) == −i` law is machine-checked in `lean/PPVM/Pauli/Matrix.lean`
/// (`star_iU`: `star iU = -iU`).
pub trait Conjugate: Coefficient {
    /// The ring involution applied to this value.
    fn conj(&self) -> Self;
}

impl ImaginaryUnit for num::Complex<f64> {
    #[inline]
    fn imaginary_unit() -> Self {
        num::Complex::new(0.0, 1.0)
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
