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
    ///
    /// # Behaviour parity
    ///
    /// The `±i` arms go through [`ImaginaryUnit::mul_i`], **not** through a bare
    /// `c * imaginary_unit()`. On `Complex<f64>` the two are *not* the same
    /// function: the ring multiply computes `(re·0 − im·1, re·1 + im·0)`, and
    /// `inf·0`/`NaN·0` are `NaN`, so `(inf + 0i)·i` is `NaN + inf·i` while the
    /// old `ppvm_traits::ComplexCoefficient::mul_phase` — which swapped the
    /// components by hand — gave `−0 + inf·i`. `mul_i` restores the old
    /// component swap (and with it the sign of zero, visible through `Display`
    /// and serialization). See `phase_apply_matches_old_mul_phase_encoding`.
    ///
    /// The whole fold is delegated to [`ImaginaryUnit::mul_i_pow`], which is an
    /// **override point**: a ring whose values carry `iᵏ` symbolically (the
    /// symbolic `Term` of `ppvm-sym-2`) folds `iᵏ` into its own representation
    /// rather than through the `±1`/`mul_i` arms, which is what old's
    /// `ComplexCoefficient::mul_phase` did. The default body *is* those arms, so
    /// nothing changes for `f64`/`Complex<f64>`/`GaussianInt`.
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

/// A key whose set carries a product — the (projective) group structure that
/// lifts `C[K]` from a module to an algebra.
///
/// The keys form a group only *up to phase*: the product is **not closed on
/// keys**, it emits an `iᵏ`, which is why `key_mul` returns `(Self, Phase)` and
/// why `C[PauliWord]` is a **2-cocycle-twisted** group algebra.
///
/// # Laws
///
/// Write `key_mul(u, v) = (u · v, i^{β(u,v)})`. Every impl must satisfy:
///
/// * the key product is associative: `(u · v) · w == u · (v · w)`;
/// * the phase exponent is a **2-cocycle**:
///   `β(u,v) + β(u·v, w) == β(v,w) + β(u, v·w)` in `ℤ/4`.
///
/// Under exactly those two hypotheses the twisted product on `C × K` is
/// associative for any commutative coefficient ring `C` with `i⁴ = 1` — proved
/// key-agnostically in `lean/PPVM/Algebra/Twisted.lean` (`gtmul_assoc`, over an
/// abstract `kmul` and `IsCocycle`), so the obligation is stated once and every
/// key discharges it. `PauliWord` does so via `Bool.xor_assoc` and
/// `lean/PPVM/Pauli/Phase.lean` (`phaseExp_cocycle`), recovered as the instance
/// in `phaseExp_isCocycle` / `tmul_assoc_of_gtmul` (with `tmul_assoc` the
/// concrete Pauli statement). A future ordered fermionic-word key must discharge
/// the same two hypotheses; it does **not** inherit associativity from the Pauli
/// proof.
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

    /// Multiply by `i`. Semantically `self * imaginary_unit()`, which is the
    /// default body — but it is an **override point**, because on a
    /// floating-point ring the generic product is not extensionally equal to the
    /// rotation it denotes.
    ///
    /// On `Complex<f64>`, `c * i` expands to
    /// `(re·0 − im·1, re·1 + im·0)`; `inf·0` and `NaN·0` are `NaN`, so a
    /// non-finite component contaminates *both* output components, and `re·0`
    /// also loses the sign of zero. Multiplication by `i` is really the
    /// component swap `(re, im) ↦ (−im, re)`, which is total and exact — this is
    /// what the old `ppvm_traits::ComplexCoefficient::mul_phase` did, and the
    /// `Complex<f64>` impl below restores it verbatim. It is also cheaper (two
    /// negations instead of four multiplies and two adds).
    #[inline]
    fn mul_i(&self) -> Self {
        self.clone() * Self::imaginary_unit()
    }

    /// Multiply by `iᵏ` (`k` taken mod 4) — the fold [`Phase::apply`] delegates
    /// to, and the second **override point**.
    ///
    /// The default body is the four-arm `{clone, mul_i, neg, neg∘mul_i}` fold,
    /// which is the only sensible spelling on a ring whose values are numbers.
    /// It is overridable because a ring whose values carry the `iᵏ` *as data* —
    /// the symbolic `Term` of `ppvm-sym-2`, whose monomials hold a `ℤ/4` phase
    /// byte — must fold the phase into that representation instead: old's
    /// `ComplexCoefficient::mul_phase` promoted `Const(c)` to `One(i⁰, c)`
    /// **unconditionally**, including at `k = 0`, and `Term`'s `PartialEq` and
    /// `Display` are representational, so taking the `clone()` arm at `k = 0`
    /// would be a user-visible divergence from old.
    ///
    /// Impls must satisfy `x.mul_i_pow(k) == iᵏ · x` denotationally, and
    /// `mul_i_pow(1) == mul_i`.
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
