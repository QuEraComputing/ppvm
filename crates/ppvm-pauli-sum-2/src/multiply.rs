// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! L4 — the **operator product** on a [`Sum`]: the twisted convolution
//!
//! ```text
//! (A · B)[k] = Σ_{p·q = k} A[p] · B[q] · i^{β(p,q)}
//! ```
//!
//! where `(p·q, i^{β(p,q)}) = p.key_mul(q)`
//! ([`KeyProduct`]). The Pauli key product is closed on
//! keys only **up to phase**, so `C[PauliWord]` is a 2-cocycle-twisted group
//! algebra and the coefficient ring must be able to hold the emitted `iᵏ` — hence
//! the [`ImaginaryUnit`] bound, which is what makes the product *unavailable* on a
//! real ring like `f64` (§"Compatibility with current names": old bounded the same
//! operator on `ComplexCoefficient`, which `f64` likewise did not implement).
//!
//! Design: `traits-2-configuration-and-hashing.md` §"The map is a graded algebra
//! over `C[K]`" (L4 `Multiply`) and §"Every gate is a producer feeding
//! `accumulate`" (the multiply row: "outer product over two operands' support,
//! one term per `(v, w)` pair"). Lean spec: `twistedConv` in
//! `lean/PPVM/Algebra/Twisted.lean`, with
//!
//! * associativity `tmul_assoc` / `gtmul_assoc` (over any commutative ring with
//!   `i⁴ = 1`, from key-product associativity plus the `phaseExp` 2-cocycle law
//!   `phaseExp_isCocycle`), lifted to whole maps as `twistedConv_assoc`,
//! * biadditivity `twistedConv_add_left` / `twistedConv_add_right` — the law old
//!   violates (below), and the step that lifts the monomial `tmul_assoc` to the
//!   whole-map `twistedConv_assoc`,
//! * the unit laws `one_tmul` / `tmul_one`, and
//! * `twistedConv_apply_id`: `(A · B)[I] = ⟨A, B⟩` — the L4↔L3 tie, so
//!   [`Sum::overlap`] is literally the identity coefficient of this product.
//!
//! # Behaviour: this is the ONE deliberate divergence from the old crate
//!
//! Old's `impl MulAssign<PauliSum<T>> for PauliSum<T>`
//! (`crates/ppvm-pauli-sum/src/sum/ops.rs:70`) loops over the rhs terms calling
//! `self.map_add(..)` once per term — but `map_add` **replaces** the support with
//! its image (it writes the re-key into a cleared aux map and swaps). So after the
//! first rhs term the support is `A·b₀P₀`, and the second iteration multiplies
//! *that* by `b₁P₁`: old computes the product **chain** `A·b₀P₀·b₁P₁` instead of
//! the bilinear sum `A·b₀P₀ + A·b₁P₁`. It is non-bilinear for any rhs with more
//! than one term, and completely untested in old (a single-term rhs, and
//! `MulAssign<PauliWord>`, are unaffected and correct).
//!
//! Bilinearity of `twistedConv` in each argument is exactly what old violates, so
//! the golden master here is the Lean semantics, not old's output: every monomial
//! product is accumulated into a **fresh accumulator**, never folded back into an
//! operand. `product_is_bilinear_not_a_chain` below pins the divergence.
//!
//! # Behaviour preserved from the old crate
//!
//! * **No truncation.** Like every other operation, the product never invokes the
//!   policy; `Sum::truncate` remains the single place `policy.truncate` runs.
//! * **No zero-dropping.** No `reduce` runs either, so an exact cancellation
//!   (`A·B` producing `0·I`) stays in the support — matching old, which has no
//!   `reduce` at all and whose exact-map equality depends on zero terms surviving.
//!   Canonicalize explicitly with [`Accumulate::reduce`] if you want it.
//! * **Right-multiplication by a single word is a bijective re-key.**
//!   [`Sum::mul_word_assign`] ports old's `MulAssign<PauliWord>` verbatim: `p ↦
//!   p·q` is injective (`p₁⊕q = p₂⊕q ⟹ p₁ = p₂`), so it takes the plain-`insert`
//!   `RekeyBijective` path — no accumulation probe, no new allocation, the phase
//!   folded onto the coefficient exactly as old's `coeff.mul_phase(phase)`.
//! * **The product exists only on a ring that can hold `i`.** Old bounded
//!   `Mul`/`MulAssign` on `ComplexCoefficient`, which `f64` does not implement, so
//!   a real-coefficient `PauliSum` had *no product method at all* — a
//!   compile-time restriction, not a runtime error. The [`ImaginaryUnit`] bound
//!   here reproduces it exactly. On `Complex<f64>` the product is available:
//!
//! ```
//! use num::Complex;
//! use ppvm_pauli_sum_2::{PauliSum, PauliWord};
//!
//! let a = PauliSum::<Complex<f64>>::from_terms(
//!     1,
//!     [(PauliWord::from("X"), Complex::new(1.0, 0.0))],
//! );
//! let b = PauliSum::<Complex<f64>>::from_terms(
//!     1,
//!     [(PauliWord::from("Y"), Complex::new(1.0, 0.0))],
//! );
//! // X·Y = +i Z.
//! assert_eq!(a.multiply(&b).get(&PauliWord::from("Z")), Some(Complex::new(0.0, 1.0)));
//! ```
//!
//! …and on the real ring `f64` it does not compile, matching old:
//!
//! ```compile_fail
//! use ppvm_pauli_sum_2::{PauliSum, PauliWord};
//!
//! // `f64: ImaginaryUnit` does not hold (as `f64: ComplexCoefficient` did not in
//! // the old crate), so `multiply` is not in scope on a real-coefficient sum.
//! let a = PauliSum::<f64>::from_terms(1, [(PauliWord::from("X"), 1.0)]);
//! let b = PauliSum::<f64>::from_terms(1, [(PauliWord::from("Y"), 1.0)]);
//! let _ = a.multiply(&b);
//! ```
//!
//! The same holds for [`Sum::mul_word_assign`] and the `*=` / `*` operators:
//!
//! ```compile_fail
//! use ppvm_pauli_sum_2::{PauliSum, PauliWord};
//!
//! let mut a = PauliSum::<f64>::from_terms(1, [(PauliWord::from("X"), 1.0)]);
//! a.mul_word_assign(&PauliWord::from("Y"));
//! ```
//!
//! ```compile_fail
//! use ppvm_pauli_sum_2::{PauliSum, PauliWord};
//!
//! let mut a = PauliSum::<f64>::from_terms(1, [(PauliWord::from("X"), 1.0)]);
//! let b = PauliSum::<f64>::from_terms(1, [(PauliWord::from("Y"), 1.0)]);
//! a *= &b;
//! ```
//!
//! ```compile_fail
//! use ppvm_pauli_sum_2::{PauliSum, PauliWord};
//!
//! let a = PauliSum::<f64>::from_terms(1, [(PauliWord::from("X"), 1.0)]);
//! let b = PauliSum::<f64>::from_terms(1, [(PauliWord::from("Y"), 1.0)]);
//! let _ = &a * &b;
//! ```

use std::hash::BuildHasher;

use ppvm_pauli_word_2::{HashFinalize, PauliStorage, PauliWord};
use ppvm_traits_2::{Accumulate, ImaginaryUnit, Indexable, KeyProduct, Multiply, Word};

use crate::policy::Policy;
use crate::store::{MultiplyInPlace, RekeyBijective, StoreAlloc};
use crate::sum::Sum;

// --- The ring product (L4): needs `Multiply` on the storage. -----------------
impl<S, P> Sum<S, P>
where
    S: Multiply,
    S::Key: Word + Indexable + KeyProduct,
    S::Coeff: ImaginaryUnit,
    P: Policy<S::Key, S::Coeff>,
{
    /// Accumulate the operator product `self · other` **into** `acc`:
    /// `acc[k] += Σ_{p·q = k} self[p]·other[q]·i^{β(p,q)}`.
    ///
    /// Three-address by design (`acc` is a third sum, never an operand): the outer
    /// product writes to key `p·q`, which may collide with a key still to be read
    /// from either operand, so the product **cannot** be folded in place. Old's
    /// `MulAssign<PauliSum>` did fold in place and lost bilinearity — see the
    /// module docs. Use [`multiply`](Self::multiply) for a fresh product or
    /// [`multiply_in_place`](Self::multiply_in_place) for a correct `A *= B` that
    /// still reuses the store's double-buffer.
    ///
    /// Because it accumulates rather than replaces, calling it twice sums the two
    /// products (`acc += A·B; acc += A·C` is `A·(B + C)`) — the biadditivity
    /// `twistedConv` has and old does not.
    ///
    /// Runs **no** `reduce` and **no** truncation: exact-zero cancellations stay
    /// in `acc`'s support (`reduce()` is first-class and runs only at finalize).
    ///
    /// Design: §"The map is a graded algebra over `C[K]`" (L4 `Multiply`). Lean:
    /// `twistedConv` (`lean/PPVM/Algebra/Twisted.lean`), biadditive by
    /// `twistedConv_add_left` / `twistedConv_add_right` and associative on whole
    /// maps by `twistedConv_assoc` (the monomial `tmul_assoc` alone does not
    /// give that — bilinearity is what lifts it), with `(A·B)[I] = ⟨A, B⟩` by
    /// `twistedConv_apply_id`.
    #[inline]
    pub fn multiply_into(&self, other: &Self, acc: &mut Self) {
        debug_assert_eq!(
            self.n_sites(),
            other.n_sites(),
            "operator product requires equal-width operands"
        );
        debug_assert_eq!(
            self.n_sites(),
            acc.n_sites(),
            "operator product accumulator must have the operands' width"
        );
        Multiply::multiply_into(self.storage(), other.storage(), acc.storage_mut());
    }
}

// --- The allocating product: additionally needs `StoreAlloc`. ----------------
impl<S, P> Sum<S, P>
where
    S: Multiply + StoreAlloc,
    S::Key: Word + Indexable + KeyProduct,
    S::Coeff: ImaginaryUnit,
    P: Policy<S::Key, S::Coeff>,
{
    /// The operator product `self · other` as a fresh sum, carrying `self`'s
    /// policy (and its capacity hint). The old crate's `Mul`, minus its
    /// non-bilinearity — see the module docs.
    ///
    /// The policy's hint sizes the fresh sum, but the *accumulator* is then
    /// grown to the convolution's own estimate by
    /// [`Multiply::multiply_into`](ppvm_traits_2::Multiply) on the storage
    /// (`store::product_capacity_hint`) — the policy hint is a function of
    /// `n_sites` alone (`n·10`, or `min(4ⁿ/2, 2¹⁴)` for `NoPolicy`) and knows
    /// nothing of `|A|·|B|`, so on its own it would leave a `10³ × 10³` product
    /// walking a doubling chain of full rehashes over a growing ~10⁶-entry map.
    ///
    /// Runs no `reduce` and no truncation.
    #[inline]
    pub fn multiply(&self, other: &Self) -> Self {
        let mut acc = Self::with_policy(self.n_sites(), self.policy().clone());
        self.multiply_into(other, &mut acc);
        acc
    }
}

// --- In-place `A ← A·B`: needs the store's aux-backed accumulator. ------------
impl<S, P> Sum<S, P>
where
    S: Accumulate + MultiplyInPlace<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
    P: Policy<S::Key, S::Coeff>,
{
    /// Replace this sum with the operator product `self · other`, reusing the
    /// store's persistent double-buffer as the accumulator (no fresh map).
    ///
    /// This is `A *= B` done *correctly*: the whole convolution is accumulated
    /// into the aux and only then swapped in, so `A` is read complete throughout
    /// and bilinearity holds. Old's `MulAssign` swapped after **each rhs term**,
    /// which is what turned its sum into a product chain
    /// (`ppvm-pauli-sum/src/sum/ops.rs:70`; module docs).
    ///
    /// Runs no `reduce` and no truncation.
    #[inline]
    pub fn multiply_in_place(&mut self, other: &Self) {
        debug_assert_eq!(
            self.n_sites(),
            other.n_sites(),
            "operator product requires equal-width operands"
        );
        self.storage_mut().multiply_in_place(other.storage());
    }
}

// --- Right-multiply by a single key: a bijective re-key. ---------------------
impl<S, P> Sum<S, P>
where
    S: Accumulate + RekeyBijective<S::Key, S::Coeff>,
    // `Sync` is what lets the re-key closure — which borrows `rhs` — satisfy the
    // `Send + Sync` bound `RekeyBijective` carries for a partitioned backend
    // (architecture feature 12). It is not a new restriction on real keys:
    // `PauliStorage` is `Send + Sync` in `-2` exactly as in old, so `PauliWord`
    // and `LossyPauliWord` are both.
    S::Key: Word + Indexable + KeyProduct + Sync,
    S::Coeff: ImaginaryUnit,
    P: Policy<S::Key, S::Coeff>,
{
    /// Right-multiply every term by the single key `rhs`: `A ← A·rhs`.
    ///
    /// Old's `impl MulAssign<PauliWord<..>> for PauliSum<T>`
    /// (`ppvm-pauli-sum/src/sum/ops.rs:95`), ported: that path is a *single*
    /// `map_add` over a bijection and is believed correct, so its behaviour is
    /// preserved exactly — the support is **replaced** by its image, the emitted
    /// `iᵏ` is folded onto the coefficient (old's `coeff.mul_phase(phase)`, here
    /// `Phase::apply`, which delegates to `ImaginaryUnit::mul_i_pow` — the
    /// coefficient's *own* `iᵏ` fold, so a ring that carries the phase as data
    /// lands in old's representation and not merely at old's value), and nothing
    /// is truncated or dropped.
    ///
    /// `p ↦ p·rhs` is injective (`p₁ ⊕ rhs = p₂ ⊕ rhs ⟹ p₁ = p₂`; the key product
    /// is a group operation up to phase), so this takes the plain-`insert`
    /// `RekeyBijective` fast path — the same provably-injective rewrite rule the
    /// Clifford re-key uses, with no accumulation probe and no reallocation.
    ///
    /// Lean: the monomial product `tmul` of `lean/PPVM/Algebra/Twisted.lean`;
    /// injectivity of `p ↦ p·rhs` — the precondition the `RekeyBijective` plain
    /// `insert` rests on, since a collision would DROP a term rather than sum it
    /// — is `mulWord_right_injective`, and invertibility of a Pauli word
    /// (`P·P = +I`) is `phaseExpN_self`, both in `lean/PPVM/Pauli/Word.lean`.
    /// That this re-key *is* the L4 product against a one-term map — i.e. that
    /// this path and `multiply_into` agree — is `twistedConv_single_right`, and
    /// that it needs no aggregation probe is `twistedConv_single_right_apply`,
    /// both in `lean/PPVM/Algebra/Twisted.lean`.
    #[inline]
    pub fn mul_word_assign(&mut self, rhs: &S::Key) {
        debug_assert_eq!(
            self.n_sites(),
            rhs.n_sites(),
            "operator product requires equal-width operands"
        );
        self.rekey_bijective(|k: S::Key, c: S::Coeff| {
            let (word, phase) = k.key_mul(rhs);
            (word, phase.apply(&c))
        });
    }
}

/// `&A * &B` — the fresh operator product, i.e. [`Sum::multiply`].
impl<S, P> std::ops::Mul<&Sum<S, P>> for &Sum<S, P>
where
    S: Multiply + StoreAlloc,
    S::Key: Word + Indexable + KeyProduct,
    S::Coeff: ImaginaryUnit,
    P: Policy<S::Key, S::Coeff>,
{
    type Output = Sum<S, P>;

    #[inline]
    fn mul(self, rhs: &Sum<S, P>) -> Sum<S, P> {
        self.multiply(rhs)
    }
}

/// `A * B` — old's by-value operator spelling, with the Lean-correct bilinear
/// product rather than old's product chain.
impl<S, P> std::ops::Mul<Sum<S, P>> for Sum<S, P>
where
    S: Multiply + StoreAlloc,
    S::Key: Word + Indexable + KeyProduct,
    S::Coeff: ImaginaryUnit,
    P: Policy<S::Key, S::Coeff>,
{
    type Output = Sum<S, P>;

    #[inline]
    fn mul(self, rhs: Sum<S, P>) -> Sum<S, P> {
        self.multiply(&rhs)
    }
}

/// `A *= &B` — the in-place operator product, i.e. [`Sum::multiply_in_place`].
/// Bilinear, unlike old's `MulAssign` (module docs).
impl<S, P> std::ops::MulAssign<&Sum<S, P>> for Sum<S, P>
where
    S: Accumulate + MultiplyInPlace<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
    P: Policy<S::Key, S::Coeff>,
{
    #[inline]
    fn mul_assign(&mut self, rhs: &Sum<S, P>) {
        self.multiply_in_place(rhs);
    }
}

/// `A *= B` — the by-value form old exposed.
impl<S, P> std::ops::MulAssign<Sum<S, P>> for Sum<S, P>
where
    S: Accumulate + MultiplyInPlace<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
    P: Policy<S::Key, S::Coeff>,
{
    #[inline]
    fn mul_assign(&mut self, rhs: Sum<S, P>) {
        self.multiply_in_place(&rhs);
    }
}

/// Right-multiply by one ordinary Pauli word, preserving old's operator surface.
impl<A, H, S, P> std::ops::MulAssign<PauliWord<A, H>> for Sum<S, P>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize + Send + Sync,
    PauliWord<A, H>: KeyProduct,
    S: Accumulate<Key = PauliWord<A, H>> + RekeyBijective<PauliWord<A, H>, S::Coeff>,
    S::Coeff: ImaginaryUnit,
    P: Policy<PauliWord<A, H>, S::Coeff>,
{
    #[inline]
    fn mul_assign(&mut self, rhs: PauliWord<A, H>) {
        self.mul_word_assign(&rhs);
    }
}

/// `A * word` — by-value single-word right multiplication.
impl<A, H, S, P> std::ops::Mul<PauliWord<A, H>> for Sum<S, P>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize + Send + Sync,
    PauliWord<A, H>: KeyProduct,
    S: Accumulate<Key = PauliWord<A, H>> + RekeyBijective<PauliWord<A, H>, S::Coeff>,
    S::Coeff: ImaginaryUnit,
    P: Policy<PauliWord<A, H>, S::Coeff>,
{
    type Output = Sum<S, P>;

    #[inline]
    fn mul(mut self, rhs: PauliWord<A, H>) -> Sum<S, P> {
        self *= rhs;
        self
    }
}

#[cfg(test)]
mod tests {
    use num::Complex;

    use crate::{PauliSum, PauliWord};

    type CSum = PauliSum<Complex<f64>>;

    fn c(re: f64) -> Complex<f64> {
        Complex::new(re, 0.0)
    }

    fn word(s: &str) -> PauliWord {
        s.into()
    }

    fn sum(n: usize, terms: &[(&str, Complex<f64>)]) -> CSum {
        CSum::from_terms(n, terms.iter().map(|(s, v)| (word(s), *v)))
    }

    fn close(a: Complex<f64>, b: Complex<f64>) -> bool {
        (a - b).norm() <= 1e-12
    }

    /// Compare two sums as key sets + per-key coefficients (never zipped
    /// iteration order — the hash backend is unordered).
    fn assert_same(a: &CSum, b: &CSum, what: &str) {
        assert_eq!(a.len(), b.len(), "{what}: support size");
        for (k, v) in a.iter() {
            let w = b.get(&k).unwrap_or_else(|| panic!("{what}: missing {k}"));
            assert!(close(v, w), "{what}: {k} -> {v} vs {w}");
        }
    }

    /// The monomial case: `tmul` — one term times one term is the key product
    /// with the emitted `iᵏ` folded onto the coefficient product.
    #[test]
    fn monomial_product_is_key_mul() {
        // X · Y = +i Z.
        let a = sum(1, &[("X", c(2.0))]);
        let b = sum(1, &[("Y", c(3.0))]);
        let p = a.multiply(&b);
        assert_eq!(p.len(), 1);
        assert!(close(p.get(&word("Z")).unwrap(), Complex::new(0.0, 6.0)));
    }

    /// `one_tmul` / `tmul_one`: the identity word is the unit of the product.
    #[test]
    fn identity_is_the_unit() {
        let a = sum(3, &[("XYZ", c(1.5)), ("ZZI", c(-0.25)), ("IIX", c(2.0))]);
        let id = sum(3, &[("III", c(1.0))]);
        assert_same(&a.multiply(&id), &a, "A·I");
        assert_same(&id.multiply(&a), &a, "I·A");
    }

    /// `phaseExpN_self`: every Pauli word squares to `+I` (no residual phase).
    #[test]
    fn word_squares_to_identity() {
        for s in ["XYZ", "ZZZ", "IXY"] {
            let p = sum(3, &[(s, c(1.0))]);
            let sq = p.multiply(&p);
            assert_eq!(sq.len(), 1, "{s}²");
            assert!(close(sq.get(&word("III")).unwrap(), c(1.0)), "{s}²");
        }
    }

    /// **Bilinearity** — the property old's `MulAssign<PauliSum>` violates.
    ///
    /// `A·(B + C) == A·B + A·C`. Old (`ppvm-pauli-sum/src/sum/ops.rs:70`) loops
    /// the rhs terms through `map_add`, which *replaces* the support each time, so
    /// it computes the chain `A·b₀P₀·b₁P₁`. The Lean oracle is `twistedConv`
    /// (`lean/PPVM/Algebra/Twisted.lean`), whose biadditivity is machine-checked
    /// as `twistedConv_add_left` / `twistedConv_add_right`; this asserts the
    /// CORRECT (Lean) value, not old's.
    #[test]
    fn product_is_bilinear_not_a_chain() {
        let a = sum(2, &[("XZ", c(1.0)), ("IY", c(-0.5))]);
        let b = sum(2, &[("ZI", c(2.0))]);
        let d = sum(2, &[("IX", c(3.0))]);
        let b_plus_d = sum(2, &[("ZI", c(2.0)), ("IX", c(3.0))]);

        // A·(B + D)
        let lhs = a.multiply(&b_plus_d);
        // A·B + A·D, accumulated into one accumulator.
        let mut rhs = CSum::new(2);
        a.multiply_into(&b, &mut rhs);
        a.multiply_into(&d, &mut rhs);
        assert_same(&lhs, &rhs, "A·(B+D) vs A·B + A·D");

        // And it is NOT old's chain `(A·B)·D` — the divergence is observable, not
        // a coincidentally-equal rewrite.
        let chain = a.multiply(&b).multiply(&d);
        let same_as_chain = lhs.len() == chain.len()
            && lhs
                .iter()
                .all(|(k, v)| chain.get(&k).is_some_and(|w| close(v, w)));
        assert!(
            !same_as_chain,
            "the bilinear product must differ from old's product chain",
        );
    }

    /// `tmul_assoc` / `gtmul_assoc`: `(A·B)·C == A·(B·C)`.
    #[test]
    fn product_is_associative() {
        let a = sum(3, &[("XYZ", c(1.0)), ("ZIX", c(-2.0))]);
        let b = sum(3, &[("YYI", c(0.5)), ("IZZ", c(1.5))]);
        let d = sum(3, &[("ZXY", c(-1.0)), ("XXX", c(0.25))]);
        let mut left = a.multiply(&b).multiply(&d);
        let mut right = a.multiply(&b.multiply(&d));
        // The product deliberately keeps exact-zero cancellations, and the two
        // association orders can zero *different* keys; canonicalize explicitly
        // (the caller-driven `reduce`) before comparing supports.
        left.reduce();
        right.reduce();
        assert_same(&left, &right, "(A·B)·C vs A·(B·C)");
    }

    /// `twistedConv_apply_id`: `(A·B)[I] == ⟨A, B⟩` — the L4↔L3 tie.
    #[test]
    fn identity_coefficient_is_the_overlap() {
        let a = sum(3, &[("XYZ", c(1.0)), ("ZIX", c(-2.0)), ("IZI", c(0.5))]);
        let b = sum(3, &[("XYZ", c(3.0)), ("IZI", c(-1.0)), ("YYY", c(2.0))]);
        let p = a.multiply(&b);
        let id = p.get(&word("III")).unwrap();
        assert!(close(id, a.overlap(&b)), "{id} vs {}", a.overlap(&b));
    }

    /// The product must not drop an exact-zero cancellation (old has no `reduce`
    /// and its exact-map equality depends on zero terms surviving).
    #[test]
    fn exact_cancellation_stays_in_the_support() {
        let a = sum(1, &[("X", c(1.0))]);
        let mut acc = CSum::new(1);
        a.multiply_into(&sum(1, &[("X", c(1.0))]), &mut acc);
        a.multiply_into(&sum(1, &[("X", c(-1.0))]), &mut acc);
        assert_eq!(acc.len(), 1, "the cancelled key must still be present");
        assert_eq!(acc.get(&word("I")), Some(c(0.0)));
    }

    /// `multiply_in_place` (the `*=` operator) agrees with the allocating form.
    #[test]
    fn in_place_matches_allocating() {
        let a = sum(3, &[("XYZ", c(1.0)), ("ZIX", c(-2.0))]);
        let b = sum(3, &[("YYI", c(0.5)), ("IZZ", c(1.5))]);
        let expected = a.multiply(&b);
        let mut got = a.clone();
        got *= &b;
        assert_same(&got, &expected, "A *= B");
    }

    /// `mul_word_assign` is right-multiplication by a single word, and agrees with
    /// the general product against that word's one-term sum.
    #[test]
    fn mul_word_assign_matches_the_general_product() {
        let a = sum(3, &[("XYZ", c(1.0)), ("ZIX", c(-2.0)), ("III", c(0.5))]);
        let q = word("YZX");
        let expected = a.multiply(&sum(3, &[("YZX", c(1.0))]));
        let mut got = a.clone();
        got.mul_word_assign(&q);
        assert_same(&got, &expected, "A·q");
        // Bijective: the support size is preserved exactly.
        assert_eq!(got.len(), a.len());
    }

    #[test]
    fn by_value_operator_spellings_match_the_methods() {
        let a = sum(2, &[("XZ", c(1.0)), ("IY", c(-0.5))]);
        let b = sum(2, &[("ZI", c(2.0)), ("IX", c(3.0))]);
        let expected = a.multiply(&b);

        assert_same(&(a.clone() * b.clone()), &expected, "A * B");
        let mut assigned = a.clone();
        assigned *= b;
        assert_same(&assigned, &expected, "A *= B");

        let q = word("ZY");
        let expected_word = a.multiply(&sum(2, &[("ZY", c(1.0))]));
        assert_same(&(a.clone() * q.clone()), &expected_word, "A * q");
        let mut assigned_word = a;
        assigned_word *= q;
        assert_same(&assigned_word, &expected_word, "A *= q");
    }
}
