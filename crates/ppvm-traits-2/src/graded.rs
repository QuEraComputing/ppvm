// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The graded map algebra over `C[K]`: the free `C`-module on a key set,
//! layered by algebraic strength. Each layer is a distinct trait justified by a
//! distinct algebraic property *and* a distinct consumer.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"The map is a graded algebra
//! over `C[K]`". The layers are `Support` (L0, the container), `Accumulate` (L1,
//! the module core), `Scale` (L2, the `C`-module action), `Pair` (L3, the trace
//! pairings), and `Multiply` (L4, the ring product). [`Retain`] sits *outside*
//! the algebra: dropping supported terms breaks module exactness, so it is a
//! non-algebraic capability that `Policy` (in `ppvm-pauli-sum-2`) — not the
//! algebra — consumes.
//!
//! The key type is bounded on `Eq + Clone`, **not** [`crate::hash::Indexable`]:
//! `C[K]` is the free module over any index set, so the algebra needs only a
//! valid map key. Hash backends re-add `Indexable` on *their* impls; Pauli
//! propagation re-adds `Word`/`PauliBits` on *its* methods.

use crate::algebra::{Conjugate, ImaginaryUnit, KeyProduct};
use crate::batch::{KeyBatch, TermBatch, TermSink};
use crate::coefficient::Coefficient;

/// L0 — the container: a finitely-supported function `K ⇀ C`.
///
/// No `&mut (K, C)` and no `&mut [C]` slot access is exposed, so a columnar
/// (structure-of-arrays) backend is expressible; `iter` is read-only export.
///
/// Design: §"The map is a graded algebra over `C[K]`" (L0 `Support`).
pub trait Support {
    /// The key type — minimal `Eq + Clone`; hash backends add `Indexable`.
    type Key: Eq + Clone;
    /// The coefficient ring.
    type Coeff: Coefficient;

    /// Number of terms in the (reduced) support.
    fn len(&self) -> usize;

    /// Whether the support is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The coefficient at `key`, if present.
    fn get(&self, key: &Self::Key) -> Option<Self::Coeff>;

    /// Read-only export of the support as `(key, coeff)` pairs. A SoA backend
    /// synthesizes the pairs from its columns.
    fn iter(&self) -> impl Iterator<Item = (Self::Key, Self::Coeff)>;
}

/// L1 — the module core: form linear combinations, then canonicalize.
///
/// Design: §"The map is a graded algebra over `C[K]`" (L1 `Accumulate`). The
/// module laws are machine-checked in `lean/PPVM/Algebra/GradedMap.lean`
/// (`accumulate_comm`, `accumulate_assoc`); `reduce` drops exactly the zero
/// coefficients (`reduce_structural`).
pub trait Accumulate: Support {
    /// Build side of the hash join: merge a produced batch, accumulating onto an
    /// existing key or inserting a new one. Columnar in.
    ///
    /// The batch is algebraically a **multiset** of terms — folding it in is a
    /// homomorphism from the free commutative monoid — so an impl is free to
    /// reorder it and to split it across partitions/threads: order-invariance
    /// (`accumulateTerms_perm`) and partition-invariance
    /// (`accumulateTerms_add`, `accumulateTerms B m` for `B = B1 + B2`) are
    /// machine-checked in `lean/PPVM/Algebra/GradedMap.lean`.
    fn accumulate_batch(&mut self, terms: &TermBatch<Self::Key, Self::Coeff>);

    /// Canonicalize to reduced finite-support form: drop every key whose
    /// coefficient `is_zero()`. First-class and run **only** at finalize — never
    /// inline during accumulation.
    fn reduce(&mut self);

    /// Scalar sugar over a batch of one — accumulate a single `(key, coeff)`.
    ///
    /// Design: §"The map is a graded algebra over `C[K]`" ("the scalar
    /// `accumulate(k, c)` is provided sugar over a batch of one"); that the
    /// singleton batch is the same operation is `accumulateTerms_singleton` in
    /// `lean/PPVM/Algebra/GradedMap.lean`.
    fn accumulate(&mut self, key: Self::Key, coeff: Self::Coeff) {
        let mut batch = TermBatch::with_capacity(1);
        batch.push(key, coeff);
        self.accumulate_batch(&batch);
    }
}

/// L2 — the `C`-module action: a pure elementwise map over the coefficients.
///
/// Design: §"The map is a graded algebra over `C[K]`" (L2 `Scale`). Machine-
/// checked in `lean/PPVM/Algebra/GradedMap.lean` (`scale_scale`,
/// `scale_accumulate`).
pub trait Scale: Support {
    /// Multiply every coefficient by `s`: `∀ k. c_k *= s`.
    fn scale(&mut self, s: &Self::Coeff);
}

/// L3 — the read side of the hash join. Two pairings live here, differing only
/// in whether the first operand is conjugated.
///
/// `overlap` is the **symmetric bilinear** Hilbert–Schmidt trace pairing
/// `⟨A, B⟩ = ∑_k a_k b_k` — bilinear, *not* conjugated. Full `C`-bilinearity
/// (biadditivity, homogeneity in each slot, and symmetry over a commutative
/// ring) is machine-checked in `lean/PPVM/Algebra/GradedMap.lean`
/// (`overlap_add_left`/`overlap_add_right`, `overlap_smul_left`/`overlap_smul_right`,
/// `overlap_comm`); Pauli-basis orthonormality is `overlap_single_single` in
/// `lean/PPVM/Algebra/Noise.lean`. The `= Tr(A B)/2ⁿ` reading is proved, not
/// asserted: `overlap_eq_trace_div` in `lean/PPVM/Pauli/Matrix.lean` gives
/// `Tr(Â B̂) = 2ⁿ · ⟨A, B⟩` for the genuine `2ⁿ×2ⁿ` operators over `ℤ[i]`
/// (via `trace_tensorPauli_mul`), and `twistedConv_apply_id` in
/// `lean/PPVM/Algebra/Twisted.lean` states the same fact inside `C[K]`:
/// `⟨A, B⟩` is the identity coefficient of the L4 product.
///
/// `hermitian_overlap` is the **sesquilinear** inner product
/// `⟨φ | ψ⟩ = ∑_k conj(a_k)·b_k`, conjugate-linear in the first argument, so it
/// requires the coefficient ring to carry a [`Conjugate`]. Conjugate symmetry,
/// sesquilinearity, and `⟨f, f⟩ ≥ 0` are machine-checked in
/// `lean/PPVM/Algebra/GradedMap.lean` (`hermitianOverlap_conj_symm`,
/// `hermitianOverlap_smul_left`/`smul_right`, `hermitianOverlap_self_nonneg`).
///
/// Design: §"The map is a graded algebra over `C[K]`" (L3 `Pair`).
pub trait Pair: Support {
    /// Read-only probe of a key column: `out[i]` is the coefficient at
    /// `keys[i]`, or `None` on a miss.
    fn probe_batch(&self, keys: &KeyBatch<Self::Key>, out: &mut [Option<Self::Coeff>]);

    /// The symmetric bilinear trace pairing `∑_k a_k b_k`.
    fn overlap(&self, other: &Self) -> Self::Coeff;

    /// The sesquilinear inner product `∑_k conj(a_k)·b_k`.
    fn hermitian_overlap(&self, other: &Self) -> Self::Coeff
    where
        Self::Coeff: Conjugate;
}

/// `tr(self · value)` against a *different* right-hand type.
///
/// The old `ppvm_traits::Trace` (`ppvm-traits/src/traits/trace.rs`), ported
/// unchanged. It is **not** subsumed by [`Pair::overlap`]: `overlap` pairs a map
/// with another map of the *same* type, while `Trace` is the heterogeneous form
/// the old crate used to trace a word or a whole sum against a `PauliPattern`
/// (`ppvm-pauli-word/src/pattern/trace.rs`, `ppvm-pauli-sum/src/sum/trace.rs`),
/// which is why the right-hand type and the numeric `Output` are both free.
///
/// Design: §"Compatibility with current names" (`Trace`). Its `-2` implementers
/// land with the pattern port; the same-type Pauli-sum case is `Pair::overlap`,
/// whose bilinearity is machine-checked in `lean/PPVM/Algebra/GradedMap.lean`.
pub trait Trace<'a, RHS: 'a> {
    /// Numeric output of the trace.
    type Output;
    /// Compute `tr(self · value)`.
    fn trace(&'a self, value: &'a RHS) -> Self::Output;
}

/// L4 — the ring product. The only layer that needs the *key* to carry a
/// product; it stays optional and is not implemented for a key type that has
/// none. The Pauli product injects powers of `i`, so the coefficient must absorb
/// phase — bounded on [`ImaginaryUnit`], the minimal requirement.
///
/// Design: §"The map is a graded algebra over `C[K]`" (L4 `Multiply`).
/// Associativity of the twisted product holds over any commutative ring with
/// `i⁴ = 1`, machine-checked in `lean/PPVM/Algebra/Twisted.lean` (`tmul_assoc`);
/// the basis-monomial product is `multiply_single` in
/// `lean/PPVM/Algebra/GradedMap.lean`. The whole-map product is `twistedConv`
/// (`lean/PPVM/Algebra/Twisted.lean`), whose identity-key coefficient is exactly
/// [`Pair::overlap`] (`twistedConv_apply_id`) — the spec tying L4 to L3.
pub trait Multiply: Accumulate
where
    Self::Key: KeyProduct,
    Self::Coeff: ImaginaryUnit,
{
    /// Accumulate the ring product `self · other` into `acc`.
    fn multiply_into(&self, other: &Self, acc: &mut Self);
}

/// The one non-algebraic map operation: keep only the terms a predicate
/// selects. Dropping supported terms breaks module exactness, so this lives
/// outside the graded algebra and is consumed by `Policy`, not the algebra.
///
/// Design: §"The map is a graded algebra over `C[K]`" and §"Truncation"
/// (`Policy::truncate` bounds on `Retain`). The truncation error incurred is
/// bounded in `lean/PPVM/Algebra/Truncation.lean` (`l1_bound` over `ℝ`,
/// `l1_bound_abv` over any coefficient ring whose
/// [`Coefficient::magnitude`] is an
/// absolute value — the law that bound consumes).
pub trait Retain<W, C> {
    /// Retain exactly the terms for which `keep(&word, &coeff)` is `true`.
    fn retain(&mut self, keep: impl Fn(&W, &C) -> bool);
}
