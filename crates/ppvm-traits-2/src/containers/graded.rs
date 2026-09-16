// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Container algebra over `C[K]`: support, accumulation, scaling, pairing, and products.
//! Keys require `Eq + Clone`; hashing and Pauli capabilities are backend-specific.
//! [`Retain`] is separate because truncation can break algebraic exactness.

use crate::algebra::{Conjugate, ImaginaryUnit, KeyProduct};
use crate::arithmetic::Coefficient;
use crate::containers::{KeyBatch, TermBatch, TermSink};

/// A finitely supported map from keys to coefficients.
/// Read-only export avoids requiring mutable pair slots in columnar backends.
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

    /// Visit each supported key and coefficient by reference, in backend order.
    /// The default materializes owned pairs through [`Self::iter`]; stored backends
    /// should override it to avoid cloning coefficients before filtering.
    fn for_each_ref(&self, mut f: impl FnMut(&Self::Key, &Self::Coeff)) {
        for (k, c) in self.iter() {
            f(&k, &c);
        }
    }
}

/// Accumulate linear combinations and explicitly reduce exact-zero coefficients.
pub trait Accumulate: Support {
    /// Merge a batch, adding to existing keys or inserting new ones.
    /// The algebraic multiset contract permits reordering and partitioning;
    /// floating-point accumulation can still differ through rounding.
    fn accumulate_batch(&mut self, terms: &TermBatch<Self::Key, Self::Coeff>);

    /// Canonicalize to reduced finite-support form: drop every key whose
    /// coefficient `is_zero()`. First-class and run **only** at finalize — never
    /// inline during accumulation.
    fn reduce(&mut self);

    /// Accumulate one term through a singleton batch.
    fn accumulate(&mut self, key: Self::Key, coeff: Self::Coeff) {
        let mut batch = TermBatch::with_capacity(1);
        batch.push(key, coeff);
        self.accumulate_batch(&batch);
    }
}

/// Scale coefficients by a scalar without changing their keys.
pub trait Scale: Support {
    /// Multiply every coefficient by `s`: `∀ k. c_k *= s`.
    fn scale(&mut self, s: &Self::Coeff);
}

/// Pair supports: `overlap` computes `Σ a_k b_k` without conjugation.
/// `hermitian_overlap` computes `Σ conj(a_k) b_k` and requires [`Conjugate`].
/// For Pauli coefficients, the bilinear pairing is the normalized trace `Tr(A B)/2ⁿ`.
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

/// Compute `tr(self · value)` against a potentially different right-hand type.
/// Unlike [`Pair::overlap`], both the right-hand type and output are unconstrained.
pub trait Trace<'a, RHS: 'a> {
    /// Numeric output of the trace.
    type Output;
    /// Compute `tr(self · value)`.
    fn trace(&'a self, value: &'a RHS) -> Self::Output;
}

/// Accumulate a ring product using [`KeyProduct`] and [`ImaginaryUnit`].
/// Key-product phases are absorbed into coefficients; types without a product
/// need not implement this layer.
pub trait Multiply: Accumulate
where
    Self::Key: KeyProduct,
    Self::Coeff: ImaginaryUnit,
{
    /// Accumulate the ring product `self · other` into `acc`.
    fn multiply_into(&self, other: &Self, acc: &mut Self);
}

/// Retain terms selected by a predicate, independently of the algebraic layers.
/// Used by truncation policies; error guarantees depend on the coefficient magnitude.
pub trait Retain<W, C> {
    /// Retain exactly the terms for which `keep(&word, &coeff)` is `true`.
    fn retain(&mut self, keep: impl Fn(&W, &C) -> bool);
}
