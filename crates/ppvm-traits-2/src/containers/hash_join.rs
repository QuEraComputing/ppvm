// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `HashMap<K, C, IdentityBuildHasher>` — the hash-join backend: `accumulate`
//! is a probe-and-merge against the bucket table, requiring `K: Indexable` (the
//! direct structural digest, consumed pass-through through
//! [`IdentityBuildHasher`]). Best for large support, e.g. `PauliSum`.
//!
//! See [`super`] for the shared design references and the orphan-rule note.

use std::collections::HashMap;

use crate::algebra::{Conjugate, ImaginaryUnit, KeyProduct};
use crate::batch::{KeyBatch, TermBatch};
use crate::coefficient::Coefficient;
use crate::graded::{Accumulate, Multiply, Pair, Retain, Scale, Support};
use crate::hash::{IdentityBuildHasher, Indexable};

impl<K, C> Support for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable,
    C: Coefficient,
{
    type Key = K;
    type Coeff = C;

    #[inline]
    fn len(&self) -> usize {
        HashMap::len(self)
    }

    #[inline]
    fn get(&self, key: &K) -> Option<C> {
        HashMap::get(self, key).cloned()
    }

    #[inline]
    fn iter(&self) -> impl Iterator<Item = (K, C)> {
        HashMap::iter(self).map(|(k, v)| (k.clone(), v.clone()))
    }

    /// The borrowing scan: hands out `(&K, &C)` straight from the buckets, so a
    /// filtering reader never clones a coefficient it is about to reject. Same
    /// order as [`Support::iter`] (the map's own bucket order).
    #[inline]
    fn for_each_ref(&self, mut f: impl FnMut(&K, &C)) {
        for (k, v) in HashMap::iter(self) {
            f(k, v);
        }
    }
}

impl<K, C> Accumulate for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable,
    C: Coefficient,
{
    /// Build side of the hash join: probe each produced term, accumulate its
    /// coefficient onto a matching key, insert on a miss.
    #[inline]
    fn accumulate_batch(&mut self, terms: &TermBatch<K, C>) {
        for (k, c) in terms.iter() {
            self.entry(k.clone())
                .and_modify(|e| *e += c.clone())
                .or_insert_with(|| c.clone());
        }
    }

    /// Drop every zero-coefficient key (`reduce_structural`).
    #[inline]
    fn reduce(&mut self) {
        self.retain(|_, v| !v.is_zero());
    }
}

impl<K, C> Scale for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn scale(&mut self, s: &C) {
        for v in self.values_mut() {
            *v *= s.clone();
        }
    }
}

impl<K, C> Pair for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn probe_batch(&self, keys: &KeyBatch<K>, out: &mut [Option<C>]) {
        debug_assert!(out.len() >= keys.keys().len());
        for (slot, k) in out.iter_mut().zip(keys.keys().iter()) {
            *slot = HashMap::get(self, k).cloned();
        }
    }

    /// `Σ_k self[k]·other[k]`, driven from the **smaller** support.
    ///
    /// The shared support is contained in both, so scanning either side and
    /// probing the other is `O(1)` per candidate and yields the same pair set;
    /// walking the smaller one makes the cost `O(min(|a|, |b|))` instead of
    /// `O(|self|)`. Against old's `data().trace(k)`-per-term
    /// (`ppvm-pauli-sum/src/sum/trace.rs`, a full linear scan of `self` per term
    /// of `other` — anti-feature 13) this is the intended asymptotic improvement;
    /// picking the smaller side is the second half of it, and it matters when the
    /// operands are deliberately unequal (`|a| = 10⁵`, `|b| = 10`).
    ///
    /// The *value* is unchanged — the pairing is symmetric and the left factor
    /// stays `self`'s coefficient either way — only the float **summation order**
    /// depends on the direction, which is why the differential bar on `overlap` is
    /// relative (`1e-12`) rather than bit-exact.
    #[inline]
    fn overlap(&self, other: &Self) -> C {
        if self.len() <= other.len() {
            HashMap::iter(self)
                .filter_map(|(k, a)| HashMap::get(other, k).map(|b| a.clone() * b.clone()))
                .sum()
        } else {
            HashMap::iter(other)
                .filter_map(|(k, b)| HashMap::get(self, k).map(|a| a.clone() * b.clone()))
                .sum()
        }
    }

    /// `Σ_k conj(self[k])·other[k]`, driven from the smaller support — see
    /// [`Pair::overlap`] for why the direction is free to differ.
    #[inline]
    fn hermitian_overlap(&self, other: &Self) -> C
    where
        C: Conjugate,
    {
        if self.len() <= other.len() {
            HashMap::iter(self)
                .filter_map(|(k, a)| HashMap::get(other, k).map(|b| a.conj() * b.clone()))
                .sum()
        } else {
            HashMap::iter(other)
                .filter_map(|(k, b)| HashMap::get(self, k).map(|a| a.conj() * b.clone()))
                .sum()
        }
    }
}

impl<K, C> Retain<K, C> for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn retain(&mut self, keep: impl Fn(&K, &C) -> bool) {
        // Inherent `HashMap::retain` shadows the trait method; no recursion.
        self.retain(|k, v| keep(k, v));
    }
}

impl<K, C> Multiply for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable + KeyProduct,
    C: ImaginaryUnit,
{
    /// The twisted convolution `(A·B)[k] = Σ_{p·q = k} A[p]·B[q]·i^{β(p,q)}`,
    /// accumulated into `acc` through the hash join — `twistedConv` of
    /// `lean/PPVM/Algebra/Twisted.lean` (monomial case `tmul`; associative by
    /// `tmul_assoc`/`gtmul_assoc`, and `(A·B)[I] = ⟨A, B⟩` by
    /// `twistedConv_apply_id`, tying L4 back to [`Pair::overlap`]).
    ///
    /// Every `(p, q)` pair contributes: the outer product is `O(|A|·|B|)` and is
    /// accumulated **into a distinct `acc`**, never folded back into an operand.
    /// (That is the bilinearity old's `MulAssign<PauliSum>` loses — see
    /// `crate::graded::Multiply` and `ppvm-pauli-sum-2::multiply`.)
    ///
    /// Neither `reduce` nor any truncation runs here: an exact-zero cancellation
    /// stays in `acc`'s support. Canonicalization is the caller's explicit
    /// [`Accumulate::reduce`], per §"`reduce()` is first-class, and runs only at
    /// finalize".
    fn multiply_into(&self, other: &Self, acc: &mut Self) {
        for (p, a) in HashMap::iter(self) {
            for (q, b) in HashMap::iter(other) {
                let (k, phase) = p.key_mul(q);
                let c = phase.apply(&(a.clone() * b.clone()));
                // Warm the fresh product key's structural digest *before* the
                // probe, for the reason `RotateInPlace` does: `key_mul` returns a
                // key with an empty hash cache, and letting the finalize fold fire
                // lazily inside `entry()` puts its mul-chain latency on the
                // bucket-index critical path with nothing to hide it. Semantic
                // no-op (identical digest).
                let _ = k.key_hash();
                acc.entry(k).and_modify(|e| *e += c.clone()).or_insert(c);
            }
        }
    }
}
