// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `Vec<(K, C)>` — the coordinate-list backend: an unsorted association list
//! scanned linearly, requiring only `K: Eq + Clone` (it never hashes). Best for
//! small support, e.g. the `GeneralizedTableau` amplitude vector.
//!
//! See [`super`] for the shared design references and the orphan-rule note.

use crate::algebra::{Conjugate, ImaginaryUnit, KeyProduct};
use crate::batch::{KeyBatch, TermBatch};
use crate::coefficient::Coefficient;
use crate::graded::{Accumulate, Multiply, Pair, Retain, Scale, Support};

impl<K, C> Support for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    type Key = K;
    type Coeff = C;

    #[inline]
    fn len(&self) -> usize {
        self.as_slice().len()
    }

    #[inline]
    fn get(&self, key: &K) -> Option<C> {
        self.as_slice()
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, c)| c.clone())
    }

    #[inline]
    fn iter(&self) -> impl Iterator<Item = (K, C)> {
        self.as_slice().iter().map(|(k, c)| (k.clone(), c.clone()))
    }

    /// The borrowing scan: no clone at all, so a reader that rejects most terms
    /// pays nothing for the ones it rejects.
    #[inline]
    fn for_each_ref(&self, mut f: impl FnMut(&K, &C)) {
        for (k, c) in self.as_slice() {
            f(k, c);
        }
    }
}

impl<K, C> Accumulate for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    /// Linear-scan hash-join: for each produced term, find the matching key and
    /// add onto it, else append. `O(n·m)` in the support size, which is the
    /// right cost model for the small support this backend targets.
    #[inline]
    fn accumulate_batch(&mut self, terms: &TermBatch<K, C>) {
        for (k, c) in terms.iter() {
            if let Some(slot) = self.iter_mut().find(|(ek, _)| ek == k) {
                slot.1 += c.clone();
            } else {
                self.push((k.clone(), c.clone()));
            }
        }
    }

    /// Drop every zero-coefficient term (`reduce_structural`): canonicalize to
    /// reduced finite support. Runs only at finalize, never inline.
    #[inline]
    fn reduce(&mut self) {
        self.retain(|(_, v)| !v.is_zero());
    }
}

impl<K, C> Scale for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    #[inline]
    fn scale(&mut self, s: &C) {
        for (_, v) in self.iter_mut() {
            *v *= s.clone();
        }
    }
}

impl<K, C> Pair for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    #[inline]
    fn probe_batch(&self, keys: &KeyBatch<K>, out: &mut [Option<C>]) {
        debug_assert!(out.len() >= keys.keys().len());
        for (slot, k) in out.iter_mut().zip(keys.keys().iter()) {
            *slot = Support::get(self, k);
        }
    }

    #[inline]
    fn overlap(&self, other: &Self) -> C {
        self.as_slice()
            .iter()
            .filter_map(|(k, a)| Support::get(other, k).map(|b| a.clone() * b))
            .sum()
    }

    #[inline]
    fn hermitian_overlap(&self, other: &Self) -> C
    where
        C: Conjugate,
    {
        self.as_slice()
            .iter()
            .filter_map(|(k, a)| Support::get(other, k).map(|b| a.conj() * b))
            .sum()
    }
}

impl<K, C> Retain<K, C> for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    #[inline]
    fn retain(&mut self, keep: impl Fn(&K, &C) -> bool) {
        // Inherent `Vec::retain` shadows the trait method (inherent-first
        // resolution), so this does not recurse.
        self.retain(|(k, v)| keep(k, v));
    }
}

impl<K, C> Multiply for Vec<(K, C)>
where
    K: KeyProduct,
    C: ImaginaryUnit,
{
    /// The twisted convolution `(A·B)[k] = Σ_{p·q = k} A[p]·B[q]·i^{β(p,q)}`,
    /// accumulated into `acc` — the coordinate-list spelling of `twistedConv`
    /// (`lean/PPVM/Algebra/Twisted.lean`), whose monomial case is `tmul`.
    ///
    /// Neither `reduce` nor any truncation runs: `acc` keeps an exact-zero
    /// cancellation, exactly as `twistedConv` (a finitely-supported map is
    /// canonicalized only by an explicit [`Accumulate::reduce`]).
    fn multiply_into(&self, other: &Self, acc: &mut Self) {
        for (p, a) in self.as_slice() {
            for (q, b) in other.as_slice() {
                let (k, phase) = p.key_mul(q);
                let c = phase.apply(&(a.clone() * b.clone()));
                if let Some(slot) = acc.iter_mut().find(|(ek, _)| *ek == k) {
                    slot.1 += c;
                } else {
                    acc.push((k, c));
                }
            }
        }
    }
}
