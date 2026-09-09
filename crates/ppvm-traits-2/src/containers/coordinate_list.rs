// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `Vec<(K, C)>` — the coordinate-list backend: an unsorted association list
//! scanned linearly, requiring only `K: Eq + Clone` (it never hashes). Best for
//! small support, e.g. the `GeneralizedTableau` amplitude vector.
//!
//! See [`super`] for the shared design references and the orphan-rule note.

use crate::algebra::{Conjugate, ImaginaryUnit, KeyProduct};
use crate::arithmetic::Coefficient;
use crate::containers::{Accumulate, Multiply, Pair, Retain, Scale, Support};
use crate::containers::{KeyBatch, TermBatch};

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::containers::TermSink;

    fn batch(terms: &[(&str, f64)]) -> TermBatch<String, f64> {
        let mut b = TermBatch::with_capacity(terms.len());
        for (k, c) in terms {
            b.push((*k).to_string(), *c);
        }
        b
    }

    #[test]
    fn vec_accumulate_combines_keys() {
        let mut v: Vec<(String, f64)> = Vec::new();
        v.accumulate_batch(&batch(&[("a", 1.0), ("b", 2.0), ("a", 3.0)]));
        assert_eq!(Support::get(&v, &"a".to_string()), Some(4.0));
        assert_eq!(Support::get(&v, &"b".to_string()), Some(2.0));
        assert_eq!(Support::len(&v), 2);
    }

    #[test]
    fn vec_reduce_drops_zero() {
        let mut v: Vec<(String, f64)> = Vec::new();
        v.accumulate_batch(&batch(&[("a", 1.0), ("a", -1.0), ("b", 2.0)]));
        v.reduce();
        assert_eq!(Support::len(&v), 1);
        assert_eq!(Support::get(&v, &"b".to_string()), Some(2.0));
    }

    #[test]
    fn vec_scale_and_overlap() {
        let mut a: Vec<(String, f64)> = Vec::new();
        a.accumulate_batch(&batch(&[("x", 2.0), ("y", 3.0)]));
        let mut b: Vec<(String, f64)> = Vec::new();
        b.accumulate_batch(&batch(&[("x", 5.0), ("z", 7.0)]));
        a.scale(&2.0);
        // overlap = (2*2)*5 = 20; y and z do not match.
        assert_eq!(Pair::overlap(&a, &b), 20.0);
    }

    #[test]
    fn for_each_ref_agrees_with_iter() {
        let terms = batch(&[("a", 1.0), ("b", 2.0), ("c", -3.0), ("a", 0.5)]);

        let mut v: Vec<(String, f64)> = Vec::new();
        v.accumulate_batch(&terms);
        let mut seen: Vec<(String, f64)> = Vec::new();
        v.for_each_ref(|k, c| seen.push((k.clone(), *c)));
        seen.sort_by(|a, b| a.0.cmp(&b.0));
        let mut want: Vec<(String, f64)> = Support::iter(&v).collect();
        want.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(seen, want);
        assert_eq!(seen.len(), 3);
    }

    #[test]
    fn vec_retain_filters() {
        let mut v: Vec<(String, f64)> = Vec::new();
        v.accumulate_batch(&batch(&[("keep", 2.0), ("drop", 0.5)]));
        Retain::retain(&mut v, |_, c| *c >= 1.0);
        assert_eq!(Support::len(&v), 1);
        assert_eq!(Support::get(&v, &"keep".to_string()), Some(2.0));
    }
}
