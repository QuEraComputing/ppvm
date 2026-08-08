// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::ColumnStore;
use super::*;

impl<K, C> Support for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    type Key = K;
    type Coeff = C;

    #[inline]
    fn len(&self) -> usize {
        self.primary.len()
    }

    #[inline]
    fn get(&self, key: &K) -> Option<C> {
        let hash = key.key_hash();
        self.primary
            .find(key, hash)
            .map(|slot| self.primary.coeffs[slot].clone())
    }

    #[inline]
    fn iter(&self) -> impl Iterator<Item = (K, C)> {
        (0..self.primary.rows())
            .filter(|&i| self.primary.is_live(i))
            .map(move |i| (self.primary.key(i), self.primary.coeffs[i].clone()))
    }

    #[inline]
    fn for_each_ref(&self, mut f: impl FnMut(&K, &C)) {
        for i in 0..self.primary.rows() {
            if !self.primary.is_live(i) {
                continue;
            }
            let key = self.primary.key(i);
            f(&key, &self.primary.coeffs[i]);
        }
    }
}

impl<K, C> Accumulate for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    #[inline]
    fn accumulate_batch(&mut self, terms: &TermBatch<K, C>) {
        for (k, c) in terms.iter() {
            self.primary.add(k.clone(), c.clone());
        }
    }

    #[inline]
    fn reduce(&mut self) {
        self.primary.retain_coeffs(|c| !c.is_zero());
    }
}

impl<K, C> Scale for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    #[inline]
    fn scale(&mut self, s: &C) {
        if self.primary.is_dense() {
            for c in &mut self.primary.coeffs {
                *c *= s.clone();
            }
        } else {
            for i in 0..self.primary.rows() {
                if self.primary.is_live(i) {
                    self.primary.coeffs[i] *= s.clone();
                }
            }
        }
    }
}

impl<K, C> Pair for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    #[inline]
    fn probe_batch(&self, keys: &KeyBatch<K>, out: &mut [Option<C>]) {
        debug_assert!(out.len() >= keys.keys().len());
        let hashes = keys.hashes();
        for (i, (slot, k)) in out.iter_mut().zip(keys.keys().iter()).enumerate() {
            let hash = if hashes.len() == keys.keys().len() {
                hashes[i]
            } else {
                k.key_hash()
            };
            *slot = self
                .primary
                .find(k, hash)
                .map(|s| self.primary.coeffs[s].clone());
        }
    }

    #[inline]
    fn overlap(&self, other: &Self) -> C {
        let mut acc = C::zero();
        for i in 0..self.primary.rows() {
            if !self.primary.is_live(i) {
                continue;
            }
            let key = self.primary.key(i);
            if let Some(slot) = other.primary.find(&key, self.primary.hashes[i]) {
                acc += self.primary.coeffs[i].clone() * other.primary.coeffs[slot].clone();
            }
        }
        acc
    }

    #[inline]
    fn hermitian_overlap(&self, other: &Self) -> C
    where
        C: Conjugate,
    {
        let mut acc = C::zero();
        for i in 0..self.primary.rows() {
            if !self.primary.is_live(i) {
                continue;
            }
            let key = self.primary.key(i);
            if let Some(slot) = other.primary.find(&key, self.primary.hashes[i]) {
                acc += self.primary.coeffs[i].conj() * other.primary.coeffs[slot].clone();
            }
        }
        acc
    }
}

impl<K, C> Retain<K, C> for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    #[inline]
    fn retain(&mut self, keep: impl Fn(&K, &C) -> bool) {
        self.primary.retain_terms(keep);
    }
}

impl<K, C> Multiply for ColumnStore<K, C>
where
    K: Columnar + KeyProduct,
    C: ImaginaryUnit,
{
    fn multiply_into(&self, other: &Self, acc: &mut Self) {
        let hint = product_capacity_hint(self.primary.len(), other.primary.len());
        acc.primary
            .reserve_for_live_len(acc.primary.len().saturating_add(hint));
        for i in 0..self.primary.rows() {
            if !self.primary.is_live(i) {
                continue;
            }
            let p = self.primary.key(i);
            for j in 0..other.primary.rows() {
                if !other.primary.is_live(j) {
                    continue;
                }
                let q = other.primary.key(j);
                let (k, phase) = p.key_mul(&q);
                let c = phase
                    .apply(&(self.primary.coeffs[i].clone() * other.primary.coeffs[j].clone()));
                acc.primary.add(k, c);
            }
        }
    }
}
