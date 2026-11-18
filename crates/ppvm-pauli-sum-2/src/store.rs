// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::collections::HashMap;

use ppvm_traits_2::{
    Accumulate, Coefficient, IdentityBuildHasher, Indexable, Pair, Retain, Scale, Support,
    TermBatch,
};

/// Allocation and reuse of sparse-sum storage.
pub trait StoreAlloc {
    /// Allocate storage for at least `capacity` terms.
    fn with_capacity(capacity: usize) -> Self;

    /// Remove all terms while retaining allocated capacity.
    fn clear(&mut self);
}

/// Hash-map storage whose keys supply their own structural digests.
/// Exact-zero coefficients remain present until [`Accumulate::reduce`] is called.
#[derive(Debug, Clone)]
pub struct HashMapStore<K, C> {
    terms: HashMap<K, C, IdentityBuildHasher>,
}

impl<K, C> StoreAlloc for HashMapStore<K, C> {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            terms: HashMap::with_capacity_and_hasher(capacity, IdentityBuildHasher),
        }
    }

    fn clear(&mut self) {
        self.terms.clear();
    }
}

impl<K, C> Support for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    type Key = K;
    type Coeff = C;

    #[inline]
    fn len(&self) -> usize {
        self.terms.len()
    }

    #[inline]
    fn get(&self, key: &K) -> Option<C> {
        Support::get(&self.terms, key)
    }

    fn iter(&self) -> impl Iterator<Item = (K, C)> {
        Support::iter(&self.terms)
    }

    fn for_each_ref(&self, f: impl FnMut(&K, &C)) {
        self.terms.for_each_ref(f);
    }
}

impl<K, C> Accumulate for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    fn accumulate_batch(&mut self, terms: &TermBatch<K, C>) {
        self.terms.accumulate_batch(terms);
    }

    #[inline]
    fn accumulate(&mut self, key: K, coeff: C) {
        self.terms.accumulate(key, coeff);
    }

    fn reduce(&mut self) {
        self.terms.reduce();
    }
}

impl<K, C> Scale for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    fn scale(&mut self, scalar: &C) {
        self.terms.scale(scalar);
    }
}

impl<K, C> Pair for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
}

impl<K: Indexable, C: Coefficient> Retain<K, C> for HashMapStore<K, C> {
    fn retain(&mut self, keep: impl Fn(&K, &C) -> bool) {
        Retain::retain(&mut self.terms, keep);
    }
}

#[cfg(test)]
mod tests;
