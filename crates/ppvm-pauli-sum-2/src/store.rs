// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use hashbrown::{HashMap, hash_map::Entry};
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

/// Replace each key with a unique new key, optionally changing its coefficient.
pub trait RekeyBijective: Support {
    /// The callback must map distinct input keys to distinct output keys.
    /// Zero coefficients remain stored; callers must preserve word widths.
    fn rekey_bijective(
        &mut self,
        transform: impl FnMut(Self::Key, Self::Coeff) -> (Self::Key, Self::Coeff),
    );
}

/// Negate selected coefficients without moving keys or rebuilding storage.
pub trait SignFlipByKey: Support {
    /// Negate a coefficient exactly when `negate` returns true for its key.
    fn flip_sign_by_key(&mut self, negate: impl Fn(&Self::Key) -> bool);
}

/// Hash-map storage whose keys supply their own structural digests.
/// Exact-zero coefficients remain present until [`Accumulate::reduce`] is called.
#[derive(Debug)]
pub struct HashMapStore<K, C> {
    terms: HashMap<K, C, IdentityBuildHasher>,
    // Empty between operations; reused when a gate changes keys.
    scratch: HashMap<K, C, IdentityBuildHasher>,
}

impl<K: Clone, C: Clone> Clone for HashMapStore<K, C> {
    fn clone(&self) -> Self {
        Self {
            terms: self.terms.clone(),
            scratch: HashMap::with_hasher(IdentityBuildHasher),
        }
    }
}

impl<K, C> StoreAlloc for HashMapStore<K, C> {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            terms: HashMap::with_capacity_and_hasher(capacity, IdentityBuildHasher),
            scratch: HashMap::with_hasher(IdentityBuildHasher),
        }
    }

    fn clear(&mut self) {
        self.terms.clear();
        self.scratch.clear();
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
        self.terms.get(key).cloned()
    }

    fn iter(&self) -> impl Iterator<Item = (K, C)> {
        self.terms
            .iter()
            .map(|(key, coeff)| (key.clone(), coeff.clone()))
    }

    fn for_each_ref(&self, mut f: impl FnMut(&K, &C)) {
        for (key, coeff) in &self.terms {
            f(key, coeff);
        }
    }
}

impl<K, C> Accumulate for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    fn accumulate_batch(&mut self, terms: &TermBatch<K, C>) {
        for (key, coeff) in terms.iter() {
            if let Some(value) = self.terms.get_mut(key) {
                *value += coeff;
            } else {
                self.terms.insert(key.clone(), coeff.clone());
            }
        }
    }

    #[inline]
    fn accumulate(&mut self, key: K, coeff: C) {
        match self.terms.entry(key) {
            Entry::Occupied(mut entry) => *entry.get_mut() += coeff,
            Entry::Vacant(entry) => {
                entry.insert(coeff);
            }
        }
    }

    fn reduce(&mut self) {
        self.terms.retain(|_, coeff| !coeff.is_zero());
    }
}

impl<K, C> Scale for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    fn scale(&mut self, scalar: &C) {
        for coeff in self.terms.values_mut() {
            *coeff *= scalar;
        }
    }
}

impl<K, C> Pair for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
}

impl<K: Indexable, C: Coefficient> RekeyBijective for HashMapStore<K, C> {
    fn rekey_bijective(&mut self, mut transform: impl FnMut(K, C) -> (K, C)) {
        self.scratch.clear();
        self.scratch.reserve(self.terms.len());
        // Changed keys must be reinserted into their new hash buckets.
        for (key, coeff) in self.terms.drain() {
            let (new_key, new_coeff) = transform(key, coeff);
            let replaced = self.scratch.insert(new_key, new_coeff);
            debug_assert!(replaced.is_none(), "rekey must preserve distinct keys");
        }
        // Keep the emptied source map for the next gate to reuse.
        std::mem::swap(&mut self.terms, &mut self.scratch);
    }
}

impl<K: Indexable, C: Coefficient> SignFlipByKey for HashMapStore<K, C> {
    fn flip_sign_by_key(&mut self, negate: impl Fn(&K) -> bool) {
        for (key, coeff) in &mut self.terms {
            if negate(key) {
                coeff.mul_sign_assign(-1);
            }
        }
    }
}

impl<K, C> Retain<K, C> for HashMapStore<K, C> {
    fn retain(&mut self, keep: impl Fn(&K, &C) -> bool) {
        self.terms.retain(|key, coeff| keep(key, coeff));
    }
}

#[cfg(test)]
mod tests;
