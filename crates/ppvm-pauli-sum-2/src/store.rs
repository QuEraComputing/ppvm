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

/// Replace each key with a unique new key, optionally changing its coefficient.
pub trait RekeyBijective: Support {
    /// The callback may update the coefficient in place and returns the new key,
    /// or `None` if the key is unchanged. Distinct keys must map to distinct keys.
    /// Zero coefficients remain stored; callers must preserve word widths.
    fn rekey_bijective(
        &mut self,
        transform: impl FnMut(&Self::Key, &mut Self::Coeff) -> Option<Self::Key>,
    );
}

/// Change coefficients without copying keys or generating terms.
pub trait UpdateCoefficients: Support {
    /// Visit each coefficient in place, preserving keys and explicit zeros.
    fn update_coefficients(&mut self, update: impl FnMut(&Self::Key, &mut Self::Coeff));
}

/// Update original coefficients, then accumulate generated terms.
pub trait BranchInPlace: Support {
    /// Visit each original term once; merge branches only after all visits finish.
    /// Collisions add coefficients, exact zeros remain stored, and widths must match.
    fn branch_in_place<I: IntoIterator<Item = (Self::Key, Self::Coeff)>>(
        &mut self,
        transform: impl FnMut(&Self::Key, &mut Self::Coeff) -> I,
    );
}

/// Hash-map storage whose keys supply their own structural digests.
/// Exact-zero coefficients remain present until [`Accumulate::reduce`] is called.
#[derive(Debug)]
pub struct HashMapStore<K, C> {
    terms: HashMap<K, C, IdentityBuildHasher>,
    // Empty between operations; buffers rekeyed or generated terms.
    buffer: Vec<(K, C)>,
}

impl<K: Clone, C: Clone> Clone for HashMapStore<K, C> {
    fn clone(&self) -> Self {
        Self {
            terms: self.terms.clone(),
            buffer: Vec::new(),
        }
    }
}

impl<K, C> StoreAlloc for HashMapStore<K, C> {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            terms: HashMap::with_capacity_and_hasher(capacity, IdentityBuildHasher),
            buffer: Vec::new(),
        }
    }

    fn clear(&mut self) {
        self.terms.clear();
        self.buffer.clear();
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

impl<K: Indexable, C: Coefficient> RekeyBijective for HashMapStore<K, C> {
    fn rekey_bijective(&mut self, mut transform: impl FnMut(&K, &mut C) -> Option<K>) {
        let mut moved = std::mem::take(&mut self.buffer);
        // Unchanged keys stay in place; only changed keys move to new buckets.
        self.terms.retain(|key, coeff| match transform(key, coeff) {
            None => true,
            Some(new_key) => {
                moved.push((new_key, coeff.clone()));
                false
            }
        });
        // Bijectivity: a changed key cannot land on a key that stayed in place.
        for (key, coeff) in moved.drain(..) {
            let replaced = self.terms.insert(key, coeff);
            debug_assert!(replaced.is_none(), "rekey must preserve distinct keys");
        }
        self.buffer = moved;
    }
}

impl<K: Indexable, C: Coefficient> UpdateCoefficients for HashMapStore<K, C> {
    fn update_coefficients(&mut self, mut update: impl FnMut(&K, &mut C)) {
        for (key, coeff) in &mut self.terms {
            update(key, coeff);
        }
    }
}

impl<K: Indexable, C: Coefficient> BranchInPlace for HashMapStore<K, C> {
    fn branch_in_place<I: IntoIterator<Item = (K, C)>>(
        &mut self,
        mut transform: impl FnMut(&K, &mut C) -> I,
    ) {
        let mut branches = std::mem::take(&mut self.buffer);
        for (key, coeff) in &mut self.terms {
            branches.extend(transform(key, coeff));
        }
        for (key, coeff) in branches.drain(..) {
            self.terms.accumulate(key, coeff);
        }
        self.buffer = branches;
    }
}

#[cfg(test)]
mod tests;
