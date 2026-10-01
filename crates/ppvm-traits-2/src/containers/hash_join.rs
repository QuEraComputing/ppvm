// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! A `HashMap<K, C, IdentityBuildHasher>` backend for large supports.
//! [`Indexable`] keys supply structural digests consumed by the pass-through hasher.

use std::collections::HashMap;
use std::collections::hash_map::Entry;

use crate::algebra::{ImaginaryUnit, KeyProduct};
use crate::arithmetic::Coefficient;
use crate::containers::{
    Accumulate, IdentityBuildHasher, Indexable, Multiply, Pair, Retain, Scale, Support, TermBatch,
};

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
            if let Some(value) = self.get_mut(k) {
                *value += c;
            } else {
                self.insert(k.clone(), c.clone());
            }
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
            *v *= s;
        }
    }
}

/// Takes every [`Pair`] default: the shared scan already drives from the
/// smaller support and probes through [`Support::get`], which is a hash lookup
/// here.
impl<K, C> Pair for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable,
    C: Coefficient,
{
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
    /// Accumulate every key-pair product and its phase into a distinct `acc`.
    /// The outer product costs `O(|A|·|B|)`; reduction and truncation are explicit.
    /// Exact-zero entries remain until [`Accumulate::reduce`] runs.
    fn multiply_into(&self, other: &Self, acc: &mut Self) {
        for (p, a) in HashMap::iter(self) {
            for (q, b) in HashMap::iter(other) {
                let (k, phase) = p.key_mul(q);
                let mut product = a.clone();
                product *= b;
                let c = phase.apply(&product);
                match acc.entry(k) {
                    Entry::Occupied(mut slot) => *slot.get_mut() += c,
                    Entry::Vacant(slot) => {
                        slot.insert(c);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A minimal [`Indexable`] key, so the hash backend is testable here (the
    /// only real one, `PauliWord`, lives downstream). `Hash` is exactly
    /// `write_u64(key_hash())`, as the contract requires.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    struct Key(u64);

    impl std::hash::Hash for Key {
        fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
            state.write_u64(self.key_hash());
        }
    }

    impl Indexable for Key {
        fn key_hash(&self) -> u64 {
            self.0.wrapping_mul(0x9E37_79B9_7F4A_7C15)
        }
    }

    #[test]
    fn borrowed_arithmetic_matches_vector_backend() {
        use crate::TermSink;
        use num::Complex;

        let mut terms = TermBatch::new();
        terms.push(Key(1), Complex::new(1.0, 2.0));
        terms.push(Key(1), Complex::new(2.0, -1.0));
        terms.push(Key(2), Complex::new(-1.0, 3.0));
        let mut map = HashMap::<Key, Complex<f64>, IdentityBuildHasher>::default();
        let mut vector = Vec::<(Key, Complex<f64>)>::new();
        map.accumulate_batch(&terms);
        vector.accumulate_batch(&terms);
        let scalar = Complex::new(2.0, -1.0);
        map.scale(&scalar);
        vector.scale(&scalar);
        assert_eq!(Support::get(&map, &Key(1)), Some(Complex::new(7.0, -1.0)));

        let other_vector = vec![(Key(1), Complex::new(1.0, 1.0))];
        let other_map = other_vector.as_slice().iter().copied().collect();
        assert_eq!(map.overlap(&other_map), Complex::new(8.0, 6.0));
        assert_eq!(map.hermitian_overlap(&other_map), Complex::new(6.0, 8.0));
        assert_eq!(map.overlap(&other_map), vector.overlap(&other_vector));
        assert_eq!(other_map.overlap(&map), other_vector.overlap(&vector));
        assert_eq!(
            map.hermitian_overlap(&other_map),
            vector.hermitian_overlap(&other_vector)
        );
        assert_eq!(
            other_map.hermitian_overlap(&map),
            other_vector.hermitian_overlap(&vector)
        );
    }

    #[test]
    fn for_each_ref_agrees_with_iter() {
        let mut m: HashMap<Key, f64, IdentityBuildHasher> = HashMap::default();
        for (k, c) in [(1u64, 1.0), (2, 2.0), (3, -3.0)] {
            m.insert(Key(k), c);
        }
        let mut seen: Vec<(Key, f64)> = Vec::new();
        m.for_each_ref(|k, c| seen.push((*k, *c)));
        seen.sort_by_key(|(k, _)| k.0);
        let mut want: Vec<(Key, f64)> = Support::iter(&m).collect();
        want.sort_by_key(|(k, _)| k.0);
        assert_eq!(seen, want);
        assert_eq!(seen.len(), 3);
    }
}
