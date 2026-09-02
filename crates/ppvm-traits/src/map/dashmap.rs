// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::BuildHasher;

use crate::traits::*;
use dashmap::DashMap;
use rayon::prelude::*;

impl<V, Hasher, T> ACMapBase for DashMap<T, V, Hasher>
where
    T: std::hash::Hash + Eq,
    V: Coefficient,
    Hasher: Clone + BuildHasher + Default,
{
    fn with_capacity(capacity: usize) -> Self {
        DashMap::with_capacity_and_hasher(capacity, Hasher::default())
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn clear(&mut self) {
        Self::clear(self);
    }
}

impl<S, V, Hasher, W> ACMapAddAssign<S, V, Hasher, W> for DashMap<W, V, Hasher>
where
    S: PauliStorage,
    V: Coefficient + Sync + Send,
    Hasher: Default + Clone + BuildHasher + Sync + Send,
    W: PauliWordTrait + Sync + Send,
{
    fn add_assign(&mut self, key: W, value: V) {
        self.entry(key)
            .and_modify(|v| *v += value.clone())
            .or_insert(value);
    }

    fn map_add_assign<F>(&self, dest: &mut Self, f: F)
    where
        F: Fn(&W, &V) -> (W, V) + Sync + Send,
    {
        self.par_iter().for_each(|entry| {
            let (new_k, new_v) = f(entry.key(), entry.value());
            dest.entry(new_k)
                .and_modify(|v| *v += new_v.clone())
                .or_insert(new_v);
        });
    }
}

impl<V, Hasher, W> ACMapMulAssign<V, Hasher> for DashMap<W, V, Hasher>
where
    V: Coefficient + Send + Sync,
    Hasher: Default + Clone + BuildHasher + Sync + Send,
    W: std::hash::Hash + std::cmp::Eq + Send + Sync,
{
    fn mul_assign(&mut self, value: V) {
        self.par_iter_mut()
            .for_each(|mut v| *v.value_mut() *= value.clone());
    }
}

impl<'a, V, Hasher, W> ACMapIter<'a> for DashMap<W, V, Hasher>
where
    V: Coefficient + 'a,
    Hasher: Default + Clone + BuildHasher + 'a,
    W: std::hash::Hash + std::cmp::Eq + 'a,
{
    type Item = dashmap::mapref::multiple::RefMulti<'a, W, V>;
    type Iter = dashmap::iter::Iter<'a, W, V, Hasher, DashMap<W, V, Hasher>>;

    fn iter(&'a self) -> Self::Iter {
        DashMap::iter(self)
    }
}

impl<'a, P, C, Hasher, W> Trace<'a, P> for DashMap<W, C, Hasher>
where
    P: Sync + 'a,
    for<'b> P: Trace<'b, W, Output = bool>,
    C: Coefficient + Send + Sync + 'a,
    Hasher: Clone + std::hash::BuildHasher + Default + Send + Sync + 'a,
    W: std::hash::Hash + std::cmp::Eq + Sync + Send + 'a,
{
    type Output = C;
    fn trace(&'a self, value: &'a P) -> Self::Output {
        self.par_iter()
            .filter(|entry| value.trace(entry.key()))
            .map(|entry| entry.value().clone())
            .sum()
    }
}

impl<C, H, W> ACMapConsume for DashMap<W, C, H>
where
    C: Coefficient + Send + Sync,
    H: Clone + BuildHasher + Send + Sync,
    W: std::hash::Hash + std::cmp::Eq + Clone + Sync + Send,
{
    fn consume(&mut self, dest: &mut Self) {
        dest.par_iter().for_each(|entry| {
            self.entry(entry.key().clone())
                .and_modify(|v| *v += entry.value().clone())
                .or_insert_with(|| entry.value().clone());
        });

        dest.clear();
    }
}

impl<S, C, H, W> ACMapInsert<S, C, H, W> for DashMap<W, C, H>
where
    S: PauliStorage,
    C: Coefficient + Send + Sync,
    H: Default + Clone + BuildHasher + Send + Sync,
    W: PauliWordTrait + Send + Sync,
{
    fn map_insert_vec<F>(&mut self, dest: &mut Vec<(W, C)>, f: F)
    where
        F: Fn(&W, &mut C) -> Option<(W, C)> + Sync + Send,
    {
        use rayon::iter::ParallelExtend;
        dest.par_extend(self.par_iter_mut().filter_map(|mut entry| {
            let (k, v) = entry.pair_mut();
            f(k, v)
        }));
    }

    fn map_insert_multiple<F>(&mut self, dest: &mut Self, f: F)
    where
        F: Fn(&W, &mut C) -> Option<Vec<(W, C)>> + Sync + Send,
    {
        self.par_iter_mut().for_each(|mut entry| {
            let (k, v) = entry.pair_mut();
            if let Some(new_entries) = f(k, v) {
                for (new_k, new_v) in new_entries {
                    dest.insert(new_k, new_v);
                }
            }
        })
    }
}

impl<S, C, H, W> ACMapContains<S, C, H, W> for DashMap<W, C, H>
where
    S: PauliStorage,
    C: Coefficient + PartialEq,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait,
{
    fn contains_with<F>(&self, key: &W, f: F) -> bool
    where
        F: Fn(&C) -> bool,
    {
        self.get(key).is_some_and(|v| f(v.value()))
    }
}

impl<S, C, H, W> ACMapScale<S, C, H, W> for DashMap<W, C, H>
where
    S: PauliStorage,
    C: Coefficient + Send + Sync,
    H: BuildHasher + Clone + Default + Sync + Send,
    W: PauliWordTrait + Send + Sync,
{
    fn scale<F>(&mut self, f: F)
    where
        F: Fn(&W, &mut C) + Sync + Send,
    {
        self.par_iter_mut().for_each(|mut entry| {
            let (k, v) = entry.pair_mut();
            f(k, v);
        });
    }
}

impl<S, C, H, W> ACMapRetain<S, C, H, W> for DashMap<W, C, H>
where
    S: PauliStorage,
    C: Coefficient,
    H: BuildHasher + Clone + Default + Sync + Send,
    W: PauliWordTrait,
{
    fn retain<F>(&mut self, mut f: F)
    where
        F: FnMut(&W, &C) -> bool,
    {
        Self::retain(self, |k, v| f(k, v));
    }
}
