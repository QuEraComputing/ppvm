use std::hash::BuildHasher;

use crate::{pattern::PauliPattern, traits::*};
use dashmap::DashMap;
use rayon::prelude::*;

impl<'a, V, Hasher, T> ACMapBase for DashMap<T, V, Hasher>
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

impl<'a, S, V, Hasher, W> ACMapAddAssign<S, V, Hasher, W> for DashMap<W, V, Hasher>
where
    S: PauliStorage + 'a,
    V: Coefficient + Sync + Send + 'a,
    Hasher: Default + Clone + BuildHasher + Sync + Send + 'a,
    W: PauliWordTrait<S, Hasher> + Sync + Send + 'a,
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

// impl<'a, S, V, State> ACMapIterMut<'a, S, V> for DashMap<PauliWord<S>, V, State>
// where
//     S: PauliStorage + 'a,
//     V: Coefficient + 'a,
//     State: Clone + BuildHasher + 'a,
// {
//     type Item = dashmap::mapref::multiple::RefMutMulti<'a, PauliWord<S>, V>;
//     type IterMut =
//         dashmap::iter::IterMut<'a, PauliWord<S>, V, State, DashMap<PauliWord<S>, V, State>>;

//     fn iter_mut(&'a mut self) -> Self::IterMut {
//         DashMap::iter_mut(self)
//     }
// }

impl<'a, C, Hasher, W> Trace<'a, W> for DashMap<W, C, Hasher>
where
    C: Coefficient + Send + Sync + 'a,
    Hasher: Clone + BuildHasher + Default + Send + Sync + 'a,
    W: std::hash::Hash + std::cmp::Eq + Sync + Send + 'a,
    for<'b> W: Trace<'b, W, Output = bool>,
{
    type Output = C;
    fn trace(&'a self, value: &'a W) -> Self::Output {
        self.par_iter()
            .filter(|entry| value.trace(entry.key()))
            .map(|entry| entry.value().clone())
            .sum()
    }
}

impl<'a, C, State, W> Trace<'a, PauliPattern> for DashMap<W, C, State>
where
    C: Coefficient + Send + Sync + 'a,
    State: Clone + BuildHasher + 'a + Send + Sync,
    W: std::hash::Hash + std::cmp::Eq + Sync + Send + 'a,
    for<'b> W: Trace<'b, PauliPattern, Output = bool>,
{
    type Output = C;
    fn trace(&'a self, value: &'a PauliPattern) -> Self::Output {
        self.par_iter()
            .filter(|entry| entry.key().trace(value))
            .map(|entry| entry.value().clone())
            .sum()
    }
}

impl<'a, C, H, W> ACMapConsume for DashMap<W, C, H>
where
    C: Coefficient + Send + Sync + 'a,
    H: Clone + BuildHasher + 'a + Send + Sync,
    W: std::hash::Hash + std::cmp::Eq + Clone + Sync + Send + 'a,
{
    fn consume(&mut self, dest: &mut Self) {
        dest.par_iter().for_each(|entry| {
            self.entry(entry.key().clone())
                .and_modify(|v| *v += entry.value().clone())
                .or_insert_with(|| entry.value().clone());
        });

        // // FIXME: clone is not very efficient when T::Coeff is an expression
        // self.par_extend(
        //     dest.par_iter()
        //         .map(|entry| (entry.key().clone(), entry.value().clone())),
        // );
        dest.clear();
    }
}

impl<'a, S, C, H, W> ACMapInsert<S, C, H, W> for DashMap<W, C, H>
where
    S: PauliStorage + 'a,
    C: Coefficient + Send + Sync + 'a,
    H: Default + Clone + BuildHasher + 'a + Send + Sync,
    W: PauliWordTrait<S, H> + Send + Sync + 'a,
{
    fn map_insert<F>(&mut self, dest: &mut Self, f: F)
    where
        F: Fn(&W, &mut C) -> Option<(W, C)> + Sync + Send,
    {
        self.par_iter_mut().for_each(|mut entry| {
            let (k, v) = entry.pair_mut();
            if let Some((new_k, new_v)) = f(k, v) {
                dest.insert(new_k, new_v);
            }
        })
    }
}

impl<S, C, H, W> ACMapContains<S, C, H, W> for DashMap<W, C, H>
where
    S: PauliStorage,
    C: Coefficient + PartialEq,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait<S, H>,
{
    fn contains_with<F>(&self, key: &W, f: F) -> bool
    where
        F: Fn(&C) -> bool,
    {
        self.get(key).map_or(false, |v| f(v.value()))
    }
}

impl<S, C, H, W> ACMapScale<S, C, H, W> for DashMap<W, C, H>
where
    S: PauliStorage,
    C: Coefficient + Send + Sync,
    H: BuildHasher + Clone + Default + Sync + Send,
    W: PauliWordTrait<S, H> + Send + Sync,
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
    W: PauliWordTrait<S, H>,
{
    fn retain<F>(&mut self, f: F)
    where
        F: Fn(&W, &C) -> bool + Sync + Send,
    {
        Self::retain(self, |k, v| f(k, v));
    }
}
