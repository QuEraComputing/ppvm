use std::hash::BuildHasher;

use crate::{pattern::PauliPattern, traits::*, word::PauliWord};
use dashmap::DashMap;
use rayon::prelude::*;

impl<'a, S, V, Hasher> ACMapBase for DashMap<PauliWord<S, Hasher>, V, Hasher>
where
    S: PauliStorage,
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

impl<'a, S, V, Hasher> ACMapAddAssign<S, V, Hasher> for DashMap<PauliWord<S, Hasher>, V, Hasher>
where
    S: PauliStorage + 'a,
    V: Coefficient + Sync + Send + 'a,
    Hasher: Default + Clone + BuildHasher + Sync + Send + 'a,
{
    fn add_assign(&mut self, key: PauliWord<S, Hasher>, value: V) {
        self.entry(key)
            .and_modify(|v| *v += value.clone())
            .or_insert(value);
    }

    fn map_add_assign<F>(&self, dest: &mut Self, f: F)
    where
        F: Fn(&PauliWord<S, Hasher>, &V) -> (PauliWord<S, Hasher>, V) + Sync + Send,
    {
        self.par_iter().for_each(|entry| {
            let (new_k, new_v) = f(entry.key(), entry.value());
            dest.entry(new_k)
                .and_modify(|v| *v += new_v.clone())
                .or_insert(new_v);
        });
    }
}

impl<'a, S, V, Hasher> ACMapMulAssign<V, Hasher> for DashMap<PauliWord<S, Hasher>, V, Hasher>
where
    S: PauliStorage + 'a,
    V: Coefficient + Send + Sync + 'a,
    Hasher: Default + Clone + BuildHasher + Sync + Send + 'a,
{
    fn mul_assign(&mut self, value: V) {
        self.par_iter_mut()
            .for_each(|mut v| *v.value_mut() *= value.clone());
    }
}

impl<'a, S, V, Hasher> ACMapIter<'a> for DashMap<PauliWord<S, Hasher>, V, Hasher>
where
    S: PauliStorage + 'a,
    V: Coefficient + 'a,
    Hasher: Default + Clone + BuildHasher + 'a,
{
    type Item = dashmap::mapref::multiple::RefMulti<'a, PauliWord<S, Hasher>, V>;
    type Iter = dashmap::iter::Iter<
        'a,
        PauliWord<S, Hasher>,
        V,
        Hasher,
        DashMap<PauliWord<S, Hasher>, V, Hasher>,
    >;

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

impl<'a, S, C, Hasher> Trace<'a, PauliWord<S, Hasher>, C>
    for DashMap<PauliWord<S, Hasher>, C, Hasher>
where
    S: PauliStorage + 'a,
    C: Coefficient + Send + Sync + 'a + From<bool>,
    Hasher: Clone + BuildHasher + Default + Send + Sync + 'a,
{
    fn trace(&'a self, value: &'a PauliWord<S, Hasher>) -> C {
        self.par_iter()
            .map(|entry| {
                let tr: C = value.trace(entry.key());
                tr * entry.value().clone()
            })
            .sum()
    }
}

impl<'a, S, C, State> Trace<'a, PauliPattern, C> for DashMap<PauliWord<S>, C, State>
where
    S: PauliStorage + 'a,
    C: Coefficient + Send + Sync + 'a + From<bool>,
    State: Clone + BuildHasher + 'a + Send + Sync,
{
    fn trace(&'a self, value: &'a PauliPattern) -> C {
        self.par_iter()
            .map(|entry| {
                let tr: C = value.trace(entry.key());
                tr * entry.value().clone()
            })
            .sum()
    }
}

impl<'a, S, C, H> ACMapConsume for DashMap<PauliWord<S, H>, C, H>
where
    S: PauliStorage + 'a,
    C: Coefficient + Send + Sync + 'a,
    H: Clone + BuildHasher + 'a + Send + Sync,
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

impl<'a, S, C, H> ACMapInsert<S, C, H> for DashMap<PauliWord<S, H>, C, H>
where
    S: PauliStorage + 'a,
    C: Coefficient + Send + Sync + 'a,
    H: Default + Clone + BuildHasher + 'a + Send + Sync,
{
    fn map_insert<F>(&mut self, dest: &mut Self, f: F)
    where
        F: Fn(&PauliWord<S, H>, &mut C) -> Option<(PauliWord<S, H>, C)> + Sync + Send,
    {
        self.par_iter_mut().for_each(|mut entry| {
            let (k, v) = entry.pair_mut();
            if let Some((new_k, new_v)) = f(k, v) {
                dest.insert(new_k, new_v);
            }
        })
    }
}

impl<S, C, H> ACMapContains<S, C, H> for DashMap<PauliWord<S, H>, C, H>
where
    S: PauliStorage,
    C: Coefficient + PartialEq,
    H: BuildHasher + Clone + Default,
{
    fn contains_with<F>(&self, key: &PauliWord<S, H>, f: F) -> bool
    where
        F: Fn(&C) -> bool,
    {
        self.get(key).map_or(false, |v| f(v.value()))
    }
}

impl<S, C, H> ACMapScale<S, C, H> for DashMap<PauliWord<S, H>, C, H>
where
    S: PauliStorage,
    C: Coefficient + Send + Sync,
    H: BuildHasher + Clone + Default + Sync + Send,
{
    fn scale<F>(&mut self, f: F)
    where
        F: Fn(&PauliWord<S, H>, &mut C) + Sync + Send,
    {
        self.par_iter_mut().for_each(|mut entry| {
            let (k, v) = entry.pair_mut();
            f(k, v);
        });
    }
}

impl<S, C, H> ACMapRetain<S, C, H> for DashMap<PauliWord<S, H>, C, H>
where
    S: PauliStorage,
    C: Coefficient,
    H: BuildHasher + Clone + Default + Sync + Send,
{
    fn retain<F>(&mut self, f: F)
    where
        F: Fn(&PauliWord<S, H>, &C) -> bool + Sync + Send,
    {
        Self::retain(self, |k, v| f(k, v));
    }
}
