use std::hash::BuildHasher;

use crate::{pattern::PauliPattern, traits::*, word::PauliWord};
use dashmap::DashMap;
use rayon::prelude::*;

impl<'a, S, V, State> ACMapBase for DashMap<PauliWord<S>, V, State>
where
    S: PauliStorage,
    V: Coefficient,
    State: Clone + BuildHasher + Default,
{
    fn with_capacity(capacity: usize) -> Self {
        DashMap::with_capacity_and_hasher(capacity, State::default())
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn clear(&mut self) {
        Self::clear(self);
    }
}

impl<'a, S, V, State> ACMapAddAssign<S, V> for DashMap<PauliWord<S>, V, State>
where
    S: PauliStorage + 'a,
    V: Coefficient + Sync + Send + 'a,
    State: Clone + BuildHasher + Sync + Send + 'a,
{
    fn add_assign(&mut self, key: PauliWord<S>, value: V) {
        self.entry(key)
            .and_modify(|v| *v += value.clone())
            .or_insert(value);
    }

    fn map_add_assign<F>(&self, dest: &mut Self, f: F)
    where
        F: Fn(&PauliWord<S>, &V) -> (PauliWord<S>, V) + Sync + Send,
    {
        self.par_iter().for_each(|entry| {
            let (new_k, new_v) = f(entry.key(), entry.value());
            dest.entry(new_k)
                .and_modify(|v| *v += new_v.clone())
                .or_insert(new_v);
        });
    }
}

impl<'a, S, V, State> ACMapMulAssign<V> for DashMap<PauliWord<S>, V, State>
where
    S: PauliStorage + 'a,
    V: Coefficient + Send + Sync + 'a,
    State: Clone + BuildHasher + 'a,
{
    fn mul_assign(&mut self, value: V) {
        self.par_iter_mut()
            .for_each(|mut v| *v.value_mut() *= value.clone());
    }
}

impl<'a, S, V, State> ACMapIter<'a> for DashMap<PauliWord<S>, V, State>
where
    S: PauliStorage + 'a,
    V: Coefficient + 'a,
    State: Clone + BuildHasher + 'a,
{
    type Item = dashmap::mapref::multiple::RefMulti<'a, PauliWord<S>, V>;
    type Iter = dashmap::iter::Iter<'a, PauliWord<S>, V, State, DashMap<PauliWord<S>, V, State>>;

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

impl<'a, S, C, State> Trace<'a, PauliWord<S>> for DashMap<PauliWord<S>, C, State>
where
    S: PauliStorage + 'a,
    C: Coefficient + Send + Sync + 'a,
    State: Clone + BuildHasher + 'a + Send + Sync,
{
    type Output = C;
    fn trace(&'a self, value: &'a PauliWord<S>) -> Self::Output {
        self.par_iter()
            .filter(|entry| value.trace(entry.key()))
            .map(|entry| entry.value().clone())
            .sum()
    }
}

impl<'a, S, C, State> Trace<'a, PauliPattern> for DashMap<PauliWord<S>, C, State>
where
    S: PauliStorage + 'a,
    C: Coefficient + Send + Sync + 'a,
    State: Clone + BuildHasher + 'a + Send + Sync,
{
    type Output = C;
    fn trace(&'a self, value: &'a PauliPattern) -> Self::Output {
        self.par_iter()
            .filter(|entry| entry.key().trace(value))
            .map(|entry| entry.value().clone())
            .sum()
    }
}

impl<'a, S, C, State> ACMapCombineUnique for DashMap<PauliWord<S>, C, State>
where
    S: PauliStorage + 'a,
    C: Coefficient + Send + Sync + 'a,
    State: Clone + BuildHasher + 'a + Send + Sync,
{
    fn combine_unique(&mut self, dest: &mut Self) {
        // FIXME: clone is not very efficient when T::Coeff is an expression
        self.par_extend(
            dest.par_iter()
                .map(|entry| (entry.key().clone(), entry.value().clone())),
        );
        dest.clear();
    }
}

impl<'a, S, C, State> ACMapInsert<S, C> for DashMap<PauliWord<S>, C, State>
where
    S: PauliStorage + 'a,
    C: Coefficient + Send + Sync + 'a,
    State: Clone + BuildHasher + 'a + Send + Sync,
{
    fn map_insert<F>(&mut self, dest: &mut Self, f: F)
    where
        F: Fn(&PauliWord<S>, &mut C) -> Option<(PauliWord<S>, C)> + Sync + Send,
    {
        self.par_iter_mut().for_each(|mut entry| {
            let (k, v) = entry.pair_mut();
            if let Some((new_k, new_v)) = f(k, v) {
                dest.insert(new_k, new_v);
            }
        })
    }
}
