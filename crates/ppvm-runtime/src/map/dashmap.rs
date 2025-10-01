use std::hash::BuildHasher;

use crate::{
    traits::{ACMap, ACMapAddAssign, ACMapIter, ACMapIterMut, ACMapMulAssign, Coefficient, PauliStorage},
    word::PauliWord,
};
use dashmap::DashMap;
use rayon::prelude::*;

impl<'a, S, V, State> ACMap<S, V> for DashMap<PauliWord<S>, V, State>
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
}

impl<'a, S, V, State> ACMapAddAssign<S, V> for DashMap<PauliWord<S>, V, State>
where
    S: PauliStorage + 'a,
    V: Coefficient + std::ops::AddAssign + 'a,
    State: Clone + BuildHasher + 'a,
{
    fn add_assign(&mut self, key: PauliWord<S>, value: V) {
        self.entry(key)
            .and_modify(|v| *v += value.clone())
            .or_insert(value);
    }
}

impl<'a, S, V, State> ACMapMulAssign<V> for DashMap<PauliWord<S>, V, State>
where
    S: PauliStorage + 'a,
    V: Coefficient + std::ops::MulAssign + Send + Sync + 'a,
    State: Clone + BuildHasher + 'a,
{
    fn mul_assign(&mut self, value: V) {
        self.par_iter_mut()
            .for_each(|mut v| *v.value_mut() *= value.clone());
    }
}

impl<'a, S, V, State> ACMapIter<'a, S, V> for DashMap<PauliWord<S>, V, State>
where
    S: PauliStorage + 'a,
    V: Coefficient + 'a,
    State: Clone + BuildHasher + 'a,
{
    type Item = dashmap::mapref::multiple::RefMulti<'a, PauliWord<S>, V>;
    type Iter =
        dashmap::iter::Iter<'a, PauliWord<S>, V, State, DashMap<PauliWord<S>, V, State>>;

    fn iter(&'a self) -> Self::Iter {
        DashMap::iter(self)
    }
}

impl<'a, S, V, State> ACMapIterMut<'a, S, V> for DashMap<PauliWord<S>, V, State>
where
    S: PauliStorage + 'a,
    V: Coefficient + 'a,
    State: Clone + BuildHasher + 'a,
{
    type Item = dashmap::mapref::multiple::RefMutMulti<'a, PauliWord<S>, V>;
    type IterMut =
        dashmap::iter::IterMut<'a, PauliWord<S>, V, State, DashMap<PauliWord<S>, V, State>>;

    fn iter_mut(&'a mut self) -> Self::IterMut {
        DashMap::iter_mut(self)
    }
}
