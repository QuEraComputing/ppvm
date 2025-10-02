use std::{collections::HashMap, hash::BuildHasher};

use crate::{
    traits::{Coefficient, PauliStorage},
    word::PauliWord,
};

macro_rules! impl_acmap {
    ($name:ident) => {
        impl<'a, S, V, State> crate::traits::ACMap<S, V> for $name<PauliWord<S>, V, State>
        where
            S: PauliStorage,
            V: Coefficient,
            State: Clone + BuildHasher + Default,
        {
            fn with_capacity(capacity: usize) -> Self {
                HashMap::with_capacity_and_hasher(capacity, State::default())
            }

            fn len(&self) -> usize {
                self.len()
            }

            fn clear(&mut self) {
                self.clear();
            }
        }
    };
}

macro_rules! impl_acmap_add_assign {
    ($name:ident) => {
        impl<'a, S, V, State> crate::traits::ACMapAddAssign<S, V> for $name<PauliWord<S>, V, State>
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

            fn map_add_assign<F>(&self, dest: &mut Self, f: F)
            where
                F: Fn(&PauliWord<S>, &V) -> (PauliWord<S>, V) + Sync + Send,
            {
                for (k, v) in self.iter() {
                    let (new_k, new_v) = f(k, v);
                    dest.add_assign(new_k, new_v);
                }
            }
        }
    };
}

macro_rules! impl_acmap_mul_assign {
    ($name:ident) => {
        impl<'a, S, V, State> crate::traits::ACMapMulAssign<V> for $name<PauliWord<S>, V, State>
        where
            S: PauliStorage + 'a,
            V: Coefficient + std::ops::MulAssign + 'a,
            State: Clone + BuildHasher + 'a,
        {
            fn mul_assign(&mut self, value: V) {
                for v in self.values_mut() {
                    *v *= value.clone();
                }
            }
        }
    };
}

macro_rules! impl_acmap_iter {
    ($name:ident) => {
        impl<'a, S, V, State> crate::traits::ACMapIter<'a, S, V> for $name<PauliWord<S>, V, State>
        where
            S: PauliStorage + 'a,
            V: Coefficient + 'a,
            State: Clone + BuildHasher + 'a,
        {
            type Item = (&'a PauliWord<S>, &'a V);
            type Iter = std::collections::hash_map::Iter<'a, PauliWord<S>, V>;

            fn iter(&'a self) -> Self::Iter {
                Self::iter(self)
            }
        }

        impl<'a, S, V, State> crate::traits::ACMapIterMut<'a, S, V>
            for $name<PauliWord<S>, V, State>
        where
            S: PauliStorage + 'a,
            V: Coefficient + 'a,
            State: Clone + BuildHasher + 'a,
        {
            type Item = (&'a PauliWord<S>, &'a mut V);
            type IterMut = std::collections::hash_map::IterMut<'a, PauliWord<S>, V>;
            fn iter_mut(&'a mut self) -> Self::IterMut {
                Self::iter_mut(self)
            }
        }
    };
}

macro_rules! impl_acmap_trace {
    ($name:ident) => {
        impl<'a, P, S, C, State> crate::traits::Trace<'a, P> for $name<PauliWord<S>, C, State>
        where
            P: crate::traits::Trace<'a, PauliWord<S>, Output = bool> + 'a,
            S: PauliStorage + 'a,
            C: Coefficient + std::ops::AddAssign + num::Zero + Clone + 'a,
            State: Clone + BuildHasher + 'a,
        {
            type Output = C;
            fn trace(&'a self, value: &'a P) -> Self::Output {
                let sum = C::zero();
                self.iter().fold(sum, |mut acc, (k, v)| {
                    value.trace(k).then(|| acc += v.clone());
                    acc
                })
            }
        }
    };
}

impl_acmap!(HashMap);
impl_acmap_add_assign!(HashMap);
impl_acmap_mul_assign!(HashMap);
impl_acmap_iter!(HashMap);
impl_acmap_trace!(HashMap);
