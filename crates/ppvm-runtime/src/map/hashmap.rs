use std::{collections::HashMap, hash::BuildHasher};

use crate::{
    traits::{ACMapCombineUnique, ACMapInsert, Coefficient, PauliStorage},
    word::PauliWord,
};

macro_rules! impl_acmap_base {
    ($($seg:ident)::+) => {
        impl<'a, S, V, State> crate::traits::ACMapBase for $($seg)::+<PauliWord<S>, V, State>
        where
            S: PauliStorage,
            V: Coefficient,
            State: Clone + BuildHasher + Default,
        {
            fn with_capacity(capacity: usize) -> Self {
                Self::with_capacity_and_hasher(capacity, State::default())
            }

            fn len(&self) -> usize {
                $($seg)::+::len(self)
            }

            fn clear(&mut self) {
                $($seg)::+::clear(self)
            }
        }
    };
}

macro_rules! impl_acmap_add_assign {
    ($($seg:ident)::+) => {
        impl<'a, S, V, State> crate::traits::ACMapAddAssign<S, V> for $($seg)::+<PauliWord<S>, V, State>
        where
            S: PauliStorage + 'a,
            V: Coefficient + 'a,
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
    ($($seg:ident)::+) => {
        impl<'a, S, V, State> crate::traits::ACMapMulAssign<V> for $($seg)::+<PauliWord<S>, V, State>
        where
            S: PauliStorage + 'a,
            V: Coefficient + 'a,
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
    ($($seg:ident)::+) => {
        impl<'a, S, V, State> crate::traits::ACMapIter<'a> for $($seg)::+<PauliWord<S>, V, State>
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

        // impl<'a, S, V, State> crate::traits::ACMapIterMut<'a, S, V>
        //     for $name<PauliWord<S>, V, State>
        // where
        //     S: PauliStorage + 'a,
        //     V: Coefficient + 'a,
        //     State: Clone + BuildHasher + 'a,
        // {
        //     type Item = (&'a PauliWord<S>, &'a mut V);
        //     type IterMut = std::collections::hash_map::IterMut<'a, PauliWord<S>, V>;
        //     fn iter_mut(&'a mut self) -> Self::IterMut {
        //         Self::iter_mut(self)
        //     }
        // }
    };
}

macro_rules! impl_acmap_trace {
    ($($seg:ident)::+) => {
        impl<'a, P, S, C, State> crate::traits::Trace<'a, P> for $($seg)::+<PauliWord<S>, C, State>
        where
            P: crate::traits::Trace<'a, PauliWord<S>, Output = bool> + 'a,
            S: PauliStorage + 'a,
            C: Coefficient + 'a,
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

macro_rules! impl_acmap_combine_unique {
    ($($seg:ident)::+) => {
        impl<S, C> ACMapCombineUnique for $($seg)::+<PauliWord<S>, C>
        where
            S: PauliStorage,
            C: Coefficient,
        {
            fn combine_unique(&mut self, dest: &mut Self) {
                for (k, v) in self.drain() {
                    dest.insert(k, v);
                }
            }
        }
    };
}

macro_rules! impl_acmap_insert {
    ($($seg:ident)::+) => {
        impl<'a, S, C, State> ACMapInsert<S, C> for $($seg)::+<PauliWord<S>, C, State>
        where
            S: PauliStorage + 'a,
            C: Coefficient + 'a,
            State: Clone + BuildHasher + 'a,
        {
            fn map_insert<F>(&mut self, dest: &mut Self, f: F)
            where
                F: Fn(&PauliWord<S>, &mut C) -> Option<(PauliWord<S>, C)> + Sync + Send,
            {
                for (k, v) in self.iter_mut() {
                    if let Some((new_k, new_v)) = f(k, v) {
                        dest.insert(new_k, new_v);
                    }
                }
            }
        }
    };
}

macro_rules! impl_acmap {
    ($($seg:ident)::+) => {
        impl_acmap_base!($($seg)::+);
        impl_acmap_add_assign!($($seg)::+);
        impl_acmap_mul_assign!($($seg)::+);
        impl_acmap_iter!($($seg)::+);
        impl_acmap_trace!($($seg)::+);
        impl_acmap_combine_unique!($($seg)::+);
        impl_acmap_insert!($($seg)::+);
    };
}

impl_acmap!(HashMap);

#[cfg(feature = "indexmap")]
mod indexmap_impl {
    use super::*;

    impl<'a, S, V, State> crate::traits::ACMapIter<'a> for indexmap::IndexMap<PauliWord<S>, V, State>
    where
        S: PauliStorage + 'a,
        V: Coefficient + 'a,
        State: Clone + BuildHasher + 'a,
    {
        type Item = (&'a PauliWord<S>, &'a V);
        type Iter = indexmap::map::Iter<'a, PauliWord<S>, V>;

        fn iter(&'a self) -> Self::Iter {
            Self::iter(self)
        }
    }

    impl<S, C> ACMapCombineUnique for indexmap::IndexMap<PauliWord<S>, C>
    where
        S: PauliStorage,
        C: Coefficient,
    {
        fn combine_unique(&mut self, dest: &mut Self) {
            for (k, v) in self.drain(..) {
                dest.insert(k, v);
            }
        }
    }

    impl_acmap_base!(indexmap::IndexMap);
    impl_acmap_add_assign!(indexmap::IndexMap);
    impl_acmap_mul_assign!(indexmap::IndexMap);
    impl_acmap_trace!(indexmap::IndexMap);
    impl_acmap_insert!(indexmap::IndexMap);
}

#[cfg(feature = "ahash")]
mod ahash_impl {
    use super::*;

    impl<'a, S, V, State> crate::traits::ACMapBase for ahash::AHashMap<PauliWord<S>, V, State>
    where
        S: PauliStorage,
        V: Coefficient,
        State: Clone + BuildHasher + Default,
    {
        fn with_capacity(capacity: usize) -> Self {
            Self::with_capacity_and_hasher(capacity, State::default())
        }

        fn len(&self) -> usize {
            HashMap::len(self)
        }

        fn clear(&mut self) {
            HashMap::clear(self)
        }
    }

    impl<'a, S, V, State> crate::traits::ACMapIter<'a> for ahash::AHashMap<PauliWord<S>, V, State>
    where
        S: PauliStorage + 'a,
        V: Coefficient + 'a,
        State: Clone + BuildHasher + 'a,
    {
        type Item = (&'a PauliWord<S>, &'a V);
        type Iter = std::collections::hash_map::Iter<'a, PauliWord<S>, V>;

        fn iter(&'a self) -> Self::Iter {
            HashMap::iter(self)
        }
    }

    impl_acmap_add_assign!(ahash::AHashMap);
    impl_acmap_mul_assign!(ahash::AHashMap);
    impl_acmap_trace!(ahash::AHashMap);
    impl_acmap_combine_unique!(ahash::AHashMap);
    impl_acmap_insert!(ahash::AHashMap);
}
