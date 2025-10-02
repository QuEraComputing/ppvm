use std::{collections::HashMap, hash::BuildHasher};

use crate::{
    traits::{ACMapConsumeUnique, ACMapContains, ACMapInsert, Coefficient, PauliStorage},
    word::PauliWord,
};

macro_rules! impl_acmap_base {
    ($($seg:ident)::+) => {
        impl<'a, S, V, Hasher> crate::traits::ACMapBase for $($seg)::+<PauliWord<S, Hasher>, V, Hasher>
        where
            S: PauliStorage,
            V: Coefficient,
            Hasher: Clone + BuildHasher + Default,
        {
            fn with_capacity(capacity: usize) -> Self {
                Self::with_capacity_and_hasher(capacity, Hasher::default())
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
        impl<'a, S, V, Hasher> crate::traits::ACMapAddAssign<S, V, Hasher> for $($seg)::+<PauliWord<S, Hasher>, V, Hasher>
        where
            S: PauliStorage + 'a,
            V: Coefficient + 'a,
            Hasher: Default + Clone + BuildHasher + 'a,
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
        impl<'a, S, V, Hasher> crate::traits::ACMapMulAssign<V, Hasher> for $($seg)::+<PauliWord<S, Hasher>, V, Hasher>
        where
            S: PauliStorage + 'a,
            V: Coefficient + 'a,
            Hasher: Default + Clone + BuildHasher + 'a,
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
        impl<'a, S, V, Hasher> crate::traits::ACMapIter<'a> for $($seg)::+<PauliWord<S, Hasher>, V, Hasher>
        where
            S: PauliStorage + 'a,
            V: Coefficient + 'a,
            Hasher: Default + Clone + BuildHasher + 'a,
        {
            type Item = (&'a PauliWord<S, Hasher>, &'a V);
            type Iter = std::collections::hash_map::Iter<'a, PauliWord<S, Hasher>, V>;

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
        impl<'a, P, S, C, Hasher> crate::traits::Trace<'a, P> for $($seg)::+<PauliWord<S, Hasher>, C, Hasher>
        where
            P: crate::traits::Trace<'a, PauliWord<S, Hasher>, Output = bool> + 'a,
            S: PauliStorage + 'a,
            C: Coefficient + 'a,
            Hasher: Default + Clone + BuildHasher + 'a,
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
        impl<S, C, Hasher> ACMapConsumeUnique for $($seg)::+<PauliWord<S, Hasher>, C, Hasher>
        where
            S: PauliStorage,
            C: Coefficient,
            Hasher: Default + Clone + BuildHasher,
        {
            fn consume_unique(&mut self, dest: &mut Self) {
                for (k, v) in dest.drain() {
                    self.entry(k)
                        .and_modify(|val| *val += v.clone())
                        .or_insert(v);
                }
            }
        }
    };
}

macro_rules! impl_acmap_insert {
    ($($seg:ident)::+) => {
        impl<'a, S, C, Hasher> ACMapInsert<S, C, Hasher> for $($seg)::+<PauliWord<S, Hasher>, C, Hasher>
        where
            S: PauliStorage + 'a,
            C: Coefficient + 'a,
            Hasher: Default + Clone + BuildHasher + 'a,
        {
            fn map_insert<F>(&mut self, dest: &mut Self, f: F)
            where
                F: Fn(&PauliWord<S, Hasher>, &mut C) -> Option<(PauliWord<S, Hasher>, C)> + Sync + Send,
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

macro_rules! impl_acmap_contains {
    ($($seg:ident)::+) => {
        impl<S, C, H> ACMapContains<S, C, H> for $($seg)::+<PauliWord<S, H>, C, H>
        where
            S: PauliStorage,
            C: Coefficient,
            H: Default + Clone + BuildHasher,
        {
            fn contains_with<F>(&self, key: &PauliWord<S, H>, f: F) -> bool
            where
                F: Fn(&C) -> bool,
            {
                match self.get(key) {
                    Some(v) => f(v),
                    None => false,
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
        impl_acmap_contains!($($seg)::+);
    };
}

impl_acmap!(HashMap);

#[cfg(feature = "indexmap")]
mod indexmap_impl {
    use super::*;

    impl<'a, S, V, H> crate::traits::ACMapIter<'a> for indexmap::IndexMap<PauliWord<S, H>, V, H>
    where
        S: PauliStorage + 'a,
        V: Coefficient + 'a,
        H: Default + Clone + BuildHasher + 'a,
    {
        type Item = (&'a PauliWord<S, H>, &'a V);
        type Iter = indexmap::map::Iter<'a, PauliWord<S, H>, V>;

        fn iter(&'a self) -> Self::Iter {
            Self::iter(self)
        }
    }

    impl<S, C, H> ACMapConsumeUnique for indexmap::IndexMap<PauliWord<S, H>, C, H>
    where
        S: PauliStorage,
        C: Coefficient,
        H: Default + Clone + BuildHasher,
    {
        fn consume_unique(&mut self, dest: &mut Self) {
            for (k, v) in dest.drain(..) {
                self.entry(k)
                    .and_modify(|val| *val += v.clone())
                    .or_insert(v);
            }
        }
    }

    impl_acmap_base!(indexmap::IndexMap);
    impl_acmap_add_assign!(indexmap::IndexMap);
    impl_acmap_mul_assign!(indexmap::IndexMap);
    impl_acmap_trace!(indexmap::IndexMap);
    impl_acmap_insert!(indexmap::IndexMap);
    impl_acmap_contains!(indexmap::IndexMap);
}

#[cfg(feature = "ahash")]
mod ahash_impl {
    use super::*;

    impl<S, V, H> crate::traits::ACMapBase for ahash::AHashMap<PauliWord<S, H>, V, H>
    where
        S: PauliStorage,
        V: Coefficient,
        H: Clone + BuildHasher + Default,
    {
        fn with_capacity(capacity: usize) -> Self {
            Self::with_capacity_and_hasher(capacity, H::default())
        }

        fn len(&self) -> usize {
            HashMap::len(self)
        }

        fn clear(&mut self) {
            HashMap::clear(self)
        }
    }

    impl<'a, S, V, H> crate::traits::ACMapIter<'a> for ahash::AHashMap<PauliWord<S, H>, V, H>
    where
        S: PauliStorage + 'a,
        V: Coefficient + 'a,
        H: Default + Clone + BuildHasher + 'a,
    {
        type Item = (&'a PauliWord<S, H>, &'a V);
        type Iter = std::collections::hash_map::Iter<'a, PauliWord<S, H>, V>;

        fn iter(&'a self) -> Self::Iter {
            HashMap::iter(self)
        }
    }

    impl_acmap_add_assign!(ahash::AHashMap);
    impl_acmap_mul_assign!(ahash::AHashMap);
    impl_acmap_trace!(ahash::AHashMap);
    impl_acmap_combine_unique!(ahash::AHashMap);
    impl_acmap_insert!(ahash::AHashMap);
    impl_acmap_contains!(ahash::AHashMap);
}
