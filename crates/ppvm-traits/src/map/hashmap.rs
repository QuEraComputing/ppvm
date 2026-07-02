// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::{collections::HashMap, hash::BuildHasher};

use crate::traits::*;

/// Below this many source entries the rayon fan-out costs more than it saves,
/// so `map_insert` / `map_insert_vec` stay on the sequential loop. Only used
/// when the `rayon` acceleration path is compiled in.
#[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
const PARALLEL_MAP_INSERT_THRESHOLD: usize = 1024;

macro_rules! impl_acmap_base {
    ($($seg:ident)::+) => {
        impl<'a, V, Hasher, W> crate::traits::ACMapBase for $($seg)::+<W, V, Hasher>
        where
            V: Coefficient,
            Hasher: Clone + BuildHasher + Default,
            W:,
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
        impl<'a, S, V, Hasher, W> crate::traits::ACMapAddAssign<S, V, Hasher, W> for $($seg)::+<W, V, Hasher>
        where
            S: PauliStorage + 'a,
            V: Coefficient + 'a,
            Hasher: Default + Clone + BuildHasher + 'a,
            W: PauliWordTrait + 'a,
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
                for (k, v) in self.iter() {
                    let (new_k, new_v) = f(k, v);
                    <Self as ACMapAddAssign<S, V, Hasher, W>>::add_assign(dest, new_k, new_v);
                }
            }
        }
    };
}

macro_rules! impl_acmap_mul_assign {
    ($($seg:ident)::+) => {
        impl<'a, V, Hasher, W> crate::traits::ACMapMulAssign<V, Hasher> for $($seg)::+<W, V, Hasher>
        where
            V: Coefficient + 'a,
            Hasher: Default + Clone + BuildHasher + 'a,
            W:,
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
        impl<'a, V, Hasher, W> crate::traits::ACMapIter<'a> for $($seg)::+<W, V, Hasher>
        where
            V: Coefficient + 'a,
            Hasher: Default + Clone + BuildHasher + 'a,
            W: 'a,
        {
            type Item = (&'a W, &'a V);
            type Iter = std::collections::hash_map::Iter<'a, W, V>;

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
        impl<'a, P, C, Hasher, W> crate::traits::Trace<'a, P> for $($seg)::+<W, C, Hasher>
        where
            P: crate::traits::Trace<'a, W, Output = bool> + 'a,
            C: Coefficient + 'a,
            Hasher: Default + Clone + BuildHasher + 'a,
            W: 'a,
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
        impl<C, Hasher, W> ACMapConsume for $($seg)::+<W, C, Hasher>
        where
            C: Coefficient,
            Hasher: Default + Clone + BuildHasher,
            W: std::hash::Hash + std::cmp::Eq,
        {
            fn consume(&mut self, dest: &mut Self) {
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
        impl<'a, S, C, Hasher, W> ACMapInsert<S, C, Hasher, W> for $($seg)::+<W, C, Hasher>
        where
            S: PauliStorage + 'a,
            C: Coefficient + 'a,
            Hasher: Default + Clone + BuildHasher + 'a,
            // `Sync + Send` on the word lets the rayon path below shuttle
            // `(&W, &mut C)` entry handles across worker threads. `C` is
            // already `Sync + Send` via `Coefficient`.
            W: PauliWordTrait + Sync + Send + 'a,
        {
            fn map_insert<F>(&mut self, dest: &mut Self, f: F)
            where
                F: Fn(&W, &mut C) -> Option<(W, C)> + Sync + Send,
            {
                // Parallel path: scale surviving values in place and collect
                // anti-commuting branch terms concurrently, then merge those
                // terms into `dest` in one sequential pass (collision keys
                // must sum). The merge is cheap relative to the closure work
                // and needs the map's `entry` API, so it stays sequential.
                #[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
                if self.len() >= PARALLEL_MAP_INSERT_THRESHOLD {
                    use rayon::prelude::*;
                    // Collect the disjoint `&mut` value handles once so rayon
                    // can split them; each entry's `*v *= cos` mutation is
                    // independent, so this is data-race free.
                    let entries: Vec<(&W, &mut C)> = self.iter_mut().collect();
                    let produced: Vec<(W, C)> = entries
                        .into_par_iter()
                        .filter_map(|(k, v)| f(k, v))
                        .collect();
                    for (new_k, new_v) in produced {
                        dest.entry(new_k)
                            .and_modify(|val| *val += new_v.clone())
                            .or_insert(new_v);
                    }
                    return;
                }
                for (k, v) in self.iter_mut() {
                    if let Some((new_k, new_v)) = f(k, v) {
                        dest.entry(new_k)
                            .and_modify(|val| *val += new_v.clone())
                            .or_insert(new_v);
                    }
                }
            }

            fn map_insert_vec<F>(&mut self, dest: &mut Vec<(W, C)>, f: F)
            where
                F: Fn(&W, &mut C) -> Option<(W, C)> + Sync + Send,
            {
                // Hot rot2 path (via `PauliSum::map_insert`). Same structure as
                // `map_insert`, but produced terms go into a plain `Vec` that
                // the caller later folds into its map. rayon's `par_extend`
                // accumulates into per-thread buffers and concatenates them.
                #[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
                if self.len() >= PARALLEL_MAP_INSERT_THRESHOLD {
                    use rayon::prelude::*;
                    let entries: Vec<(&W, &mut C)> = self.iter_mut().collect();
                    dest.par_extend(entries.into_par_iter().filter_map(|(k, v)| f(k, v)));
                    return;
                }
                for (k, v) in self.iter_mut() {
                    if let Some(entry) = f(k, v) {
                        dest.push(entry);
                    }
                }
            }

            fn map_insert_multiple<F>(&mut self, dest: &mut Self, f: F)
            where
                F: Fn(&W, &mut C) -> Option<Vec<(W, C)>> + Sync + Send,
            {
                for (k, v) in self.iter_mut() {
                    if let Some(entries) = f(k, v) {
                        for (new_k, new_v) in entries {
                            dest.entry(new_k)
                                .and_modify(|val| *val += new_v.clone())
                                .or_insert(new_v);
                        }
                    }
                }
            }
        }
    };
}

macro_rules! impl_acmap_contains {
    ($($seg:ident)::+) => {
        impl<S, C, H, W> ACMapContains<S, C, H, W> for $($seg)::+<W, C, H>
        where
            S: PauliStorage,
            C: Coefficient,
            H: Default + Clone + BuildHasher,
            W: PauliWordTrait,
        {
            fn contains_with<F>(&self, key: &W, f: F) -> bool
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

macro_rules! impl_acmap_scale {
    ($($seg:ident)::+) => {
        impl<S, V, H, W> ACMapScale<S, V, H, W> for $($seg)::+<W, V, H>
        where
            S: PauliStorage,
            V: Coefficient,
            H: BuildHasher + Clone + Default,
            W: PauliWordTrait,
        {
            fn scale<F>(&mut self, f: F)
            where
                F: Fn(&W, &mut V) + Sync + Send,
            {
                for (k, v) in self.iter_mut() {
                    f(k, v);
                }
            }
        }
    };
}

macro_rules! impl_acmap_retain {
    ($($seg:ident)::+) => {
        impl<S, V, H, W> ACMapRetain<S, V, H, W> for $($seg)::+<W, V, H>
        where
            S: PauliStorage,
            V: Coefficient,
            H: BuildHasher + Clone + Default,
            W: PauliWordTrait,
        {
            fn retain<F>(&mut self, mut f: F)
            where
                F: FnMut(&W, &V) -> bool,
            {
                Self::retain(self, |k, v| f(k, v));
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
        impl_acmap_scale!($($seg)::+);
        impl_acmap_retain!($($seg)::+);
    };
}

impl_acmap!(HashMap);

#[cfg(feature = "indexmap")]
mod indexmap_impl {
    use super::*;

    impl<'a, V, H, W> crate::traits::ACMapIter<'a> for indexmap::IndexMap<W, V, H>
    where
        V: Coefficient + 'a,
        H: Default + Clone + BuildHasher + 'a,
        W: 'a,
    {
        type Item = (&'a W, &'a V);
        type Iter = indexmap::map::Iter<'a, W, V>;

        fn iter(&'a self) -> Self::Iter {
            Self::iter(self)
        }
    }

    impl<C, H, W> ACMapConsume for indexmap::IndexMap<W, C, H>
    where
        C: Coefficient,
        H: Default + Clone + BuildHasher,
        W: std::hash::Hash + std::cmp::Eq,
    {
        fn consume(&mut self, dest: &mut Self) {
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
    impl_acmap_scale!(indexmap::IndexMap);
    impl_acmap_retain!(indexmap::IndexMap);
}

#[cfg(all(feature = "ahash", not(target_arch = "wasm32")))]
mod ahash_impl {
    use super::*;

    impl<V, H, W> crate::traits::ACMapBase for ahash::AHashMap<W, V, H>
    where
        V: Coefficient,
        H: Clone + BuildHasher + Default,
        W:,
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

    impl<'a, V, H, W> crate::traits::ACMapIter<'a> for ahash::AHashMap<W, V, H>
    where
        V: Coefficient + 'a,
        H: Default + Clone + BuildHasher + 'a,
        W: 'a,
    {
        type Item = (&'a W, &'a V);
        type Iter = std::collections::hash_map::Iter<'a, W, V>;

        fn iter(&'a self) -> Self::Iter {
            HashMap::iter(self)
        }
    }

    impl<S, V, H, W> ACMapRetain<S, V, H, W> for ahash::AHashMap<W, V, H>
    where
        S: PauliStorage,
        V: Coefficient,
        H: BuildHasher + Clone + Default,
        W: PauliWordTrait,
    {
        fn retain<F>(&mut self, mut f: F)
        where
            F: FnMut(&W, &V) -> bool,
        {
            HashMap::retain(self, |k, v| f(k, v));
        }
    }

    impl_acmap_add_assign!(ahash::AHashMap);
    impl_acmap_mul_assign!(ahash::AHashMap);
    impl_acmap_trace!(ahash::AHashMap);
    impl_acmap_combine_unique!(ahash::AHashMap);
    impl_acmap_insert!(ahash::AHashMap);
    impl_acmap_contains!(ahash::AHashMap);
    impl_acmap_scale!(ahash::AHashMap);
}
