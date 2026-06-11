// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::collections::HashMap;

/// Implement `ACMapBase` for a map type like `HashMap<W, V, Hasher>`.
#[macro_export]
#[doc(hidden)]
macro_rules! impl_acmap_base {
    ($($seg:ident)::+) => {
        impl<'a, V, Hasher, W> $crate::traits::ACMapBase for $($seg)::+<W, V, Hasher>
        where
            V: $crate::traits::Coefficient,
            Hasher: Clone + std::hash::BuildHasher + Default,
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

/// Implement `ACMapAddAssign` for a map type.
#[macro_export]
#[doc(hidden)]
macro_rules! impl_acmap_add_assign {
    ($($seg:ident)::+) => {
        impl<'a, S, V, Hasher, W> $crate::traits::ACMapAddAssign<S, V, Hasher, W> for $($seg)::+<W, V, Hasher>
        where
            S: $crate::traits::PauliStorage + 'a,
            V: $crate::traits::Coefficient + 'a,
            Hasher: Default + Clone + std::hash::BuildHasher + 'a,
            W: $crate::traits::PauliWordTrait + 'a,
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
                    <Self as $crate::traits::ACMapAddAssign<S, V, Hasher, W>>::add_assign(dest, new_k, new_v);
                }
            }
        }
    };
}

/// Implement `ACMapMulAssign` for a map type.
#[macro_export]
#[doc(hidden)]
macro_rules! impl_acmap_mul_assign {
    ($($seg:ident)::+) => {
        impl<'a, V, Hasher, W> $crate::traits::ACMapMulAssign<V, Hasher> for $($seg)::+<W, V, Hasher>
        where
            V: $crate::traits::Coefficient + 'a,
            Hasher: Default + Clone + std::hash::BuildHasher + 'a,
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

/// Implement `ACMapIter` for a map type using std `hash_map::Iter` (for
/// `HashMap` / `AHashMap`). IndexMap supplies its own impl.
#[macro_export]
#[doc(hidden)]
macro_rules! impl_acmap_iter {
    ($($seg:ident)::+) => {
        impl<'a, V, Hasher, W> $crate::traits::ACMapIter<'a> for $($seg)::+<W, V, Hasher>
        where
            V: $crate::traits::Coefficient + 'a,
            Hasher: Default + Clone + std::hash::BuildHasher + 'a,
            W: 'a,
        {
            type Item = (&'a W, &'a V);
            type Iter = std::collections::hash_map::Iter<'a, W, V>;

            fn iter(&'a self) -> Self::Iter {
                Self::iter(self)
            }
        }
    };
}

/// Implement `Trace` for a map type.
#[macro_export]
#[doc(hidden)]
macro_rules! impl_acmap_trace {
    ($($seg:ident)::+) => {
        impl<'a, P, C, Hasher, W> $crate::traits::Trace<'a, P> for $($seg)::+<W, C, Hasher>
        where
            P: $crate::traits::Trace<'a, W, Output = bool> + 'a,
            C: $crate::traits::Coefficient + 'a,
            Hasher: Default + Clone + std::hash::BuildHasher + 'a,
            W: 'a,
        {
            type Output = C;
            fn trace(&'a self, value: &'a P) -> Self::Output {
                let sum = <C as num::Zero>::zero();
                self.iter().fold(sum, |mut acc, (k, v)| {
                    value.trace(k).then(|| acc += v.clone());
                    acc
                })
            }
        }
    };
}

/// Implement `ACMapConsume` for a map type whose `drain()` returns owning
/// `(K, V)` tuples (HashMap-like). IndexMap needs a `drain(..)` argument
/// and is implemented locally.
#[macro_export]
#[doc(hidden)]
macro_rules! impl_acmap_combine_unique {
    ($($seg:ident)::+) => {
        impl<C, Hasher, W> $crate::traits::ACMapConsume for $($seg)::+<W, C, Hasher>
        where
            C: $crate::traits::Coefficient,
            Hasher: Default + Clone + std::hash::BuildHasher,
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

/// Implement `ACMapInsert` for a map type.
#[macro_export]
#[doc(hidden)]
macro_rules! impl_acmap_insert {
    ($($seg:ident)::+) => {
        impl<'a, S, C, Hasher, W> $crate::traits::ACMapInsert<S, C, Hasher, W> for $($seg)::+<W, C, Hasher>
        where
            S: $crate::traits::PauliStorage + 'a,
            C: $crate::traits::Coefficient + 'a,
            Hasher: Default + Clone + std::hash::BuildHasher + 'a,
            W: $crate::traits::PauliWordTrait + 'a,
        {
            fn map_insert<F>(&mut self, dest: &mut Self, f: F)
            where
                F: Fn(&W, &mut C) -> Option<(W, C)> + Sync + Send,
            {
                for (k, v) in self.iter_mut() {
                    if let Some((new_k, new_v)) = f(k, v) {
                        dest.entry(new_k)
                            .and_modify(|val| *val += new_v.clone())
                            .or_insert(new_v);
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

/// Implement `ACMapContains` for a map type.
#[macro_export]
#[doc(hidden)]
macro_rules! impl_acmap_contains {
    ($($seg:ident)::+) => {
        impl<S, C, H, W> $crate::traits::ACMapContains<S, C, H, W> for $($seg)::+<W, C, H>
        where
            S: $crate::traits::PauliStorage,
            C: $crate::traits::Coefficient,
            H: Default + Clone + std::hash::BuildHasher,
            W: $crate::traits::PauliWordTrait,
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

/// Implement `ACMapScale` for a map type.
#[macro_export]
#[doc(hidden)]
macro_rules! impl_acmap_scale {
    ($($seg:ident)::+) => {
        impl<S, V, H, W> $crate::traits::ACMapScale<S, V, H, W> for $($seg)::+<W, V, H>
        where
            S: $crate::traits::PauliStorage,
            V: $crate::traits::Coefficient,
            H: std::hash::BuildHasher + Clone + Default,
            W: $crate::traits::PauliWordTrait,
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

/// Implement `ACMapRetain` for a map type whose `retain` takes `(&K, &mut V) -> bool`
/// (HashMap). AHashMap has the same signature; IndexMap differs in a way the
/// macro handles uniformly.
#[macro_export]
#[doc(hidden)]
macro_rules! impl_acmap_retain {
    ($($seg:ident)::+) => {
        impl<S, V, H, W> $crate::traits::ACMapRetain<S, V, H, W> for $($seg)::+<W, V, H>
        where
            S: $crate::traits::PauliStorage,
            V: $crate::traits::Coefficient,
            H: std::hash::BuildHasher + Clone + Default,
            W: $crate::traits::PauliWordTrait,
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

// Apply all impls to std `HashMap`.
impl_acmap_base!(HashMap);
impl_acmap_add_assign!(HashMap);
impl_acmap_mul_assign!(HashMap);
impl_acmap_iter!(HashMap);
impl_acmap_trace!(HashMap);
impl_acmap_combine_unique!(HashMap);
impl_acmap_insert!(HashMap);
impl_acmap_contains!(HashMap);
impl_acmap_scale!(HashMap);
impl_acmap_retain!(HashMap);
