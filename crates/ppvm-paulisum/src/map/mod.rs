// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! ACMap impls for optional map backends (DashMap, IndexMap, AHashMap).
//!
//! The std `HashMap` impl lives in `ppvm_runtime::map`; the helper macros
//! are re-used from there via `$crate`-qualified paths.

#[cfg(feature = "dashmap")]
mod dashmap;

#[cfg(feature = "indexmap")]
mod indexmap_impl {
    use ppvm_runtime::traits::*;
    use ppvm_runtime::{
        impl_acmap_add_assign, impl_acmap_base, impl_acmap_contains, impl_acmap_insert,
        impl_acmap_mul_assign, impl_acmap_retain, impl_acmap_scale, impl_acmap_trace,
    };
    use std::hash::BuildHasher;

    impl<'a, V, H, W> ACMapIter<'a> for indexmap::IndexMap<W, V, H>
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

#[cfg(feature = "ahash")]
mod ahash_impl {
    use ppvm_runtime::traits::*;
    use ppvm_runtime::{
        impl_acmap_add_assign, impl_acmap_combine_unique, impl_acmap_contains, impl_acmap_insert,
        impl_acmap_mul_assign, impl_acmap_scale, impl_acmap_trace,
    };
    use std::collections::HashMap;
    use std::hash::BuildHasher;

    impl<V, H, W> ACMapBase for ahash::AHashMap<W, V, H>
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

    impl<'a, V, H, W> ACMapIter<'a> for ahash::AHashMap<W, V, H>
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
