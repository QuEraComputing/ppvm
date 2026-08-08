// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Insertion-ordered sparse storage for ABI surfaces whose term order is
//! observable.
//!
//! The four buffers mirror [`crate::HashMapStore`].  Only `primary` is value
//! state; the other buffers are persistent workspaces and are empty between
//! operations.  IndexMap replacement keeps an existing key's position, while a
//! first insertion appends it.  The gate implementations preserve the ordering
//! of the legacy `config::indexmap` backend.

mod algebra;
mod branching;
mod gates;
mod lifecycle;

use indexmap::IndexMap;
use ppvm_traits_2::{IdentityBuildHasher, TermBatch};

/// An insertion-ordered hash-join backend with persistent gate workspaces.
#[derive(Debug)]
pub struct IndexMapStore<K, C> {
    pub(super) primary: IndexMap<K, C, IdentityBuildHasher>,
    pub(super) aux: IndexMap<K, C, IdentityBuildHasher>,
    pub(super) scratch: Vec<(K, C)>,
    pub(super) batch: TermBatch<K, C>,
}

impl<K: Clone + Eq + std::hash::Hash, C: Clone> Clone for IndexMapStore<K, C> {
    fn clone(&self) -> Self {
        let mut primary =
            IndexMap::with_capacity_and_hasher(self.primary.capacity(), IdentityBuildHasher);
        primary.extend(self.primary.iter().map(|(k, c)| (k.clone(), c.clone())));
        Self {
            primary,
            aux: IndexMap::with_capacity_and_hasher(self.aux.capacity(), IdentityBuildHasher),
            scratch: Vec::with_capacity(self.scratch.capacity()),
            batch: TermBatch::with_capacity(self.batch.capacity()),
        }
    }
}

impl<K, C> PartialEq for IndexMapStore<K, C>
where
    K: Eq + std::hash::Hash,
    C: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.primary == other.primary
    }
}

#[cfg(test)]
mod tests;
