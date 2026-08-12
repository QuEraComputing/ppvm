// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::columns::Columns;
use super::*;

/// The structure-of-arrays storage backend for [`Sum`](crate::Sum): parallel key
/// / coefficient / digest columns behind an open-addressed index, plus the same
/// persistent workspace [`HashMapStore`](crate::HashMapStore) carries.
///
/// Observationally identical to [`HashMapStore`](crate::HashMapStore) — see the
/// module docs — so `Sum<ColumnStore<K, C>, P>` is a drop-in swap. The domain
/// alias is [`ColumnPauliSum`](crate::ColumnPauliSum).
pub struct ColumnStore<K: Columnar, C> {
    pub(super) primary: Columns<K, C>,
    pub(super) aux: Columns<K, C>,
    pub(super) scratch: Vec<(K, C)>,
    pub(super) visited: Vec<bool>,
    pub(super) batch: TermBatch<K, C>,
}

impl<K, C> Clone for Columns<K, C>
where
    K: Columnar,
    C: Clone,
{
    fn clone(&self) -> Self {
        let mut keys = K::Column::with_capacity(self.keys.capacity());
        for i in 0..self.rows() {
            keys.push(self.keys.get(i));
        }
        let mut coeffs = Vec::with_capacity(self.coeffs.capacity());
        coeffs.extend(self.coeffs.iter().cloned());
        let mut hashes = Vec::with_capacity(self.hashes.capacity());
        hashes.extend_from_slice(&self.hashes);
        let mut live = Vec::with_capacity(self.live.capacity());
        live.extend_from_slice(&self.live);
        Self {
            keys,
            coeffs,
            hashes,
            live,
            live_len: self.live_len,
            sparse_rows: self.sparse_rows.clone(),
            live_runs: self.live_runs.clone(),
            sparse_cache_dirty: self.sparse_cache_dirty,
            index: self.index.clone(),
        }
    }
}

impl<K, C> ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    #[inline]
    pub(super) fn support_len(&self) -> usize {
        self.primary.len()
    }
}

impl<K, C> Clone for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    fn clone(&self) -> Self {
        Self {
            primary: self.primary.clone(),
            aux: Columns::with_capacity(
                self.aux
                    .coeffs
                    .capacity()
                    .max(self.primary.coeffs.capacity()),
            ),
            scratch: Vec::with_capacity(self.scratch.capacity()),
            visited: Vec::with_capacity(self.visited.capacity()),
            batch: TermBatch::with_capacity(self.batch.capacity()),
        }
    }
}

impl<K, C> std::fmt::Debug for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.primary.debug_assert_valid();
        f.debug_struct("ColumnStore")
            .field("len", &self.support_len())
            .field("rows", &self.primary.rows())
            .field("index_capacity", &self.primary.index.capacity())
            .finish()
    }
}

impl<K, C> StoreAlloc for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    fn with_capacity(cap: usize) -> Self {
        Self {
            primary: Columns::with_capacity(cap),
            // `aux` is multiplication-only and `scratch` is branching-only.
            // Keep their capacity lazy so constructing a sum does not allocate
            // and zero a second full index before either path is used.
            aux: Columns::with_capacity(0),
            scratch: Vec::new(),
            visited: Vec::new(),
            batch: TermBatch::new(),
        }
    }

    #[inline]
    fn reset(&mut self) {
        self.primary.clear();
    }
}

impl<K, C> AddTerm<K, C> for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    #[inline]
    fn add_term(&mut self, key: K, coeff: C) {
        self.primary.add(key, coeff);
    }
}

impl<K, C> InsertTerm<K, C> for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    #[inline]
    fn insert_term(&mut self, key: K, coeff: C) {
        self.primary.insert(key, coeff);
    }
}

impl<K, C> ApplyProducer<K, C> for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    /// Produce into the persistent batch, then replace the primary support.
    fn apply_producer<TP>(&mut self, producer: TP)
    where
        TP: TermProducer<K, C>,
    {
        self.batch.clear();
        for i in 0..self.primary.rows() {
            if !self.primary.is_live(i) {
                continue;
            }
            let key = self.primary.key(i);
            let coeff = self.primary.coeffs[i].clone();
            producer.produce(&key, &coeff, &mut self.batch);
        }
        self.primary.clear();
        for (k, c) in self.batch.iter() {
            self.primary.add(k.clone(), c.clone());
        }
        self.batch.clear();
    }
}

/// Equality observes the primary support only and is insertion-order independent.
impl<K, C> PartialEq for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    fn eq(&self, other: &Self) -> bool {
        if self.primary.len() != other.primary.len() {
            return false;
        }
        (0..self.primary.rows()).all(|i| {
            if !self.primary.is_live(i) {
                return true;
            }
            let key = self.primary.key(i);
            match other.primary.find(&key, self.primary.hashes[i]) {
                Some(slot) => other.primary.coeffs[slot] == self.primary.coeffs[i],
                None => false,
            }
        })
    }
}
