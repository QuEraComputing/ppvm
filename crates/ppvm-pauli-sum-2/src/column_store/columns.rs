// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;
use hashbrown::HashTable;

pub(super) struct Columns<K: Columnar, C> {
    pub(super) keys: K::Column,
    pub(super) coeffs: Vec<C>,
    pub(super) hashes: Vec<u64>,
    pub(super) live: Vec<u8>,
    pub(super) live_len: usize,
    pub(super) sparse_rows: Vec<u32>,
    pub(super) live_runs: Vec<(u32, u32)>,
    pub(super) sparse_cache_dirty: bool,
    pub(super) index: HashTable<(u32, u64)>,
}

impl<K: Columnar, C> Columns<K, C> {
    #[inline]
    pub(super) fn rows(&self) -> usize {
        self.coeffs.len()
    }

    #[inline(always)]
    pub(super) fn is_live(&self, row: usize) -> bool {
        self.live[row] != 0
    }

    #[inline]
    pub(super) fn is_dense(&self) -> bool {
        self.live_len == self.rows()
    }
}

impl<K, C> Columns<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    pub(super) fn with_capacity(cap: usize) -> Self {
        Self {
            keys: K::Column::with_capacity(cap),
            coeffs: Vec::with_capacity(cap),
            hashes: Vec::with_capacity(cap),
            live: Vec::with_capacity(cap),
            live_len: 0,
            sparse_rows: Vec::new(),
            live_runs: Vec::new(),
            sparse_cache_dirty: false,
            index: HashTable::with_capacity(cap),
        }
    }

    #[inline]
    pub(super) fn len(&self) -> usize {
        self.live_len
    }

    pub(super) fn clear(&mut self) {
        self.keys.clear();
        self.coeffs.clear();
        self.hashes.clear();
        self.live.clear();
        self.live_len = 0;
        self.sparse_rows.clear();
        self.live_runs.clear();
        self.sparse_cache_dirty = false;
        self.index.clear();
    }

    #[inline(always)]
    pub(super) fn find_any(&self, key: &K, hash: u64) -> Option<usize> {
        self.index
            .find(hash, |&(slot, stored_hash)| {
                let slot = slot as usize;
                stored_hash == hash && self.keys.key_eq(slot, key)
            })
            .map(|&(slot, _)| slot as usize)
    }

    #[inline(always)]
    pub(super) fn find(&self, key: &K, hash: u64) -> Option<usize> {
        self.find_any(key, hash).filter(|&row| self.is_live(row))
    }

    pub(super) fn reindex(&mut self) {
        self.index.clear();
        self.index.reserve(self.live_len, |&(_, hash)| hash);
        for slot in 0..self.rows() {
            if !self.is_live(slot) {
                continue;
            }
            let hash = self.hashes[slot];
            #[cfg(debug_assertions)]
            {
                let key = self.key(slot);
                debug_assert!(
                    self.find(&key, hash).is_none(),
                    "ColumnStore re-key produced duplicate keys"
                );
            }
            self.index
                .insert_unique(hash, (slot as u32, hash), |&(_, stored_hash)| stored_hash);
        }
    }

    #[inline]
    pub(super) fn add(&mut self, key: K, coeff: C) {
        let hash = key.key_hash();
        self.add_prehashed::<true>(key, hash, coeff);
    }

    #[inline(always)]
    pub(super) fn add_dense(&mut self, key: K, coeff: C) {
        let hash = key.key_hash();
        self.add_prehashed::<false>(key, hash, coeff);
    }

    #[inline]
    fn add_prehashed<const TOMBSTONES: bool>(&mut self, key: K, hash: u64, coeff: C) {
        let vacant_slot = self.rows();
        match self.index.entry(
            hash,
            |&(slot, stored_hash)| {
                let slot = slot as usize;
                stored_hash == hash && self.keys.key_eq(slot, &key)
            },
            |&(_, stored_hash)| stored_hash,
        ) {
            hashbrown::hash_table::Entry::Occupied(mut entry) => {
                let row = entry.get().0 as usize;
                if TOMBSTONES && self.live[row] == 0 {
                    debug_assert!(
                        u32::try_from(vacant_slot).is_ok(),
                        "ColumnStore slot overflow"
                    );
                    self.keys.push(key);
                    self.hashes.push(hash);
                    self.coeffs.push(coeff);
                    self.live.push(1);
                    self.live_len += 1;
                    if !self.sparse_cache_dirty {
                        self.sparse_rows.push(vacant_slot as u32);
                        Self::push_live_run_row(&mut self.live_runs, vacant_slot);
                    }
                    entry.get_mut().0 = vacant_slot as u32;
                } else {
                    self.coeffs[row] += coeff;
                }
            }
            hashbrown::hash_table::Entry::Vacant(entry) => {
                debug_assert!(
                    u32::try_from(vacant_slot).is_ok(),
                    "ColumnStore slot overflow"
                );
                self.keys.push(key);
                self.hashes.push(hash);
                self.coeffs.push(coeff);
                self.live.push(1);
                self.live_len += 1;
                if !self.sparse_cache_dirty && !self.sparse_rows.is_empty() {
                    self.sparse_rows.push(vacant_slot as u32);
                    Self::push_live_run_row(&mut self.live_runs, vacant_slot);
                }
                entry.insert((vacant_slot as u32, hash));
            }
        }
        self.debug_assert_valid();
    }

    /// Merge a branch that is expected to hit the closed support. This keeps
    /// hashbrown's insertion/growth machinery off the occupied fast path; a
    /// genuine miss falls back to the general accumulating insert.
    #[inline(always)]
    pub(super) fn add_likely_present(&mut self, key: K, coeff: C) {
        let hash = key.key_hash();
        if let Some(slot) = self.find(&key, hash) {
            self.coeffs[slot] += coeff;
        } else {
            self.add_prehashed::<true>(key, hash, coeff);
        }
    }

    /// Dense-only counterpart of [`Self::add_likely_present`]. The caller
    /// guarantees every physical row is live, so the occupied probe needs no
    /// tombstone filter.
    #[inline(always)]
    pub(super) fn add_likely_present_dense(&mut self, key: K, coeff: C) {
        let hash = key.key_hash();
        if let Some(slot) = self.find_any(&key, hash) {
            self.coeffs[slot] += coeff;
        } else {
            self.add_prehashed::<false>(key, hash, coeff);
        }
    }

    #[inline]
    pub(super) fn insert(&mut self, key: K, coeff: C) {
        let hash = key.key_hash();
        let vacant_slot = self.rows();
        match self.index.entry(
            hash,
            |&(slot, stored_hash)| {
                let slot = slot as usize;
                stored_hash == hash && self.keys.key_eq(slot, &key)
            },
            |&(_, stored_hash)| stored_hash,
        ) {
            hashbrown::hash_table::Entry::Occupied(mut entry) => {
                let row = entry.get().0 as usize;
                if self.live[row] != 0 {
                    self.coeffs[row] = coeff;
                } else {
                    debug_assert!(
                        u32::try_from(vacant_slot).is_ok(),
                        "ColumnStore slot overflow"
                    );
                    self.keys.push(key);
                    self.hashes.push(hash);
                    self.coeffs.push(coeff);
                    self.live.push(1);
                    self.live_len += 1;
                    if !self.sparse_cache_dirty {
                        self.sparse_rows.push(vacant_slot as u32);
                        Self::push_live_run_row(&mut self.live_runs, vacant_slot);
                    }
                    entry.get_mut().0 = vacant_slot as u32;
                }
            }
            hashbrown::hash_table::Entry::Vacant(entry) => {
                debug_assert!(
                    u32::try_from(vacant_slot).is_ok(),
                    "ColumnStore slot overflow"
                );
                self.keys.push(key);
                self.hashes.push(hash);
                self.coeffs.push(coeff);
                self.live.push(1);
                self.live_len += 1;
                if !self.sparse_cache_dirty && !self.sparse_rows.is_empty() {
                    self.sparse_rows.push(vacant_slot as u32);
                    Self::push_live_run_row(&mut self.live_runs, vacant_slot);
                }
                entry.insert((vacant_slot as u32, hash));
            }
        }
        self.debug_assert_valid();
    }

    #[inline]
    pub(super) fn key(&self, i: usize) -> K {
        self.keys.get(i)
    }

    pub(super) fn retain_coeffs(&mut self, keep: impl Fn(&C) -> bool) {
        let mut changed = false;
        for row in 0..self.rows() {
            if self.is_live(row) && !keep(&self.coeffs[row]) {
                self.live[row] = 0;
                self.live_len -= 1;
                changed = true;
            }
        }
        self.sparse_cache_dirty |= changed;
        self.compact_if_needed();
    }

    pub(super) fn retain_terms(&mut self, keep: impl Fn(&K, &C) -> bool) {
        let mut changed = false;
        for row in 0..self.rows() {
            if self.is_live(row) {
                let key = self.keys.get(row);
                if !keep(&key, &self.coeffs[row]) {
                    self.live[row] = 0;
                    self.live_len -= 1;
                    changed = true;
                }
            }
        }
        self.sparse_cache_dirty |= changed;
        self.compact_if_needed();
    }

    fn compact_if_needed(&mut self) {
        let rows = self.rows();
        let dead = rows - self.live_len;
        if dead != 0 && dead >= rows.div_ceil(8) {
            self.compact();
        }
        self.debug_assert_valid();
    }

    fn compact(&mut self) {
        let rows = self.rows();
        let mut write = 0;
        for read in 0..rows {
            if !self.is_live(read) {
                continue;
            }
            if write != read {
                self.keys.set(write, self.keys.get(read));
                self.hashes[write] = self.hashes[read];
                self.coeffs.swap(write, read);
            }
            self.live[write] = 1;
            write += 1;
        }
        debug_assert_eq!(write, self.live_len);
        self.keys.truncate(write);
        self.hashes.truncate(write);
        self.coeffs.truncate(write);
        self.live.truncate(write);
        self.sparse_rows.clear();
        self.live_runs.clear();
        self.sparse_cache_dirty = false;
        self.reindex();
    }

    pub(super) fn ensure_sparse_cache(&mut self) {
        if !self.sparse_cache_dirty {
            return;
        }
        self.sparse_rows.clear();
        self.sparse_rows.reserve(self.live_len);
        self.live_runs.clear();
        for row in 0..self.rows() {
            if self.is_live(row) {
                self.sparse_rows.push(row as u32);
                Self::push_live_run_row(&mut self.live_runs, row);
            }
        }
        self.sparse_cache_dirty = false;
    }

    fn push_live_run_row(live_runs: &mut Vec<(u32, u32)>, row: usize) {
        let row = row as u32;
        if let Some((_, end)) = live_runs.last_mut()
            && *end == row
        {
            *end += 1;
        } else {
            live_runs.push((row, row + 1));
        }
    }

    pub(super) fn reserve_for_live_len(&mut self, target: usize) {
        let additional = target.saturating_sub(self.live_len);
        self.coeffs.reserve(additional);
        self.hashes.reserve(additional);
        self.live.reserve(additional);
        if !self.sparse_cache_dirty && !self.sparse_rows.is_empty() {
            self.sparse_rows.reserve(additional);
        }
        self.keys.reserve(additional);
        self.index.reserve(additional, |&(_, hash)| hash);
    }

    #[inline]
    pub(super) fn debug_assert_valid(&self) {
        debug_assert_eq!(self.keys.len(), self.rows());
        debug_assert_eq!(self.hashes.len(), self.rows());
        debug_assert_eq!(self.live.len(), self.rows());
        debug_assert_eq!(
            self.live.iter().filter(|&&state| state != 0).count(),
            self.live_len
        );
        debug_assert!(if self.is_dense() {
            !self.sparse_cache_dirty && self.sparse_rows.is_empty() && self.live_runs.is_empty()
        } else if self.sparse_cache_dirty {
            true
        } else {
            self.sparse_rows.len() == self.live_len
                && self
                    .live_runs
                    .iter()
                    .map(|&(start, end)| (end - start) as usize)
                    .sum::<usize>()
                    == self.live_len
                && self
                    .sparse_rows
                    .iter()
                    .all(|&row| self.is_live(row as usize))
        });
    }
}
