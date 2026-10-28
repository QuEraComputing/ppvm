// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Candidate payload stored once, with a Swiss-table index of u32 positions.

use crate::word::{Word, word_hash};
use hashbrown::HashTable;

pub(crate) struct IndexedCandidates<const C: usize> {
    entries: Vec<(Word<C>, f64)>,
    table: HashTable<u32>,
    removed: Vec<u32>,
    removed_bits: Vec<u64>,
}

impl<const C: usize> IndexedCandidates<C> {
    pub(crate) fn new() -> Self {
        Self {
            entries: Vec::new(),
            table: HashTable::new(),
            removed: Vec::new(),
            removed_bits: Vec::new(),
        }
    }

    pub(crate) fn add(&mut self, word: Word<C>, value: f64) {
        let hash = word_hash(&word);
        if let Some(&index) = self
            .table
            .find(hash, |&index| self.entries[index as usize].0 == word)
        {
            self.entries[index as usize].1 += value;
            return;
        }
        if self.entries.len() == self.entries.capacity() {
            self.entries
                .reserve_exact((self.entries.capacity() / 4).max(32));
        }
        let index = u32::try_from(self.entries.len()).expect("candidate count exceeds u32");
        self.entries.push((word, value));
        self.table.insert_unique(hash, index, |&index| {
            word_hash(&self.entries[index as usize].0)
        });
    }

    pub(crate) fn cap(&mut self, room: usize, bits: &mut Vec<u64>) {
        if self.entries.len() <= room {
            return;
        }
        if room == 0 {
            self.entries.clear();
            self.table.clear();
            return;
        }
        bits.clear();
        if bits.capacity() < self.entries.len() {
            bits.reserve_exact(self.entries.len());
        }
        bits.extend(
            self.table
                .iter()
                .map(|&index| self.entries[index as usize].1.abs().to_bits()),
        );
        let position = bits.len() - room;
        let cutoff = *bits.select_nth_unstable(position).1;
        let above = bits.iter().filter(|&&value| value > cutoff).count();
        let mut ties = room - above;
        self.removed.clear();
        self.removed_bits.clear();
        self.removed_bits.resize(self.entries.len().div_ceil(64), 0);
        self.table.retain(|index| {
            let magnitude = self.entries[*index as usize].1.abs().to_bits();
            if magnitude > cutoff {
                return true;
            }
            let keep = magnitude == cutoff && ties > 0;
            ties -= keep as usize;
            if !keep {
                self.removed.push(*index);
                self.removed_bits[*index as usize / 64] |= 1 << (*index as usize % 64);
            }
            keep
        });
        for &index in &self.removed {
            while let Some(last) = self.entries.len().checked_sub(1) {
                if self.removed_bits[last / 64] & (1 << (last % 64)) == 0 {
                    break;
                }
                self.entries.pop();
            }
            if index as usize >= self.entries.len() {
                continue;
            }
            if self.removed_bits[index as usize / 64] & (1 << (index as usize % 64)) == 0 {
                continue;
            }
            let last = self.entries.len() - 1;
            self.entries.swap_remove(index as usize);
            self.removed_bits[index as usize / 64] &= !(1 << (index as usize % 64));
            if index as usize != last {
                let word = &self.entries[index as usize].0;
                *self
                    .table
                    .find_mut(word_hash(word), |&candidate| candidate as usize == last)
                    .expect("moved candidate must remain indexed") = index;
            }
        }
    }

    pub(crate) fn into_vec(self, threshold: f64) -> Vec<(Word<C>, f64)> {
        let Self {
            mut entries, table, ..
        } = self;
        let mut order: Vec<u32> = table.iter().copied().collect();
        drop(table);
        // Apply the table's visitation permutation without a second payload vector.
        for start in 0..order.len() {
            if order[start] as usize == start {
                continue;
            }
            let saved = entries[start];
            let mut destination = start;
            loop {
                let source = order[destination] as usize;
                order[destination] = destination as u32;
                if source == start {
                    entries[destination] = saved;
                    break;
                }
                entries[destination] = entries[source];
                destination = source;
            }
        }
        entries.retain(|(_, value)| value.abs() > threshold);
        entries
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::truncate::cap_map_to_room;
    use crate::word_from_codes;
    use fxhash::FxHashMap;

    #[test]
    fn indexed_payload_preserves_reference_caps_and_visitation() {
        let mut reference = FxHashMap::default();
        let mut indexed = IndexedCandidates::new();
        for chunk in 0..12 {
            for i in 0..500 {
                let code = (i * 37 + chunk * 11) % 4096;
                let word = word_from_codes::<2>(
                    &(0..6)
                        .map(|q| ((code >> (2 * q)) & 3) as u8)
                        .collect::<Vec<_>>(),
                )
                .unwrap();
                let value = ((i + chunk) % 9) as f64 - 4.0;
                *reference.entry(word).or_insert(0.0) += value;
                indexed.add(word, value);
            }
            cap_map_to_room(&mut reference, 300);
            indexed.cap(300, &mut Vec::new());
        }
        let expected: Vec<_> = reference
            .into_iter()
            .filter(|(_, value)| value.abs() > 0.0)
            .collect();
        assert_eq!(indexed.into_vec(0.0), expected);
    }
}
