// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Candidate payload stored once, with a Swiss-table index of u32 positions.

use crate::word::{Word, word_hash};
use hashbrown::HashTable;

pub(crate) struct IndexedCandidates<const C: usize> {
    entries: Vec<(Word<C>, f64)>,
    table: HashTable<u32>,
    free: Vec<u32>,
}

impl<const C: usize> IndexedCandidates<C> {
    pub(crate) fn new() -> Self {
        Self {
            entries: Vec::new(),
            table: HashTable::new(),
            free: Vec::new(),
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
        let index = if let Some(index) = self.free.pop() {
            self.entries[index as usize] = (word, 0.0 + value);
            index
        } else {
            if self.entries.len() == self.entries.capacity() {
                self.entries
                    .reserve_exact((self.entries.capacity() / 4).max(32));
            }
            let index = u32::try_from(self.entries.len()).expect("candidate count exceeds u32");
            self.entries.push((word, 0.0 + value));
            index
        };
        self.table.insert_unique(hash, index, |&index| {
            word_hash(&self.entries[index as usize].0)
        });
    }

    pub(crate) fn cap(&mut self, room: usize, bits: &mut Vec<u64>) {
        if self.table.len() <= room {
            return;
        }
        if room == 0 {
            self.entries.clear();
            self.table.clear();
            self.free.clear();
            return;
        }
        bits.clear();
        if bits.capacity() < self.table.len() {
            bits.reserve_exact(self.table.len());
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
        self.table.retain(|index| {
            let magnitude = self.entries[*index as usize].1.abs().to_bits();
            if magnitude > cutoff {
                return true;
            }
            let keep = magnitude == cutoff && ties > 0;
            ties -= keep as usize;
            if !keep {
                self.free.push(*index);
            }
            keep
        });
    }

    pub(crate) fn into_vec(self, threshold: f64) -> Vec<(Word<C>, f64)> {
        let Self {
            mut entries,
            table,
            free,
        } = self;
        let live = table.len();
        let mut order = Vec::with_capacity(entries.len());
        order.extend(table.iter().copied());
        drop(table);
        order.extend(free);
        assert_eq!(order.len(), entries.len());
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
        entries.truncate(live);
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
