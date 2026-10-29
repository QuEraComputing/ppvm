// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Stable candidate positions and chunk-capped ranking metadata.

use crate::candidates::CandidateBuffer;
use crate::word::{Word, word_hash};
use hashbrown::HashTable;

#[derive(Default)]
pub(crate) struct IndexedCandidates {
    pub(crate) state: Vec<u8>,
    table: HashTable<u32>,
    pub(crate) free: Vec<u32>,
}

impl IndexedCandidates {
    pub(crate) fn len(&self) -> usize {
        self.table.len()
    }

    pub(crate) fn reset_reusing_table(&mut self) {
        self.state.clear();
        self.free.clear();
        self.table.clear();
    }

    pub(crate) fn reset(&mut self) {
        self.state.clear();
        self.free.clear();
        // Table growth history determines admission order, including magnitude ties.
        self.table = HashTable::new();
    }

    pub(crate) fn add<const C: usize>(
        &mut self,
        buffer: &mut CandidateBuffer<'_, C>,
        word: Word<C>,
        value: f64,
    ) {
        let hash = word_hash(&word);
        if let Some(&index) = self.table.find(hash, |&index| {
            buffer.words[buffer.prefix + index as usize] == word
        }) {
            buffer.numeric[index as usize] += value;
            return;
        }
        let index = if let Some(index) = self.free.pop() {
            buffer.words[buffer.prefix + index as usize] = word;
            buffer.numeric[index as usize] = 0.0 + value;
            self.state[index as usize] = 1;
            index
        } else {
            if buffer.words.len() == buffer.words.capacity() {
                buffer
                    .words
                    .reserve_exact((buffer.words.capacity() / 4).max(32));
                buffer
                    .numeric
                    .reserve_exact(buffer.words.capacity() - buffer.prefix - buffer.numeric.len());
                self.state
                    .reserve_exact(buffer.words.capacity() - buffer.prefix - self.state.len());
            }
            let index = u32::try_from(buffer.numeric.len()).expect("candidate count exceeds u32");
            buffer.words.push(word);
            buffer.numeric.push(0.0 + value);
            self.state.push(1);
            index
        };
        self.table.insert_unique(hash, index, |&index| {
            word_hash(&buffer.words[buffer.prefix + index as usize])
        });
    }

    pub(crate) fn cap<const C: usize>(&mut self, buffer: &mut CandidateBuffer<'_, C>, room: usize) {
        if self.table.len() <= room {
            return;
        }
        if room == 0 {
            buffer.words.truncate(buffer.prefix);
            buffer.numeric.clear();
            self.state.clear();
            self.table.clear();
            self.free.clear();
            return;
        }
        let count = buffer.numeric.len();
        buffer.numeric.reserve_exact(self.table.len());
        // The rate prefix remains live while the tail holds nonnegative ranking keys.
        for index in 0..count {
            if self.state[index] != 0 {
                buffer.numeric.push(buffer.numeric[index].abs());
            }
        }
        let (rates, keys) = buffer.numeric.split_at_mut(count);
        let position = keys.len() - room;
        let cutoff = keys
            .select_nth_unstable_by_key(position, |v| v.to_bits())
            .1
            .to_bits();
        let above = keys.iter().filter(|v| v.to_bits() > cutoff).count();
        let mut ties = room - above;
        for (rate, state) in rates.iter().zip(&mut self.state) {
            if *state != 0 {
                let magnitude = rate.abs().to_bits();
                *state = if magnitude > cutoff {
                    3
                } else if magnitude == cutoff {
                    2
                } else {
                    1
                };
            }
        }
        self.table.retain(|index| {
            let state = &mut self.state[*index as usize];
            if *state == 3 {
                return true;
            }
            let keep = *state == 2 && ties > 0;
            ties -= keep as usize;
            if !keep {
                *state = 0;
                self.free.push(*index);
            }
            keep
        });
        buffer.numeric.truncate(count);
    }

    pub(crate) fn selected_rates<'a, const C: usize>(
        &'a self,
        buffer: &'a CandidateBuffer<'_, C>,
        threshold: f64,
    ) -> impl Iterator<Item = (Word<C>, f64)> + 'a {
        self.table.iter().filter_map(move |&index| {
            let index = index as usize;
            (buffer.numeric[index].abs() > threshold)
                .then_some((buffer.words[buffer.prefix + index], buffer.numeric[index]))
        })
    }

    pub(crate) fn select<const C: usize>(
        &mut self,
        buffer: &mut CandidateBuffer<'_, C>,
        threshold: f64,
    ) -> usize {
        self.table
            .retain(|index| buffer.numeric[*index as usize].abs() > threshold);
        let slots = buffer.numeric.len();
        buffer.numeric.clear();
        buffer.numeric.resize(2 * slots, 0.0);
        let (positions, labels) = buffer.numeric.split_at_mut(slots);
        // Integer positions use lossless bit patterns in the shared numeric buffer.
        // No floating arithmetic reads these temporary permutation entries.
        for index in 0..slots {
            let code = f64::from_bits(index as u64);
            positions[index] = code;
            labels[index] = code;
        }
        let mut count = 0;
        for &original in &self.table {
            if self.state[original as usize] == 0 {
                continue;
            }
            let source = positions[original as usize].to_bits() as usize;
            buffer
                .words
                .swap(buffer.prefix + count, buffer.prefix + source);
            let displaced = labels[count].to_bits() as usize;
            positions[displaced] = f64::from_bits(source as u64);
            labels[source] = f64::from_bits(displaced as u64);
            positions[original as usize] = f64::from_bits(count as u64);
            labels[count] = f64::from_bits(original as u64);
            count += 1;
        }
        count
    }
}
