// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Split candidate words and rates with stable Swiss-table position indices.

use crate::word::{Word, word_hash};
use hashbrown::HashTable;

pub(crate) struct IndexedCandidates<const C: usize> {
    words: Vec<Word<C>>,
    rates: Vec<f64>,
    state: Vec<u8>,
    table: HashTable<u32>,
    free: Vec<u32>,
    prefix: usize,
    tail: bool,
}

impl<const C: usize> IndexedCandidates<C> {
    pub(crate) fn new() -> Self {
        Self {
            words: Vec::new(),
            rates: Vec::new(),
            state: Vec::new(),
            table: HashTable::new(),
            free: Vec::new(),
            prefix: 0,
            tail: false,
        }
    }

    pub(crate) fn allocation_bytes(&self) -> [(&'static str, usize); 4] {
        [
            (
                "candidate_words",
                self.words.capacity() * std::mem::size_of::<Word<C>>(),
            ),
            ("candidate_rates", self.rates.capacity() * 8),
            ("candidate_state", self.state.capacity()),
            ("candidate_free", self.free.capacity() * 4),
        ]
    }

    pub(crate) fn take_numeric(&mut self) -> Vec<f64> {
        std::mem::take(&mut self.rates)
    }

    pub(crate) fn set_numeric(&mut self, mut numeric: Vec<f64>) {
        numeric.clear();
        self.rates = numeric;
    }

    pub(crate) fn reset(&mut self) {
        self.words.clear();
        self.prefix = 0;
        self.tail = false;
        self.rates.clear();
        self.state.clear();
        self.free.clear();
        // Table growth history determines admission order, including magnitude ties.
        self.table = HashTable::new();
    }

    pub(crate) fn attach_basis(&mut self, basis: Vec<Word<C>>) {
        self.prefix = basis.len();
        self.words = basis;
        self.tail = true;
    }

    pub(crate) fn source_basis(&self) -> &[Word<C>] {
        &self.words[..self.prefix]
    }

    pub(crate) fn add(&mut self, word: Word<C>, value: f64) {
        let hash = word_hash(&word);
        if let Some(&index) = self.table.find(hash, |&index| {
            self.words[self.prefix + index as usize] == word
        }) {
            self.rates[index as usize] += value;
            return;
        }
        let index = if let Some(index) = self.free.pop() {
            self.words[self.prefix + index as usize] = word;
            self.rates[index as usize] = 0.0 + value;
            self.state[index as usize] = 1;
            index
        } else {
            if self.words.len() == self.words.capacity() {
                self.words
                    .reserve_exact((self.words.capacity() / 4).max(32));
                self.rates
                    .reserve_exact(self.words.capacity() - self.prefix - self.rates.len());
                self.state
                    .reserve_exact(self.words.capacity() - self.prefix - self.state.len());
            }
            let index = u32::try_from(self.rates.len()).expect("candidate count exceeds u32");
            self.words.push(word);
            self.rates.push(0.0 + value);
            self.state.push(1);
            index
        };
        self.table.insert_unique(hash, index, |&index| {
            word_hash(&self.words[self.prefix + index as usize])
        });
    }

    pub(crate) fn cap(&mut self, room: usize) {
        if self.table.len() <= room {
            return;
        }
        if room == 0 {
            self.words.truncate(self.prefix);
            self.rates.clear();
            self.state.clear();
            self.table.clear();
            self.free.clear();
            return;
        }
        let count = self.rates.len();
        self.rates.reserve_exact(self.table.len());
        // The rate prefix remains live while the tail holds nonnegative ranking keys.
        for index in 0..count {
            if self.state[index] != 0 {
                self.rates.push(self.rates[index].abs());
            }
        }
        let (rates, keys) = self.rates.split_at_mut(count);
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
        self.rates.truncate(count);
    }

    pub(crate) fn into_vec(self, threshold: f64) -> Vec<(Word<C>, f64)> {
        self.table
            .iter()
            .filter_map(|&index| {
                let index = index as usize;
                (self.rates[index].abs() > threshold)
                    .then_some((self.words[self.prefix + index], self.rates[index]))
            })
            .collect()
    }

    pub(crate) fn admit(
        &mut self,
        basis: &mut Vec<Word<C>>,
        coefficients: &mut Vec<f64>,
        maximum: usize,
        threshold: f64,
    ) {
        let retained = if self.tail { self.prefix } else { basis.len() };
        assert!(self.table.len() <= maximum.saturating_sub(retained));
        if self.tail {
            self.table
                .retain(|index| self.rates[*index as usize].abs() > threshold);
            let slots = self.rates.len();
            self.rates.clear();
            self.rates.resize(2 * slots, 0.0);
            let (positions, labels) = self.rates.split_at_mut(slots);
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
                self.words.swap(self.prefix + count, self.prefix + source);
                let displaced = labels[count].to_bits() as usize;
                positions[displaced] = f64::from_bits(source as u64);
                labels[source] = f64::from_bits(displaced as u64);
                positions[original as usize] = f64::from_bits(count as u64);
                labels[count] = f64::from_bits(original as u64);
                count += 1;
            }
            self.words.truncate(self.prefix + count);
            *basis = std::mem::take(&mut self.words);
            coefficients.resize(basis.len(), 0.0);
        } else {
            for &index in &self.table {
                let index = index as usize;
                if self.rates[index].abs() > threshold {
                    basis.push(self.words[index]);
                    coefficients.push(0.0);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::truncate::cap_map_to_room;
    use crate::word_from_codes;
    use fxhash::FxHashMap;

    #[test]
    fn tail_admission_preserves_nonfinite_threshold_rejection() {
        let seed = word_from_codes::<2>(&[1, 0]).unwrap();
        let candidate = word_from_codes::<2>(&[0, 2]).unwrap();
        let mut tail = IndexedCandidates::new();
        tail.attach_basis(vec![seed]);
        tail.add(candidate, 0.5);
        let mut basis = Vec::new();
        let mut coeffs = vec![1.0];
        tail.admit(&mut basis, &mut coeffs, 10, f64::NAN);
        assert_eq!(basis, vec![seed]);
        assert_eq!(coeffs, vec![1.0]);
    }

    #[test]
    fn indexed_payload_preserves_reference_caps_and_visitation() {
        let mut reference = FxHashMap::default();
        let mut indexed = IndexedCandidates::new();
        let mut tail = IndexedCandidates::new();
        let seed = word_from_codes::<2>(&[1; 6]).unwrap();
        tail.attach_basis(vec![seed]);
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
                tail.add(word, value);
            }
            cap_map_to_room(&mut reference, 300);
            indexed.cap(300);
            tail.cap(300);
        }
        let expected: Vec<_> = reference
            .into_iter()
            .filter(|(_, value)| value.abs() > 0.0)
            .collect();
        let mut basis = Vec::new();
        let mut coefficients = vec![0.5];
        tail.admit(&mut basis, &mut coefficients, 301, 1.5);
        assert_eq!(
            basis,
            std::iter::once(seed)
                .chain(
                    expected
                        .iter()
                        .filter(|(_, v)| v.abs() > 1.5)
                        .map(|(w, _)| *w)
                )
                .collect::<Vec<_>>()
        );
        assert_eq!(coefficients[0], 0.5);
        assert!(coefficients[1..].iter().all(|v| *v == 0.0));
        assert_eq!(indexed.into_vec(0.0), expected);
    }
}
