// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Collision-safe basis lookup storing positions rather than duplicate words.

use crate::word::{Word, word_hash};
use hashbrown::{HashTable, hash_table::Entry};

pub(crate) struct CompactBasisIndex<'a, const C: usize> {
    basis: &'a [Word<C>],
    table: HashTable<u32>,
}

impl<'a, const C: usize> CompactBasisIndex<'a, C> {
    pub(crate) fn new(basis: &'a [Word<C>]) -> Self {
        assert!(
            basis.len() < u32::MAX as usize,
            "basis exceeds u32 positions"
        );
        let mut table = HashTable::with_capacity(basis.len());
        for (position, word) in basis.iter().enumerate() {
            match table.entry(
                word_hash(word),
                |&index| basis[index as usize] == *word,
                |&index| word_hash(&basis[index as usize]),
            ) {
                Entry::Vacant(entry) => {
                    entry.insert(position as u32);
                }
                Entry::Occupied(_) => panic!("duplicate Pauli word in basis"),
            }
        }
        Self { basis, table }
    }

    /// [`Self::get`] for a word whose [`word_hash`] the caller already has.
    pub(crate) fn get_hashed(&self, word: &Word<C>, hash: u64) -> Option<u32> {
        self.table
            .find(hash, |&index| self.basis[index as usize] == *word)
            .copied()
    }

    #[cfg(test)]
    pub(crate) fn get(&self, word: &Word<C>) -> Option<u32> {
        // Control-byte fingerprints filter probes before loading full borrowed words.
        self.table
            .find(word_hash(word), |&index| {
                self.basis[index as usize] == *word
            })
            .copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::word_from_codes;

    #[test]
    fn compact_lookup_checks_words_and_wide_registers() {
        let basis: Vec<_> = (0..192)
            .map(|site| {
                let mut codes = vec![0; 128];
                codes[site / 3] = (site % 3 + 1) as u8;
                word_from_codes::<2>(&codes).unwrap()
            })
            .collect();
        let index = CompactBasisIndex::new(&basis);
        for (i, word) in basis.iter().enumerate() {
            assert_eq!(index.get(word), Some(i as u32));
        }
        assert_eq!(index.get(&word_from_codes::<2>(&[2; 128]).unwrap()), None);
        let wide = vec![word_from_codes::<8>(&vec![1; 512]).unwrap()];
        assert_eq!(CompactBasisIndex::new(&wide).get(&wide[0]), Some(0));
    }
}
