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

    pub(crate) fn get(&self, word: &Word<C>) -> Option<u32> {
        // Control-byte fingerprints filter probes before loading full borrowed words.
        self.table
            .find(word_hash(word), |&index| {
                self.basis[index as usize] == *word
            })
            .copied()
    }
}
