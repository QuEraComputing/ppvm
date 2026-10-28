// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Collision-safe basis lookup storing positions rather than duplicate words.

use crate::word::{Word, word_hash};

pub(crate) struct CompactBasisIndex<'a, const C: usize> {
    basis: &'a [Word<C>],
    slots: Vec<u32>,
}

impl<'a, const C: usize> CompactBasisIndex<'a, C> {
    pub(crate) fn new(basis: &'a [Word<C>]) -> Self {
        assert!(
            basis.len() < u32::MAX as usize,
            "basis exceeds u32 positions"
        );
        // Keep an empty probe terminator with at most three quarters of slots occupied.
        let capacity = basis
            .len()
            .saturating_mul(4)
            .div_ceil(3)
            .max(2)
            .next_power_of_two();
        let mut index = Self {
            basis,
            slots: vec![0; capacity],
        };
        for (position, word) in basis.iter().enumerate() {
            let mut slot = word_hash(word) as usize & (capacity - 1);
            while index.slots[slot] != 0 {
                let other = index.slots[slot] as usize - 1;
                assert_ne!(basis[other], *word, "duplicate Pauli word in basis");
                slot = (slot + 1) & (capacity - 1);
            }
            index.slots[slot] = position as u32 + 1;
        }
        index
    }

    pub(crate) fn get(&self, word: &Word<C>) -> Option<u32> {
        let mask = self.slots.len() - 1;
        let mut slot = word_hash(word) as usize & mask;
        loop {
            let entry = self.slots[slot];
            if entry == 0 {
                return None;
            }
            if self.basis[entry as usize - 1] == *word {
                return Some(entry - 1);
            }
            slot = (slot + 1) & mask;
        }
    }
}
