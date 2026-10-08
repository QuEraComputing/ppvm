// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Stores Pauli words in separate X/Z plane arrays with a shared qubit width.
//! Each plane entry holds one word’s packed storage; `ColumnStore` manages hash columns.
//! Implements column reading, mutation, and Pauli access; selected via `Columnar`.
//! Uses the shared structural hash so column and individual-word hashes agree.

use std::marker::PhantomData;

use bitvec::array::BitArray;
use ppvm_traits_2::{Columnar, KeyColumn, KeyColumnMut, PauliColumn};

use crate::data::PauliWord;
use crate::hash::structural_hash;
use crate::{HashFinalize, PauliStorage};
use std::hash::BuildHasher;

/// A structure-of-arrays column of [`PauliWord`]s: parallel X and Z plane blocks
/// plus the shared qubit width.
pub struct PauliKeyColumn<A: PauliStorage, H = fxhash::FxBuildHasher> {
    xplanes: Vec<A>,
    zplanes: Vec<A>,
    nqubits: usize,
    _hasher: PhantomData<fn() -> H>,
}

impl<A, H> PauliKeyColumn<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    /// Flips one selected bit directly, without constructing or hashing a word.
    #[inline(always)]
    fn toggle_plane(plane: &mut A, qubit: usize, toggle: bool) {
        if !toggle {
            return;
        }
        #[cfg(target_endian = "little")]
        {
            bytemuck::bytes_of_mut(plane)[qubit >> 3] ^= 1 << (qubit & 7);
        }
        #[cfg(target_endian = "big")]
        {
            let mut bits = BitArray::new(*plane);
            let value = bits[qubit];
            bits.set(qubit, !value);
            *plane = bits.data;
        }
    }

    #[inline(always)]
    fn plane_bit(plane: &A, qubit: usize) -> bool {
        #[cfg(target_endian = "little")]
        {
            let bytes = bytemuck::bytes_of(plane);
            bytes[qubit >> 3] & (1 << (qubit & 7)) != 0
        }
        #[cfg(target_endian = "big")]
        {
            BitArray::<A>::new(*plane)[qubit]
        }
    }
}

impl<A: PauliStorage, H> Default for PauliKeyColumn<A, H> {
    #[inline]
    fn default() -> Self {
        Self {
            xplanes: Vec::new(),
            zplanes: Vec::new(),
            nqubits: 0,
            _hasher: PhantomData,
        }
    }
}

impl<A: PauliStorage, H> Clone for PauliKeyColumn<A, H> {
    #[inline]
    fn clone(&self) -> Self {
        Self {
            xplanes: self.xplanes.clone(),
            zplanes: self.zplanes.clone(),
            nqubits: self.nqubits,
            _hasher: PhantomData,
        }
    }
}

impl<A, H> KeyColumn for PauliKeyColumn<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    type Key = PauliWord<A, H>;

    #[inline]
    fn len(&self) -> usize {
        self.xplanes.len()
    }

    #[inline]
    fn capacity(&self) -> usize {
        self.xplanes.capacity().min(self.zplanes.capacity())
    }

    /// Hashes each pair of planes so `out[i]` equals `self.get(i).key_hash()`.
    #[inline]
    fn hash_into(&self, out: &mut [u64]) {
        debug_assert_eq!(out.len(), self.len(), "hash column length mismatch");
        for (i, slot) in out.iter_mut().enumerate() {
            *slot = structural_hash::<A, H>(&self.xplanes[i], &self.zplanes[i], self.nqubits);
        }
    }

    #[inline]
    fn key_eq(&self, i: usize, other: &Self::Key) -> bool {
        self.nqubits == other.nqubits
            && self.xplanes[i] == other.xbits.data
            && self.zplanes[i] == other.zbits.data
    }

    #[inline]
    fn gather(&self, indices: &[u32]) -> Self {
        let mut xplanes = Vec::with_capacity(indices.len());
        let mut zplanes = Vec::with_capacity(indices.len());
        for &idx in indices {
            xplanes.push(self.xplanes[idx as usize]);
            zplanes.push(self.zplanes[idx as usize]);
        }
        Self {
            xplanes,
            zplanes,
            nqubits: self.nqubits,
            _hasher: PhantomData,
        }
    }

    #[inline]
    fn get(&self, i: usize) -> Self::Key {
        PauliWord::from_planes(
            BitArray::new(self.xplanes[i]),
            BitArray::new(self.zplanes[i]),
            self.nqubits,
        )
    }
}

impl<A, H> KeyColumnMut for PauliKeyColumn<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn with_capacity(n: usize) -> Self {
        Self {
            xplanes: Vec::with_capacity(n),
            zplanes: Vec::with_capacity(n),
            nqubits: 0,
            _hasher: PhantomData,
        }
    }

    /// Appends a key’s planes. The first key sets the width; later keys must match.
    #[inline]
    fn push(&mut self, key: Self::Key) {
        if self.xplanes.is_empty() {
            self.nqubits = key.nqubits;
        } else {
            assert_eq!(self.nqubits, key.nqubits, "column width mismatch");
        }
        self.xplanes.push(key.xbits.data);
        self.zplanes.push(key.zbits.data);
    }

    /// Reserves room in both plane arrays to reduce reallocations during batch insertion.
    #[inline]
    fn reserve(&mut self, additional: usize) {
        self.xplanes.reserve(additional);
        self.zplanes.reserve(additional);
    }

    /// Clears both plane arrays while retaining their allocations and stored width.
    /// The next push sets the width again because the column is empty.
    #[inline]
    fn clear(&mut self) {
        self.xplanes.clear();
        self.zplanes.clear();
    }

    /// Replaces one key’s planes without moving other entries or reallocating.
    /// Used by `ColumnStore` for in-place Clifford updates.
    #[inline]
    fn set(&mut self, i: usize, key: Self::Key) {
        assert_eq!(self.nqubits, key.nqubits, "column width mismatch");
        self.xplanes[i] = key.xbits.data;
        self.zplanes[i] = key.zbits.data;
    }

    #[inline]
    fn truncate(&mut self, len: usize) {
        self.xplanes.truncate(len);
        self.zplanes.truncate(len);
    }

    #[inline]
    fn swap_remove(&mut self, i: usize) -> Self::Key {
        let x = self.xplanes.swap_remove(i);
        let z = self.zplanes.swap_remove(i);
        PauliWord::from_planes(BitArray::new(x), BitArray::new(z), self.nqubits)
    }
}

impl<A, H> PauliColumn for PauliKeyColumn<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline(always)]
    fn x_bit(&self, row: usize, qubit: usize) -> bool {
        Self::plane_bit(&self.xplanes[row], qubit)
    }

    #[inline(always)]
    fn z_bit(&self, row: usize, qubit: usize) -> bool {
        Self::plane_bit(&self.zplanes[row], qubit)
    }

    #[inline(always)]
    fn toggled_bits(&self, row: usize, qubit: usize, toggle_x: bool, toggle_z: bool) -> Self::Key {
        assert!(qubit < self.nqubits, "qubit {qubit} out of bounds");
        let mut x = self.xplanes[row];
        let mut z = self.zplanes[row];
        Self::toggle_plane(&mut x, qubit, toggle_x);
        Self::toggle_plane(&mut z, qubit, toggle_z);
        PauliWord::from_planes(BitArray::new(x), BitArray::new(z), self.nqubits)
    }

    /// Copies the planes once, applies both masks, and hashes only the final word.
    #[inline(always)]
    fn toggled_bits2(
        &self,
        row: usize,
        i: usize,
        [xi, zi]: [bool; 2],
        j: usize,
        [xj, zj]: [bool; 2],
    ) -> Self::Key {
        assert!(i < self.nqubits, "qubit {i} out of bounds");
        assert!(j < self.nqubits, "qubit {j} out of bounds");
        let mut x = self.xplanes[row];
        let mut z = self.zplanes[row];
        Self::toggle_plane(&mut x, i, xi);
        Self::toggle_plane(&mut z, i, zi);
        Self::toggle_plane(&mut x, j, xj);
        Self::toggle_plane(&mut z, j, zj);
        PauliWord::from_planes(BitArray::new(x), BitArray::new(z), self.nqubits)
    }
}

impl<A, H> Columnar for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    type Column = PauliKeyColumn<A, H>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_traits_2::{Indexable, PauliBits};

    #[test]
    fn mismatched_width_mutations_leave_column_unchanged() {
        use std::panic::{AssertUnwindSafe, catch_unwind};

        for (original, mismatched) in [("I", "II"), ("II", "I")] {
            let mut col = column(&[original]);
            let expected = PauliWord::from(original);
            let capacity = col.capacity();

            assert!(
                catch_unwind(AssertUnwindSafe(|| {
                    col.push(PauliWord::from(mismatched));
                }))
                .is_err()
            );
            assert_eq!(col.len(), 1);
            assert_eq!(col.get(0), expected);
            assert_eq!(col.capacity(), capacity);

            assert!(
                catch_unwind(AssertUnwindSafe(|| {
                    col.set(0, PauliWord::from(mismatched));
                }))
                .is_err()
            );
            assert_eq!(col.len(), 1);
            assert_eq!(col.get(0), expected);
            assert_eq!(col.capacity(), capacity);

            let mut hashes = [0];
            col.hash_into(&mut hashes);
            assert_eq!(hashes[0], expected.key_hash());
        }
    }

    #[test]
    fn toggles_reject_indices_beyond_logical_width() {
        let col = column(&["I"]);
        for invalid in [1, 7, 63] {
            for mask in 0..4 {
                let toggle = [mask & 1 != 0, mask & 2 != 0];
                assert!(
                    std::panic::catch_unwind(|| {
                        col.toggled_bits(0, invalid, toggle[0], toggle[1])
                    })
                    .is_err()
                );
                assert!(
                    std::panic::catch_unwind(|| {
                        col.toggled_bits2(0, invalid, toggle, 0, [false; 2])
                    })
                    .is_err()
                );
                assert!(
                    std::panic::catch_unwind(|| {
                        col.toggled_bits2(0, 0, [false; 2], invalid, toggle)
                    })
                    .is_err()
                );
            }
        }
        assert_eq!(col.get(0), PauliWord::from("I"));
    }

    fn check_toggles<A: PauliStorage>() {
        let word = PauliWord::<A>::from("XYZI".repeat(32).as_str());
        let mut col = PauliKeyColumn::<A>::default();
        col.push(word);
        for i in 0..128 {
            assert_eq!(col.x_bit(0, i), word.x_bit(i));
            assert_eq!(col.z_bit(0, i), word.z_bit(i));
            for mask in 0..4 {
                let (x, z) = (mask & 1 != 0, mask & 2 != 0);
                let actual = col.toggled_bits(0, i, x, z);
                let expected = word.toggled_bits(i, x, z);
                assert_eq!(actual, expected);
                assert_eq!(actual.key_hash(), expected.key_hash());
            }
        }
        for (i, j) in [(0, 0), (2, 5), (7, 8), (63, 64), (127, 0)] {
            for mask in 0..16 {
                let left = [mask & 1 != 0, mask & 2 != 0];
                let right = [mask & 4 != 0, mask & 8 != 0];
                let actual = col.toggled_bits2(0, i, left, j, right);
                let expected = word.toggled_bits2(i, left, j, right);
                assert_eq!(actual, expected);
                assert_eq!(actual.key_hash(), expected.key_hash());
            }
        }
        assert_eq!(col.get(0), word);
    }

    #[test]
    fn packed_toggles_match_words_across_storage_widths() {
        check_toggles::<[u8; 16]>();
        check_toggles::<[usize; 4]>();
    }

    #[test]
    fn mutations_preserve_planes_hashes_and_capacity() {
        let mut col = column(&["XYZI", "IIIZ", "YYYY"]);
        col.reserve(8);
        let capacity = col.capacity();
        col.set(1, PauliWord::from("ZZZZ"));
        assert_eq!(col.swap_remove(0), PauliWord::from("XYZI"));
        assert_eq!(col.get(0), PauliWord::from("YYYY"));
        assert_eq!(col.get(1), PauliWord::from("ZZZZ"));
        let picked = col.gather(&[1, 0, 1]);
        let mut hashes = vec![0; picked.len()];
        picked.hash_into(&mut hashes);
        for (i, hash) in hashes.into_iter().enumerate() {
            assert_eq!(hash, picked.get(i).key_hash());
        }
        col.truncate(1);
        assert_eq!(col.len(), 1);
        col.clear();
        assert!(col.is_empty());
        assert_eq!(col.capacity(), capacity);
        col.push(PauliWord::from("XY"));
        assert_eq!(col.get(0), PauliWord::from("XY"));
    }

    fn column(words: &[&str]) -> PauliKeyColumn<u64> {
        let mut col = PauliKeyColumn::<u64>::with_capacity(words.len());
        for w in words {
            col.push(PauliWord::from(*w));
        }
        col
    }

    #[test]
    fn roundtrip_and_len() {
        let words = ["XYZI", "IIIZ", "YYYY"];
        let col = column(&words);
        assert_eq!(col.len(), 3);
        assert!(!col.is_empty());
        for (i, w) in words.iter().enumerate() {
            assert_eq!(col.get(i), PauliWord::from(*w));
        }
    }

    #[test]
    fn hash_into_matches_scalar_key_hash() {
        let words = ["XYZI", "IIIZ", "YYYY", "ZXZX"];
        let col = column(&words);
        let mut out = vec![0u64; col.len()];
        col.hash_into(&mut out);
        for (i, w) in words.iter().enumerate() {
            assert_eq!(out[i], PauliWord::<u64>::from(*w).key_hash(), "key {i}");
        }
    }

    #[test]
    fn key_eq_and_gather() {
        let col = column(&["XYZI", "IIIZ", "YYYY"]);
        assert!(col.key_eq(1, &PauliWord::from("IIIZ")));
        assert!(!col.key_eq(1, &PauliWord::from("XYZI")));

        let picked = col.gather(&[2, 0]);
        assert_eq!(picked.len(), 2);
        assert_eq!(picked.get(0), PauliWord::from("YYYY"));
        assert_eq!(picked.get(1), PauliWord::from("XYZI"));
    }
}
