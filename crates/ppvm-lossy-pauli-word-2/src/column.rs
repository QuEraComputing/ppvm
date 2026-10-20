// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Column storage for [`LossyPauliWord`], with separate X, Z, and loss vectors.
//! Each vector holds one packed `A` value per word; all words share a qubit width.

use std::hash::BuildHasher;
use std::marker::PhantomData;

use bitvec::array::BitArray;
use ppvm_pauli_word_2::{HashFinalize, PauliStorage};
use ppvm_traits_2::{Columnar, KeyColumn, KeyColumnMut, LossColumn, PauliColumn};

use crate::data::LossyPauliWord;
use crate::hash::structural_hash_lossy;

/// Packed X, Z, and loss planes for words with the same qubit width.
pub struct LossyPauliKeyColumn<A: PauliStorage, H = fxhash::FxBuildHasher> {
    xplanes: Vec<A>,
    zplanes: Vec<A>,
    lplanes: Vec<A>,
    nqubits: usize,
    _hasher: PhantomData<fn() -> H>,
}

impl<A: PauliStorage, H> LossyPauliKeyColumn<A, H> {
    #[inline(always)]
    fn plane_bit(plane: &A, qubit: usize) -> bool {
        #[cfg(target_endian = "little")]
        {
            let bytes = bytemuck::bytes_of(plane);
            let byte_index = qubit / 8;
            let bit_index = qubit % 8;
            let bit_mask = 1u8 << bit_index;

            (bytes[byte_index] & bit_mask) != 0
        }
        #[cfg(target_endian = "big")]
        {
            BitArray::<A>::new(*plane)[qubit]
        }
    }

    /// Flip one bit in a packed plane; callers handle bounds and loss.
    #[inline]
    fn toggle_plane_bit(plane: &mut A, qubit: usize) {
        #[cfg(target_endian = "little")]
        {
            let byte_index = qubit / 8;
            let bit_index = qubit % 8;
            let bit_mask = 1u8 << bit_index;
            bytemuck::bytes_of_mut(plane)[byte_index] ^= bit_mask;
        }
        #[cfg(target_endian = "big")]
        {
            let mut bits = BitArray::<A>::new(*plane);
            let toggled = !bits[qubit];
            bits.set(qubit, toggled);
            *plane = bits.data;
        }
    }

    #[inline]
    fn reconstruct(&self, i: usize) -> LossyPauliWord<A, H> {
        LossyPauliWord::from_planes(
            BitArray::new(self.xplanes[i]),
            BitArray::new(self.zplanes[i]),
            BitArray::new(self.lplanes[i]),
            self.nqubits,
        )
    }
}

impl<A: PauliStorage, H> Default for LossyPauliKeyColumn<A, H> {
    #[inline]
    fn default() -> Self {
        Self {
            xplanes: Vec::new(),
            zplanes: Vec::new(),
            lplanes: Vec::new(),
            nqubits: 0,
            _hasher: PhantomData,
        }
    }
}

impl<A: PauliStorage, H> Clone for LossyPauliKeyColumn<A, H> {
    #[inline]
    fn clone(&self) -> Self {
        Self {
            xplanes: self.xplanes.clone(),
            zplanes: self.zplanes.clone(),
            lplanes: self.lplanes.clone(),
            nqubits: self.nqubits,
            _hasher: PhantomData,
        }
    }
}

impl<A, H> KeyColumn for LossyPauliKeyColumn<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    type Key = LossyPauliWord<A, H>;

    #[inline]
    fn len(&self) -> usize {
        self.xplanes.len()
    }

    #[inline]
    fn capacity(&self) -> usize {
        self.xplanes
            .capacity()
            .min(self.zplanes.capacity())
            .min(self.lplanes.capacity())
    }

    #[inline]
    fn hash_into(&self, out: &mut [u64]) {
        debug_assert_eq!(out.len(), self.len(), "hash column length mismatch");
        for (i, slot) in out.iter_mut().enumerate() {
            *slot =
                structural_hash_lossy::<A, H>(&self.xplanes[i], &self.zplanes[i], &self.lplanes[i]);
        }
    }

    #[inline]
    fn key_eq(&self, i: usize, other: &Self::Key) -> bool {
        self.nqubits == other.nqubits
            && self.xplanes[i] == other.xbits.data
            && self.zplanes[i] == other.zbits.data
            && self.lplanes[i] == other.lbits.data
    }

    #[inline]
    fn gather(&self, indices: &[u32]) -> Self {
        let mut xplanes = Vec::with_capacity(indices.len());
        let mut zplanes = Vec::with_capacity(indices.len());
        let mut lplanes = Vec::with_capacity(indices.len());
        for &idx in indices {
            xplanes.push(self.xplanes[idx as usize]);
            zplanes.push(self.zplanes[idx as usize]);
            lplanes.push(self.lplanes[idx as usize]);
        }
        Self {
            xplanes,
            zplanes,
            lplanes,
            nqubits: self.nqubits,
            _hasher: PhantomData,
        }
    }

    #[inline]
    fn get(&self, i: usize) -> Self::Key {
        self.reconstruct(i)
    }
}

impl<A, H> KeyColumnMut for LossyPauliKeyColumn<A, H>
where
    A: PauliStorage,
    H: HashFinalize + Default + BuildHasher,
{
    fn with_capacity(n: usize) -> Self {
        Self {
            xplanes: Vec::with_capacity(n),
            zplanes: Vec::with_capacity(n),
            lplanes: Vec::with_capacity(n),
            nqubits: 0,
            _hasher: PhantomData,
        }
    }

    fn push(&mut self, key: Self::Key) {
        if self.xplanes.is_empty() {
            self.nqubits = key.nqubits;
        } else {
            debug_assert_eq!(self.nqubits, key.nqubits, "column width mismatch");
        }
        self.xplanes.push(key.xbits.data);
        self.zplanes.push(key.zbits.data);
        self.lplanes.push(key.lbits.data);
    }

    /// Reserve space for additional words in all three planes.
    #[inline]
    fn reserve(&mut self, additional: usize) {
        self.xplanes.reserve(additional);
        self.zplanes.reserve(additional);
        self.lplanes.reserve(additional);
    }

    /// Empty all three plane blocks
    #[inline]
    fn clear(&mut self) {
        self.xplanes.clear();
        self.zplanes.clear();
        self.lplanes.clear();
    }

    /// Overwrite slot `i`'s planes in place: three plane stores, no move of any
    /// other element and no reallocation.
    #[inline]
    fn set(&mut self, i: usize, key: Self::Key) {
        debug_assert_eq!(self.nqubits, key.nqubits, "column width mismatch");
        self.xplanes[i] = key.xbits.data;
        self.zplanes[i] = key.zbits.data;
        self.lplanes[i] = key.lbits.data;
    }

    #[inline]
    fn truncate(&mut self, len: usize) {
        self.xplanes.truncate(len);
        self.zplanes.truncate(len);
        self.lplanes.truncate(len);
    }

    #[inline]
    fn swap_remove(&mut self, i: usize) -> Self::Key {
        let x = self.xplanes.swap_remove(i);
        let z = self.zplanes.swap_remove(i);
        let l = self.lplanes.swap_remove(i);
        LossyPauliWord::from_planes(
            BitArray::new(x),
            BitArray::new(z),
            BitArray::new(l),
            self.nqubits,
        )
    }
}

impl<A, H> PauliColumn for LossyPauliKeyColumn<A, H>
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

    #[inline]
    fn toggled_bits(&self, row: usize, qubit: usize, toggle_x: bool, toggle_z: bool) -> Self::Key {
        if qubit >= self.nqubits {
            panic!("qubit {qubit} out of bounds");
        }
        let lost = Self::plane_bit(&self.lplanes[row], qubit);
        let toggle_x = toggle_x && !lost;
        let toggle_z = toggle_z && !lost;
        let mut x = self.xplanes[row];
        let mut z = self.zplanes[row];

        if toggle_x {
            Self::toggle_plane_bit(&mut x, qubit);
        }
        if toggle_z {
            Self::toggle_plane_bit(&mut z, qubit);
        }

        LossyPauliWord::from_planes(
            BitArray::new(x),
            BitArray::new(z),
            BitArray::new(self.lplanes[row]),
            self.nqubits,
        )
    }

    #[inline]
    fn toggled_bits2(
        &self,
        row: usize,
        i: usize,
        toggle_i: [bool; 2],
        j: usize,
        toggle_j: [bool; 2],
    ) -> Self::Key {
        if i >= self.nqubits || j >= self.nqubits {
            panic!("qubit index out of bounds");
        }
        let lost_i = Self::plane_bit(&self.lplanes[row], i);
        let lost_j = Self::plane_bit(&self.lplanes[row], j);
        let mut x = self.xplanes[row];
        let mut z = self.zplanes[row];
        let toggle_x_i = toggle_i[0] && !lost_i;
        let toggle_z_i = toggle_i[1] && !lost_i;
        let toggle_x_j = toggle_j[0] && !lost_j;
        let toggle_z_j = toggle_j[1] && !lost_j;

        if toggle_x_i {
            Self::toggle_plane_bit(&mut x, i);
        }
        if toggle_z_i {
            Self::toggle_plane_bit(&mut z, i);
        }
        if toggle_x_j {
            Self::toggle_plane_bit(&mut x, j);
        }
        if toggle_z_j {
            Self::toggle_plane_bit(&mut z, j);
        }

        LossyPauliWord::from_planes(
            BitArray::new(x),
            BitArray::new(z),
            BitArray::new(self.lplanes[row]),
            self.nqubits,
        )
    }
}

impl<A, H> LossColumn for LossyPauliKeyColumn<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn is_lost(&self, row: usize, qubit: usize) -> bool {
        Self::plane_bit(&self.lplanes[row], qubit)
    }
}

impl<A, H> Columnar for LossyPauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    type Column = LossyPauliKeyColumn<A, H>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_traits_2::{Indexable, PauliBits};

    #[test]
    fn column_toggles_match_scalar_toggles_and_preserve_loss() {
        let masks = [[false, false], [true, false], [false, true], [true, true]];
        for input in ["LL", "LI", "IL", "XY"] {
            let word = LossyPauliWord::<u8>::from(input);
            let mut column = LossyPauliKeyColumn::<u8>::default();
            column.push(word.clone());

            for i in 0..2 {
                for mask_i in masks {
                    assert_eq!(
                        column.toggled_bits(0, i, mask_i[0], mask_i[1]),
                        word.toggled_bits(i, mask_i[0], mask_i[1])
                    );
                    for j in 0..2 {
                        for mask_j in masks {
                            let branch = column.toggled_bits2(0, i, mask_i, j, mask_j);
                            assert_eq!(branch, word.toggled_bits2(i, mask_i, j, mask_j));
                            for q in 0..2 {
                                assert_eq!(branch.is_lost(q), word.is_lost(q));
                                if word.is_lost(q) {
                                    assert!(!branch.x_bit(q) && !branch.z_bit(q));
                                }
                            }
                        }
                    }
                }
            }
            assert_eq!(column.get(0), word);
        }
    }

    #[test]
    fn column_toggles_reject_indices_outside_logical_width() {
        use std::panic::{AssertUnwindSafe, catch_unwind};

        let mut column = LossyPauliKeyColumn::<u8>::default();
        column.push(LossyPauliWord::from("I"));
        for invalid in [1, 7, 8] {
            assert!(
                catch_unwind(AssertUnwindSafe(|| {
                    column.toggled_bits(0, invalid, false, false)
                }))
                .is_err()
            );
            for (i, j) in [(0, invalid), (invalid, 0)] {
                assert!(
                    catch_unwind(AssertUnwindSafe(|| {
                        column.toggled_bits2(0, i, [true, true], j, [false, false])
                    }))
                    .is_err()
                );
            }
        }
        assert_eq!(column.get(0).to_string(), "I");
    }

    #[test]
    fn push_preserves_words_and_site_bits() {
        let mut column = LossyPauliKeyColumn::<u16>::default();
        let word = LossyPauliWord::from("IIIIIIIIXYL");
        column.push(word.clone());

        assert_eq!(column.len(), 1);
        assert_eq!(column.get(0), word);
        assert!(column.x_bit(0, 8));
        assert!(!column.z_bit(0, 8));
        assert!(column.x_bit(0, 9));
        assert!(column.z_bit(0, 9));
        assert!(column.is_lost(0, 10));
    }

    #[test]
    fn gather_selects_rows_in_requested_order() {
        let mut column = LossyPauliKeyColumn::<u8>::default();
        column.push(LossyPauliWord::from("XI"));
        column.push(LossyPauliWord::from("YL"));
        column.push(LossyPauliWord::from("IZ"));

        let selected = column.gather(&[2, 0, 2]);

        assert_eq!(selected.len(), 3);
        assert_eq!(selected.get(0).to_string(), "IZ");
        assert_eq!(selected.get(1).to_string(), "XI");
        assert_eq!(selected.get(2).to_string(), "IZ");
        assert_eq!(column.get(1).to_string(), "YL");
    }

    #[test]
    fn column_hashes_match_individual_word_hashes() {
        let mut column = LossyPauliKeyColumn::<u8>::default();
        column.push(LossyPauliWord::from("XLZ"));
        column.push(LossyPauliWord::from("YII"));
        let mut hashes = [0; 2];

        column.hash_into(&mut hashes);

        assert_eq!(hashes, [column.get(0).key_hash(), column.get(1).key_hash()]);
    }
}
