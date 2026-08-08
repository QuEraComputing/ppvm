// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The structure-of-arrays key column for [`PauliWord`]: [`PauliKeyColumn`]
//! stores the X and Z plane blocks in parallel and a shared width, and
//! implements [`KeyColumn`]; [`PauliWord`] implements [`Columnar`] pointing at
//! it.
//!
//! Design: `word-data-structures.md` §"Key columns (structure-of-arrays
//! batches)" — the column is "two plane blocks and a shared width". This first
//! cut keeps each key's X/Z blob as one contiguous `A` slot per plane (the
//! naive-but-correct fallback the batch contract explicitly permits); the
//! SIMD/alignment plane packing and the parallel hash-column ownership are the
//! `ColumnStore` backend's job (implementation-plan Phase 6, and
//! `word-data-structures.md` open questions 2/3).
//!
//! `hash_into` reuses [`crate::hash::structural_hash`], the exact fold
//! `Indexable::key_hash` uses, so the plane-parallel hash and the scalar `Hash`
//! agree bit for bit (Design: `word-data-structures.md` §"Key columns").

use std::marker::PhantomData;

use bitvec::array::BitArray;
use ppvm_traits_2::{Columnar, KeyColumn};

use crate::data::PauliWord;
use crate::hash::structural_hash;
use crate::storage::{HashFinalize, PauliStorage};
use std::hash::BuildHasher;

/// A structure-of-arrays column of [`PauliWord`]s: parallel X and Z plane blocks
/// plus the shared qubit width.
pub struct PauliKeyColumn<A: PauliStorage, H = fxhash::FxBuildHasher> {
    xplanes: Vec<A>,
    zplanes: Vec<A>,
    nqubits: usize,
    _hasher: PhantomData<fn() -> H>,
}

impl<A: PauliStorage, H> PauliKeyColumn<A, H> {
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

    #[inline]
    fn reconstruct(&self, i: usize) -> PauliWord<A, H> {
        let mut x = BitArray::<A>::ZERO;
        let mut z = BitArray::<A>::ZERO;
        x.data = self.xplanes[i];
        z.data = self.zplanes[i];
        PauliWord::from_planes(x, z, self.nqubits)
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

    #[inline]
    fn with_capacity(n: usize) -> Self {
        Self {
            xplanes: Vec::with_capacity(n),
            zplanes: Vec::with_capacity(n),
            nqubits: 0,
            _hasher: PhantomData,
        }
    }

    /// Append one key's planes. The first key fixes the column width; later keys
    /// must match it (they always do inside one sum — see `Sum`'s width
    /// invariant in the design).
    #[inline]
    fn push(&mut self, key: Self::Key) {
        if self.xplanes.is_empty() {
            self.nqubits = key.nqubits;
        } else {
            debug_assert_eq!(self.nqubits, key.nqubits, "column width mismatch");
        }
        self.xplanes.push(key.xbits.data);
        self.zplanes.push(key.zbits.data);
    }

    /// Pre-size both plane blocks — one reallocation instead of a doubling chain
    /// when the caller knows how many keys are about to be appended.
    #[inline]
    fn reserve(&mut self, additional: usize) {
        self.xplanes.reserve(additional);
        self.zplanes.reserve(additional);
    }

    /// Fill `out[i]` with the `i`-th key's `key_hash()`, computed from the planes
    /// via the same fold — so `out[i] == self.get(i).key_hash()` bit for bit.
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
        self.reconstruct(i)
    }

    #[inline(always)]
    fn x_bit(&self, row: usize, qubit: usize) -> bool {
        Self::plane_bit(&self.xplanes[row], qubit)
    }

    #[inline(always)]
    fn z_bit(&self, row: usize, qubit: usize) -> bool {
        Self::plane_bit(&self.zplanes[row], qubit)
    }

    #[inline(always)]
    fn is_lost(&self, _row: usize, _qubit: usize) -> bool {
        false
    }

    #[inline(always)]
    fn toggled_bits(&self, row: usize, qubit: usize, toggle_x: bool, toggle_z: bool) -> Self::Key {
        let mut x = self.xplanes[row];
        let mut z = self.zplanes[row];
        #[cfg(target_endian = "little")]
        {
            let mask = 1 << (qubit & 7);
            if toggle_x {
                bytemuck::bytes_of_mut(&mut x)[qubit >> 3] ^= mask;
            }
            if toggle_z {
                bytemuck::bytes_of_mut(&mut z)[qubit >> 3] ^= mask;
            }
        }
        #[cfg(target_endian = "big")]
        {
            let mut x_bits = BitArray::<A>::new(x);
            let mut z_bits = BitArray::<A>::new(z);
            if toggle_x {
                x_bits.set(qubit, !x_bits[qubit]);
            }
            if toggle_z {
                z_bits.set(qubit, !z_bits[qubit]);
            }
            x = x_bits.data;
            z = z_bits.data;
        }
        PauliWord::from_planes(BitArray::new(x), BitArray::new(z), self.nqubits)
    }

    #[inline(always)]
    fn toggled_bits2(
        &self,
        row: usize,
        i: usize,
        toggle_x_i: bool,
        toggle_z_i: bool,
        j: usize,
        toggle_x_j: bool,
        toggle_z_j: bool,
    ) -> Self::Key {
        let mut x = self.xplanes[row];
        let mut z = self.zplanes[row];
        #[cfg(target_endian = "little")]
        {
            if toggle_x_i {
                bytemuck::bytes_of_mut(&mut x)[i >> 3] ^= 1 << (i & 7);
            }
            if toggle_z_i {
                bytemuck::bytes_of_mut(&mut z)[i >> 3] ^= 1 << (i & 7);
            }
            if toggle_x_j {
                bytemuck::bytes_of_mut(&mut x)[j >> 3] ^= 1 << (j & 7);
            }
            if toggle_z_j {
                bytemuck::bytes_of_mut(&mut z)[j >> 3] ^= 1 << (j & 7);
            }
        }
        #[cfg(target_endian = "big")]
        {
            let mut x_bits = BitArray::<A>::new(x);
            let mut z_bits = BitArray::<A>::new(z);
            if toggle_x_i {
                x_bits.set(i, !x_bits[i]);
            }
            if toggle_z_i {
                z_bits.set(i, !z_bits[i]);
            }
            if toggle_x_j {
                x_bits.set(j, !x_bits[j]);
            }
            if toggle_z_j {
                z_bits.set(j, !z_bits[j]);
            }
            x = x_bits.data;
            z = z_bits.data;
        }
        PauliWord::from_planes(BitArray::new(x), BitArray::new(z), self.nqubits)
    }

    /// Empty both plane blocks, keeping their allocations (and the column's
    /// width, so a cleared-and-refilled column of the same sum keeps its
    /// `nqubits` even while empty — the `ColumnStore`'s aux buffer relies on the
    /// allocation surviving, not on the width).
    #[inline]
    fn clear(&mut self) {
        self.xplanes.clear();
        self.zplanes.clear();
    }

    /// Overwrite slot `i`'s planes in place — two plane stores, no move of any
    /// other element and no reallocation. The write side of the `ColumnStore`'s
    /// in-place Clifford re-key.
    #[inline]
    fn set(&mut self, i: usize, key: Self::Key) {
        debug_assert_eq!(self.nqubits, key.nqubits, "column width mismatch");
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
    use ppvm_traits_2::Indexable;

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
