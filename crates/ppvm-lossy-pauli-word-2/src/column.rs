// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The structure-of-arrays key column for [`LossyPauliWord`]:
//! [`LossyPauliKeyColumn`] stores parallel X, Z, **and loss** plane blocks and a
//! shared width, and implements [`KeyColumn`]; [`LossyPauliWord`] implements
//! [`Columnar`] pointing at it.
//!
//! Design: `word-data-structures.md` §"Key columns (structure-of-arrays
//! batches)" — "The flattened `LossyPauliWord` extends its column with a parallel
//! loss-bit plane, matching its X/Z/loss structural identity." Mirrors
//! `ppvm-pauli-word-2`'s `PauliKeyColumn` with the extra loss plane; the first
//! cut keeps each key's blob as one contiguous `A` slot per plane (the
//! naive-but-correct fallback the batch contract permits), leaving SIMD/alignment
//! plane packing to the `ColumnStore` backend (implementation-plan Phase 6).
//!
//! `hash_into` reuses [`crate::hash::structural_hash_lossy`], the exact fold
//! `Indexable::key_hash` uses, so the plane-parallel hash and the scalar `Hash`
//! agree bit for bit (`word-data-structures.md` §"Key columns").

use std::hash::BuildHasher;
use std::marker::PhantomData;

use bitvec::array::BitArray;
use ppvm_pauli_word_2::{HashFinalize, PauliStorage};
use ppvm_traits_2::{Columnar, KeyColumn};

use crate::data::LossyPauliWord;
use crate::hash::structural_hash_lossy;

/// A structure-of-arrays column of [`LossyPauliWord`]s: parallel X, Z, and loss
/// plane blocks plus the shared qubit width.
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
            bytes[qubit >> 3] & (1 << (qubit & 7)) != 0
        }
        #[cfg(target_endian = "big")]
        {
            BitArray::<A>::new(*plane)[qubit]
        }
    }

    #[inline]
    fn reconstruct(&self, i: usize) -> LossyPauliWord<A, H> {
        let mut x = BitArray::<A>::ZERO;
        let mut z = BitArray::<A>::ZERO;
        let mut l = BitArray::<A>::ZERO;
        x.data = self.xplanes[i];
        z.data = self.zplanes[i];
        l.data = self.lplanes[i];
        LossyPauliWord::from_planes(x, z, l, self.nqubits)
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
    fn with_capacity(n: usize) -> Self {
        Self {
            xplanes: Vec::with_capacity(n),
            zplanes: Vec::with_capacity(n),
            lplanes: Vec::with_capacity(n),
            nqubits: 0,
            _hasher: PhantomData,
        }
    }

    /// Append one key's planes. The first key fixes the column width; later keys
    /// must match it (they always do inside one sum).
    #[inline]
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

    /// Pre-size all three plane blocks — one reallocation instead of a doubling
    /// chain when the caller knows how many keys are about to be appended.
    #[inline]
    fn reserve(&mut self, additional: usize) {
        self.xplanes.reserve(additional);
        self.zplanes.reserve(additional);
        self.lplanes.reserve(additional);
    }

    /// Fill `out[i]` with the `i`-th key's `key_hash()`, computed from the planes
    /// via the same fold — so `out[i] == self.get(i).key_hash()` bit for bit.
    #[inline]
    fn hash_into(&self, out: &mut [u64]) {
        debug_assert_eq!(out.len(), self.len(), "hash column length mismatch");
        for (i, slot) in out.iter_mut().enumerate() {
            *slot = structural_hash_lossy::<A, H>(
                &self.xplanes[i],
                &self.zplanes[i],
                &self.lplanes[i],
                self.nqubits,
            );
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

    #[inline(always)]
    fn x_bit(&self, row: usize, qubit: usize) -> bool {
        Self::plane_bit(&self.xplanes[row], qubit)
    }

    #[inline(always)]
    fn z_bit(&self, row: usize, qubit: usize) -> bool {
        Self::plane_bit(&self.zplanes[row], qubit)
    }

    #[inline(always)]
    fn is_lost(&self, row: usize, qubit: usize) -> bool {
        Self::plane_bit(&self.lplanes[row], qubit)
    }

    #[inline]
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
        LossyPauliWord::from_planes(
            BitArray::new(x),
            BitArray::new(z),
            BitArray::new(self.lplanes[row]),
            self.nqubits,
        )
    }

    /// Empty all three plane blocks, keeping their allocations — the in-place
    /// surface a column that is a live *store* needs (`ColumnStore`, Phase 6).
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
    use ppvm_traits_2::Indexable;

    fn column(words: &[&str]) -> LossyPauliKeyColumn<u64> {
        let mut col = LossyPauliKeyColumn::<u64>::with_capacity(words.len());
        for w in words {
            col.push(LossyPauliWord::from(*w));
        }
        col
    }

    #[test]
    fn roundtrip_and_len() {
        let words = ["XLZI", "IIIL", "YYLY"];
        let col = column(&words);
        assert_eq!(col.len(), 3);
        assert!(!col.is_empty());
        for (i, w) in words.iter().enumerate() {
            assert_eq!(col.get(i), LossyPauliWord::from(*w));
        }
    }

    #[test]
    fn hash_into_matches_scalar_key_hash() {
        let words = ["XLZI", "IIIL", "YYLY", "ZXLX"];
        let col = column(&words);
        let mut out = vec![0u64; col.len()];
        col.hash_into(&mut out);
        for (i, w) in words.iter().enumerate() {
            assert_eq!(
                out[i],
                LossyPauliWord::<u64>::from(*w).key_hash(),
                "key {i}"
            );
        }
    }

    #[test]
    fn key_eq_and_gather() {
        let col = column(&["XLZI", "IIIL", "YYLY"]);
        assert!(col.key_eq(1, &LossyPauliWord::from("IIIL")));
        assert!(!col.key_eq(1, &LossyPauliWord::from("XLZI")));

        let picked = col.gather(&[2, 0]);
        assert_eq!(picked.len(), 2);
        assert_eq!(picked.get(0), LossyPauliWord::from("YYLY"));
        assert_eq!(picked.get(1), LossyPauliWord::from("XLZI"));
    }
}
