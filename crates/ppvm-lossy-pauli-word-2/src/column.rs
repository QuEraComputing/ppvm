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

impl<A, H> KeyColumn for LossyPauliKeyColumn<A, H> where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize {
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

        
}

impl<A, H> KeyColumnMut for LossyPauliKeyColumn<A, H> where 
    A: PauliStorage,
    H: HashFinalize + Default + BuildHasher {
        fn with_capacity(n: usize) -> Self {
            Self {
                    xplanes: Vec::with_capacity(n),
                    zplanes: Vec::with_capacity(n),
                    lplanes: Vec::with_capacity(n),
                    nqubits: 0,
                    _hasher: PhantomData
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

        // Pre-size all three plane blocks — one reallocation instead of a doubling
        /// chain when the caller knows how many keys are about to be 
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

impl<A, H> PauliColumn for LossyPauliKeyColumn<A, H> where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize {
        
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
            let mut x = self.xplanes[row];
            let mut z = self.zplanes[row];

            #[cfg(target_endian = "little")]
            {
                let byte_index = qubit / 8;
                let bit_index = qubit % 8;
                let bit_mask = 1u8 << bit_index;

                if toggle_x {
                    bytemuck::bytes_of_mut(&mut x)[byte_index] ^= bit_mask;
                }
                if toggle_z {
                    bytemuck::bytes_of_mut(&mut z)[byte_index] ^= bit_mask;
                }
            }

            #[cfg(target_endian = "big")]
            {
                let mut x_bits = BitArray::<A>::new(x);
                let mut z_bits = BitArray::<A>::new(z);

                if toggle_x {
                    let toggled = !x_bits[qubit];
                    x_bits.set(qubit, toggled);
                }
                if toggle_z {
                    let toggled = !z_bits[qubit];
                    z_bits.set(qubit, toggled);
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
            toggle_i: [bool; 2],
            j: usize,
            toggle_j: [bool; 2],
        ) -> Self::Key {
            let mut x = self.xplanes[row];
            let mut z = self.zplanes[row];
            let [toggle_x_i, toggle_z_i] = toggle_i;
            let [toggle_x_j, toggle_z_j] = toggle_j;

            #[cfg(target_endian = "little")]
            {
                let byte_index_i = i / 8;
                let bit_index_i = i % 8;
                let bit_mask_i = 1u8 << bit_index_i;
                let byte_index_j = j / 8;
                let bit_index_j = j % 8;
                let bit_mask_j = 1u8 << bit_index_j;

                if toggle_x_i {
                    bytemuck::bytes_of_mut(&mut x)[byte_index_i] ^= bit_mask_i;
                }
                if toggle_z_i {
                    bytemuck::bytes_of_mut(&mut z)[byte_index_i] ^= bit_mask_i;
                }
                if toggle_x_j {
                    bytemuck::bytes_of_mut(&mut x)[byte_index_j] ^= bit_mask_j;
                }
                if toggle_z_j {
                    bytemuck::bytes_of_mut(&mut z)[byte_index_j] ^= bit_mask_j;
                }
            }

            #[cfg(target_endian = "big")]
            {
                let mut x_bits = BitArray::<A>::new(x);
                let mut z_bits = BitArray::<A>::new(z);

                if toggle_x_i {
                    let toggled = !x_bits[i];
                    x_bits.set(i, toggled);
                }
                if toggle_z_i {
                    let toggled = !z_bits[i];
                    z_bits.set(i, toggled);
                }
                if toggle_x_j {
                    let toggled = !x_bits[j];
                    x_bits.set(j, toggled);
                }
                if toggle_z_j {
                    let toggled = !z_bits[j];
                    z_bits.set(j, toggled);
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
}

impl<A, H> LossColumn for LossyPauliKeyColumn<A, H> where 
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize {
    
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
