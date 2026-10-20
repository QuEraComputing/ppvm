// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::{BuildHasher, Hash, Hasher};
use std::sync::atomic::Ordering;

use ppvm_pauli_word_2::{HashFinalize, PauliStorage};
use ppvm_traits_2::Indexable;

use crate::LossyPauliWord;

#[inline]
pub(crate) fn structural_hash_lossy<A, H>(x: &A, z: &A, l: &A, _nqubits: usize) -> u64
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    let xz_hash = hash_xz::<A, H>(x, z);
    let loss_hash = H::default().hash_one(l);
    combine_hashes::<H>(xz_hash, loss_hash)
}

#[inline]
fn hash_xz<A: PauliStorage, H: BuildHasher + Default>(x: &A, z: &A) -> u64 {
    let mut hasher = H::default().build_hasher();
    x.hash(&mut hasher);
    z.hash(&mut hasher);
    hasher.finish()
}

/// Combine component hashes identically for individual words and column rows.
#[inline]
fn combine_hashes<H: BuildHasher + Default + HashFinalize>(xz_hash: u64, loss_hash: u64) -> u64 {
    let mut hasher = H::default().build_hasher();
    hasher.write_u64(xz_hash);
    hasher.write_u64(loss_hash);
    let word_hash = H::finalize_hash(hasher.finish(), std::mem::size_of::<[u64; 2]>());
    H::index_hash(word_hash)
}

/// Write the finalized word digest as one `u64` for `IdentityHasher`
/// to pass through unchanged.
impl<A, H> Hash for LossyPauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn hash<S: Hasher>(&self, state: &mut S) {
        state.write_u64(self.key_hash());
    }
}

/// Cache the X/Z and loss digests independently, computing each on demand.
/// Combine both digests and finalize the word hash on every call.
/// Zero marks an empty cache; a zero digest is recomputed on the next call.
impl<A, H> Indexable for LossyPauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn key_hash(&self) -> u64 {
        let xz_hash = match self.xz_hash_cache.load(Ordering::Relaxed) {
            0 => {
                let digest = hash_xz::<A, H>(&self.xbits.data, &self.zbits.data);
                self.xz_hash_cache.store(digest, Ordering::Relaxed);
                digest
            }
            cached => cached,
        };

        let loss_hash = match self.loss_hash_cache.load(Ordering::Relaxed) {
            0 => {
                let digest = H::default().hash_one(self.lbits.data);
                self.loss_hash_cache.store(digest, Ordering::Relaxed);
                digest
            }
            cached => cached,
        };

        combine_hashes::<H>(xz_hash, loss_hash)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_traits_2::{IdentityBuildHasher, LossySite, Pauli};
    use std::collections::HashMap;

    #[test]
    fn equal_words_find_the_same_map_entry() {
        let mut map = HashMap::with_hasher(IdentityBuildHasher);
        map.insert(LossyPauliWord::<u8>::from("XLZ"), 7);

        assert_eq!(map.get(&LossyPauliWord::from("X_L_Z")), Some(&7));
        assert_eq!(map.get(&LossyPauliWord::from("XIZ")), None);
    }

    #[test]
    fn hashing_after_mutation_matches_a_fresh_word() {
        let mut word = LossyPauliWord::<u8>::from("XLZ");
        let _ = word.key_hash();

        word.clear_loss(1);
        assert_eq!(
            word.key_hash(),
            LossyPauliWord::<u8>::from("XIZ").key_hash()
        );

        word.set(0, LossySite::Present(Pauli::Y));
        assert_eq!(
            word.key_hash(),
            LossyPauliWord::<u8>::from("YIZ").key_hash()
        );
    }
}
