// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::{BuildHasher, Hash, Hasher};
use std::sync::atomic::Ordering;

use ppvm_pauli_word_2::{HashFinalize, PauliStorage};
use ppvm_traits_2::Indexable;

use crate::LossyPauliWord;

#[inline]
pub(crate) fn structural_hash_lossy<A, H>(x: &A, z: &A, l: &A) -> u64
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    let xz_hash = hash_xz::<A, H>(x, z);
    let loss_hash = H::word_hash(&[bytemuck::bytes_of(l)], size_of::<A>());
    combine_hashes::<H>(xz_hash, loss_hash)
}

#[inline]
fn hash_xz<A: PauliStorage, H: HashFinalize + BuildHasher + Default>(x: &A, z: &A) -> u64 {
    let x_bytes = bytemuck::bytes_of(x);
    let z_bytes = bytemuck::bytes_of(z);
    H::word_hash(&[x_bytes, z_bytes], size_of::<A>())
}

/// Combine component hashes identically for individual words and column rows.
#[inline]
fn combine_hashes<H: BuildHasher + Default + HashFinalize>(xz_hash: u64, loss_hash: u64) -> u64 {
    let word_hash = H::word_hash(&[&xz_hash.to_le_bytes(), &loss_hash.to_le_bytes()], 8);
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
/// Validity flags distinguish an empty cache from a cached zero digest.
impl<A, H> Indexable for LossyPauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn key_hash(&self) -> u64 {
        let xz_hash = if self.xz_hash_valid.load(Ordering::Acquire) {
            self.xz_hash_cache.load(Ordering::Relaxed)
        } else {
            let digest = hash_xz::<A, H>(&self.xbits.data, &self.zbits.data);
            self.xz_hash_cache.store(digest, Ordering::Relaxed);
            self.xz_hash_valid.store(true, Ordering::Release);
            digest
        };

        let loss_hash = if self.loss_hash_valid.load(Ordering::Acquire) {
            self.loss_hash_cache.load(Ordering::Relaxed)
        } else {
            let digest = H::word_hash(&[bytemuck::bytes_of(&self.lbits.data)], size_of::<A>());
            self.loss_hash_cache.store(digest, Ordering::Relaxed);
            self.loss_hash_valid.store(true, Ordering::Release);
            digest
        };

        combine_hashes::<H>(xz_hash, loss_hash)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_traits_2::{IdentityBuildHasher, LossySite, Pauli, PauliBits};
    use std::collections::HashMap;

    #[test]
    fn zero_hashes_are_cached_and_reused() {
        use std::cell::Cell;

        thread_local! {
            static BUILDS: Cell<usize> = const { Cell::new(0) };
        }
        #[derive(Default)]
        struct ZeroBuildHasher;
        struct ZeroHasher;

        impl Hasher for ZeroHasher {
            fn write(&mut self, _: &[u8]) {}
            fn finish(&self) -> u64 {
                0
            }
        }
        impl BuildHasher for ZeroBuildHasher {
            type Hasher = ZeroHasher;
            fn build_hasher(&self) -> ZeroHasher {
                BUILDS.with(|count| count.set(count.get() + 1));
                ZeroHasher
            }
        }
        impl HashFinalize for ZeroBuildHasher {
            fn index_hash(raw: u64) -> u64 {
                raw
            }
        }

        let mut word = LossyPauliWord::<u8, ZeroBuildHasher>::from("II");
        BUILDS.with(|count| count.set(0));
        assert_eq!(word.key_hash(), 0);
        assert!(word.xz_hash_valid.load(Ordering::Acquire));
        assert!(word.loss_hash_valid.load(Ordering::Acquire));
        BUILDS.with(|count| assert_eq!(count.get(), 3));

        assert_eq!(word.key_hash(), 0);
        BUILDS.with(|count| assert_eq!(count.get(), 4)); // Only combine cached hashes.
        assert_eq!(word.clone().key_hash(), 0);
        BUILDS.with(|count| assert_eq!(count.get(), 5));

        word.set_lost(0);
        assert_eq!(word.key_hash(), 0);
        BUILDS.with(|count| assert_eq!(count.get(), 7)); // Recompute loss and combine.
        word.set_x_bit(1, true);
        assert_eq!(word.key_hash(), 0);
        BUILDS.with(|count| assert_eq!(count.get(), 9)); // Recompute X/Z and combine.
    }

    #[test]
    fn no_loss_fxhash_is_valid_even_when_zero() {
        let word = LossyPauliWord::<u8>::from("XI");
        word.key_hash();
        assert!(word.loss_hash_valid.load(Ordering::Acquire));
        assert_eq!(word.loss_hash_cache.load(Ordering::Relaxed), 0);
    }

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
