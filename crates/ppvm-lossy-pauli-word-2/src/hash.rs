// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The lazy structural hash: one ordered fold over `(X, Z, loss)`, the
//! [`Hash`](std::hash::Hash) impl that writes exactly `key_hash()`, and the
//! [`Indexable`] impl that caches the finalized digest.
//!
//! Design: `word-data-structures.md` §"Component hashes" and §"Lazy hashing".
//! Fixed-size planes are written in a fixed order, so the fold is
//! domain-separated by position rather than by allocating three intermediate
//! hashers. Only the consumed final value is cached, because retaining
//! intermediate atomics regressed clone-and-mutate key paths;
//! `traits-2-configuration-and-hashing.md` §"Indexable values" (the `key_hash()`
//! value contract) and §"Concrete word hashing" (the private, per-hasher/per-width
//! finalization fold, reused here from `ppvm-pauli-word-2`'s [`HashFinalize`]).

use std::hash::{BuildHasher, Hash, Hasher};
use std::sync::atomic::Ordering;

use ppvm_pauli_word_2::{HashFinalize, PauliStorage};
use ppvm_traits_2::Indexable;

use crate::data::LossyPauliWord;

/// The finalized structural digest of `(nqubits, x, z, loss)`.
///
/// Factored out so [`Indexable::key_hash`] and `KeyColumn::hash_into` compute the
/// **bit-for-bit identical** digest from a scalar word and from a column plane
/// (the agreement the batch contract requires; `word-data-structures.md`
/// §"Key columns").
#[inline]
pub(crate) fn structural_hash_lossy<A, H>(x: &A, z: &A, l: &A, _nqubits: usize) -> u64
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    let mut hasher = H::default().build_hasher();
    x.hash(&mut hasher);
    z.hash(&mut hasher);
    l.hash(&mut hasher);
    H::finalize_hash(hasher.finish(), std::mem::size_of::<A>())
}

/// `Hash` writes exactly the finalized `key_hash()` as a single `u64`, so the
/// digest reaches hashbrown untouched through `IdentityHasher`
/// (`traits-2-configuration-and-hashing.md` §"Indexable values").
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

/// The finalized structural digest, computed lazily and cached in one
/// `AtomicU64`, so a warm read (the map-lookup hot path) is one field load. The
/// one ordered hasher fold avoids both intermediate atomic cells and three
/// separate hasher initializations. Structurally equal keys return equal digests.
impl<A, H> Indexable for LossyPauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn key_hash(&self) -> u64 {
        let cached = self.hash_cache.load(Ordering::Relaxed);
        if cached != 0 {
            return cached;
        }
        let digest = structural_hash_lossy::<A, H>(
            &self.xbits.data,
            &self.zbits.data,
            &self.lbits.data,
            self.nqubits,
        );
        self.hash_cache.store(digest, Ordering::Relaxed);
        digest
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_traits_2::{IdentityBuildHasher, LossySite, Pauli};
    use std::collections::HashMap;

    #[test]
    fn equal_words_equal_digest() {
        let a: LossyPauliWord = "XLZL".into();
        let b: LossyPauliWord = "XLZL".into();
        assert_eq!(a.key_hash(), b.key_hash());
    }

    #[test]
    fn loss_distinguishes_digest() {
        let present: LossyPauliWord = "XIZI".into();
        let lossy: LossyPauliWord = "XLZL".into();
        assert_ne!(present.key_hash(), lossy.key_hash());
    }

    #[test]
    fn hash_writes_key_hash() {
        let w: LossyPauliWord = "XLZL".into();
        let bh = IdentityBuildHasher;
        assert_eq!(bh.hash_one(&w), w.key_hash());
    }

    #[test]
    fn loss_only_mutation_invalidates_digest() {
        let mut w: LossyPauliWord = "XIZI".into();
        let before = w.key_hash();
        assert_ne!(w.hash_cache.load(Ordering::Relaxed), 0, "digest cached");
        w.set_lost(1); // identity site -> lost (X/Z unchanged)
        assert_eq!(w.hash_cache.load(Ordering::Relaxed), 0);
        assert_ne!(w.key_hash(), before);
    }

    #[test]
    fn cache_warm_read_and_invalidation() {
        let mut w: LossyPauliWord = "XIZI".into();
        assert_eq!(w.hash_cache.load(Ordering::Relaxed), 0, "starts cold");
        let h0 = w.key_hash();
        assert_eq!(
            w.hash_cache.load(Ordering::Relaxed),
            h0,
            "warm read caches the finalized combined digest",
        );
        w.set_lost(1);
        assert_eq!(w.hash_cache.load(Ordering::Relaxed), 0);
        assert_ne!(w.key_hash(), h0, "recompute reflects the loss mutation");
    }

    #[test]
    fn nonidentity_loss_invalidates_digest() {
        let mut w: LossyPauliWord = "XIZI".into();
        let _ = w.key_hash();
        w.set_lost(0); // X -> lost, nonidentity: both components change
        assert_eq!(w.hash_cache.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn cache_stable_across_clone_and_mutation() {
        let w: LossyPauliWord = "XLZL".into();
        let h0 = w.key_hash();
        let c = w.clone();
        assert_eq!(c.key_hash(), h0, "clone copies the cached digest");

        let mut m = w.clone();
        m.set(0, LossySite::Present(Pauli::Y)); // X -> Y, structural change
        assert_ne!(m.key_hash(), h0, "mutation invalidates and recomputes");
    }

    #[test]
    #[allow(clippy::mutable_key_type)]
    fn usable_as_identity_hashmap_key() {
        let mut map: HashMap<LossyPauliWord, i32, IdentityBuildHasher> = HashMap::default();
        map.insert("XLZL".into(), 7);
        assert_eq!(map.get(&LossyPauliWord::from("XLZL")), Some(&7));
    }
}
