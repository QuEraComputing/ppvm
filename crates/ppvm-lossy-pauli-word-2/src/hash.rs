// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The lazy, component-split structural hash: the raw `(nqubits, X, Z)` and loss
//! components, their ordered domain-separated [`combine_components`] fold, the
//! [`Hash`](std::hash::Hash) impl that writes exactly `key_hash()`, and the
//! [`Indexable`] impl that finalizes and caches per component.
//!
//! Design: `word-data-structures.md` §"Component hashes"
//! (`lossy hash = combine(Pauli hash, loss hash)`, "`combine` must be ordered and
//! domain-separated … not an unqualified XOR"; the loss component is cached
//! separately so a loss-only mutation avoids rehashing X/Z) and §"Lazy hashing";
//! `traits-2-configuration-and-hashing.md` §"Indexable values" (the `key_hash()`
//! value contract) and §"Concrete word hashing" (the private, per-hasher/per-width
//! finalization fold, reused here from `ppvm-pauli-word-2`'s [`HashFinalize`]).

use std::hash::{BuildHasher, Hash, Hasher};
use std::sync::atomic::Ordering;

use ppvm_pauli_word_2::{HashFinalize, PauliStorage};
use ppvm_traits_2::Indexable;

use crate::data::LossyPauliWord;

/// Raw digest of the `(nqubits, X, Z)` component (no finalization — that happens
/// once in [`combine_components`]). Mirrors `ppvm-pauli-word-2`'s ordinary
/// structural hash so an all-present lossy word and the corresponding
/// `PauliWord` fold the same X/Z bytes.
#[inline]
pub(crate) fn xz_component<A, H>(x: &A, z: &A, nqubits: usize) -> u64
where
    A: PauliStorage,
    H: BuildHasher + Default,
{
    let mut hasher = H::default().build_hasher();
    hasher.write_usize(nqubits);
    hasher.write(bytemuck::bytes_of(x));
    hasher.write(bytemuck::bytes_of(z));
    hasher.finish()
}

/// Raw digest of the loss component.
#[inline]
pub(crate) fn loss_component<A, H>(l: &A) -> u64
where
    A: PauliStorage,
    H: BuildHasher + Default,
{
    let mut hasher = H::default().build_hasher();
    hasher.write(bytemuck::bytes_of(l));
    hasher.finish()
}

/// Ordered, domain-separated combination of the two raw components into the
/// finalized digest. The two components are written at distinct positions
/// (`xz` first, `loss` second) through the internal hasher — an ordered,
/// domain-separated fold, **not** an unqualified XOR — then finalized once by the
/// per-hasher/per-width [`HashFinalize`] fold so the low bits (hashbrown's
/// bucket) are avalanche-quality (`word-data-structures.md` §"Component hashes";
/// `traits-2-configuration-and-hashing.md` §"Concrete word hashing").
#[inline]
pub(crate) fn combine_components<H>(xz: u64, loss: u64, storage_bytes: usize) -> u64
where
    H: BuildHasher + Default + HashFinalize,
{
    let mut hasher = H::default().build_hasher();
    hasher.write_u64(xz);
    hasher.write_u64(loss);
    H::finalize_hash(hasher.finish(), storage_bytes)
}

/// The finalized structural digest of the planes `(x, z, loss)`. Width remains
/// part of equality but, as for the ordinary word, not of the digest.
///
/// Factored out so [`Indexable::key_hash`] and `KeyColumn::hash_into` compute the
/// **bit-for-bit identical** digest from a scalar word and from a column plane
/// (the agreement the batch contract requires; `word-data-structures.md`
/// §"Key columns").
#[inline]
pub(crate) fn structural_hash_lossy<A, H>(x: &A, z: &A, l: &A, nqubits: usize) -> u64
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    combine_components::<H>(
        xz_component::<A, H>(x, z, nqubits),
        loss_component::<A, H>(l),
        std::mem::size_of::<A>(),
    )
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

/// The finalized structural digest, computed lazily and cached **per component**:
/// the `(nqubits, X, Z)` component and the loss component each populate their own
/// `AtomicU64`, and `key_hash()` combines them into a third `AtomicU64`
/// so a warm read (the map-lookup hot path) is a single field load. A loss-only
/// mutation clears only the loss cell and the combined cell, so the X/Z digest is
/// reused (`word-data-structures.md` §"Component hashes";
/// `traits-2-configuration-and-hashing.md` §"Indexable values"). Structurally
/// equal keys return equal digests.
impl<A, H> Indexable for LossyPauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn key_hash(&self) -> u64 {
        let cached = self.combined_hash_cache.load(Ordering::Relaxed);
        if cached != 0 {
            return cached;
        }
        let mut xz = self.xz_hash_cache.load(Ordering::Relaxed);
        if xz == 0 {
            xz = xz_component::<A, H>(&self.xbits.data, &self.zbits.data, self.nqubits);
            self.xz_hash_cache.store(xz, Ordering::Relaxed);
        }
        let mut loss = self.loss_hash_cache.load(Ordering::Relaxed);
        if loss == 0 {
            loss = loss_component::<A, H>(&self.lbits.data);
            self.loss_hash_cache.store(loss, Ordering::Relaxed);
        }
        let digest = combine_components::<H>(xz, loss, std::mem::size_of::<A>());
        self.combined_hash_cache.store(digest, Ordering::Relaxed);
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
    fn loss_only_mutation_preserves_xz_component() {
        // Marking an *identity* site lost must not rehash X/Z: the cached X/Z
        // component cell is preserved across the loss write.
        let mut w: LossyPauliWord = "XIZI".into();
        let _ = w.key_hash(); // populate both component caches
        let xz_before = w.xz_hash_cache.load(Ordering::Relaxed);
        assert_ne!(xz_before, 0, "xz cached");
        w.set_lost(1); // identity site -> lost (X/Z unchanged)
        assert_eq!(
            w.xz_hash_cache.load(Ordering::Relaxed),
            xz_before,
            "loss-only mutation must preserve the X/Z hash component",
        );
        assert_eq!(w.loss_hash_cache.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn combined_cache_warm_read_and_invalidation() {
        // Warm read populates the combined digest cell; a mutation clears it.
        let mut w: LossyPauliWord = "XIZI".into();
        assert_eq!(
            w.combined_hash_cache.load(Ordering::Relaxed),
            0,
            "starts cold"
        );
        let h0 = w.key_hash();
        assert_eq!(
            w.combined_hash_cache.load(Ordering::Relaxed),
            h0,
            "warm read caches the finalized combined digest",
        );
        w.set_lost(1); // loss-only mutation invalidates the combined cell
        assert_eq!(w.combined_hash_cache.load(Ordering::Relaxed), 0);
        assert_ne!(w.xz_hash_cache.load(Ordering::Relaxed), 0);
        assert_ne!(w.key_hash(), h0, "recompute reflects the loss mutation");
    }

    #[test]
    fn nonidentity_loss_invalidates_both() {
        let mut w: LossyPauliWord = "XIZI".into();
        let _ = w.key_hash();
        w.set_lost(0); // X -> lost, nonidentity: both components change
        assert_eq!(w.xz_hash_cache.load(Ordering::Relaxed), 0);
        assert_eq!(w.loss_hash_cache.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn cache_stable_across_clone_and_mutation() {
        let w: LossyPauliWord = "XLZL".into();
        let h0 = w.key_hash();
        let c = w.clone();
        assert_eq!(c.key_hash(), h0, "clone copies the cached components");

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
