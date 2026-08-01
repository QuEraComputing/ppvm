// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The lazy structural hash: the shared [`structural_hash`] fold, the
//! [`Hash`](std::hash::Hash) impl that writes exactly `key_hash()`, and the
//! [`Indexable`] impl that finalizes and caches the avalanche-quality digest.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Indexable values" (the
//! `key_hash()` value contract) and §"Concrete word hashing" (the private,
//! per-algorithm/per-width finalization fold); `word-data-structures.md`
//! §"Component hashes" (`packed Pauli hash = hash(nqubits, X bits, Z bits)`) and
//! §"Lazy hashing and interior mutability" (a lazy `AtomicU64` sentinel cache
//! realizing the design's interior-mutable lazy-cache contract).

use std::hash::{BuildHasher, Hash, Hasher};
use std::sync::atomic::Ordering;

use ppvm_traits_2::Indexable;

use crate::data::{HASH_UNCACHED, PauliWord};
use crate::storage::{HashFinalize, PauliStorage};

/// The finalized structural digest of the planes `(nqubits, x, z)`.
///
/// Factored out so `Indexable::key_hash` and `KeyColumn::hash_into` compute the
/// **bit-for-bit identical** digest from, respectively, a scalar word and a
/// column plane — the agreement the batch contract requires (Design:
/// §"Concrete word hashing"; `word-data-structures.md` §"Key columns"). The raw
/// digest is folded per-hasher/per-width by [`HashFinalize`] so the low bits
/// (hashbrown's bucket) are avalanche-quality even for a short key consumed
/// directly.
#[inline]
pub(crate) fn structural_hash<A, H>(x: &A, z: &A, nqubits: usize) -> u64
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    let mut hasher = H::default().build_hasher();
    // Domain-separate the width from the planes so words of different widths
    // (structurally distinct) do not collide. The planes go through the
    // byte-slice path so the internal hasher consumes machine words, not bytes.
    hasher.write_usize(nqubits);
    hasher.write(bytemuck::bytes_of(x));
    hasher.write(bytemuck::bytes_of(z));
    H::finalize_hash(hasher.finish(), std::mem::size_of::<A>())
}

/// `Hash` writes exactly the finalized `key_hash()` as a single `u64`, so the
/// digest reaches hashbrown untouched through `IdentityHasher` (Design:
/// §"Indexable values": "`Hash for Self` is exactly
/// `state.write_u64(self.key_hash())`").
impl<A, H> Hash for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn hash<S: Hasher>(&self, state: &mut S) {
        state.write_u64(self.key_hash());
    }
}

/// The finalized structural digest, lazily computed once and cached.
///
/// Design: §"Indexable values". Structurally equal keys return equal digests
/// (equality and this hash both read `(nqubits, X, Z)`); the value is
/// avalanche-quality via the private [`HashFinalize`] fold.
impl<A, H> Indexable for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn key_hash(&self) -> u64 {
        // Relaxed is sufficient: the cache is a pure function of the immutable
        // structural fields, so any thread that observes a non-sentinel value
        // observes *the* digest. A racing miss just recomputes the same value.
        let cached = self.hash_cache.load(Ordering::Relaxed);
        if cached != HASH_UNCACHED {
            return cached;
        }
        let digest = structural_hash::<A, H>(&self.xbits.data, &self.zbits.data, self.nqubits);
        self.hash_cache.store(digest, Ordering::Relaxed);
        digest
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_traits_2::{IdentityBuildHasher, PauliBits};
    use std::collections::HashMap;

    #[test]
    fn equal_words_equal_digest() {
        let a: PauliWord = "XYZI".into();
        let b: PauliWord = "XYZI".into();
        assert_eq!(a.key_hash(), b.key_hash());
    }

    #[test]
    fn hash_writes_key_hash() {
        // The `Hash` impl must reproduce `key_hash()` exactly through the
        // identity build-hasher.
        let w: PauliWord = "XYZI".into();
        let bh = IdentityBuildHasher;
        assert_eq!(bh.hash_one(&w), w.key_hash());
    }

    #[test]
    fn cache_is_stable_across_clone_and_mutation() {
        let w: PauliWord = "XYZI".into();
        let h0 = w.key_hash();
        let c = w.clone();
        assert_eq!(c.key_hash(), h0, "clone copies the cached digest");

        let mut m = w.clone();
        m.set_x_bit(3, true); // I -> X on qubit 3, a structural change
        assert_ne!(m.key_hash(), h0, "mutation invalidates and recomputes");
    }

    // The `AtomicU64` hash cache is interior-mutable, but it is excluded from
    // `Eq`/`Hash` (only `(nqubits, X, Z)` participate), so the digest a stored
    // key hashes under is stable — precisely the lazy-cache pattern the design
    // sanctions (§"Lazy hashing and interior mutability"). Clippy's
    // `mutable_key_type` cannot see that exclusion.
    #[test]
    #[allow(clippy::mutable_key_type)]
    fn usable_as_identity_hashmap_key() {
        let mut map: HashMap<PauliWord, i32, IdentityBuildHasher> = HashMap::default();
        map.insert("XYZI".into(), 7);
        assert_eq!(map.get(&PauliWord::from("XYZI")), Some(&7));
    }

    #[test]
    fn avalanche_low_bits_distribute() {
        // A weak distribution property test (Design's stated contract, not a
        // type-level guarantee): enumerating single-qubit-different words, the
        // low 8 bits of the digest should not collapse into a few buckets.
        use std::collections::HashSet;
        let mut buckets = HashSet::new();
        for i in 0..8usize {
            let mut w: PauliWord<u64> = PauliWord::new(16);
            w.set_x_bit(i, true);
            buckets.insert(w.key_hash() & 0xff);
            let mut z: PauliWord<u64> = PauliWord::new(16);
            z.set_z_bit(i, true);
            buckets.insert(z.key_hash() & 0xff);
        }
        assert!(buckets.len() >= 12, "low bits collapsed: {}", buckets.len());
    }
}
