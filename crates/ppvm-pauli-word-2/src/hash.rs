// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Packed-word hashing, cached key lookup, and hash-table distribution helpers.

use std::hash::{BuildHasher, Hash, Hasher};

use ppvm_traits_2::Indexable;

use crate::{PauliStorage, PauliWord};

/// How a hasher computes a word's structural hash and its map-lookup hash.
/// `hashbrown` uses low bits for buckets and the top seven bits for control tags.
/// The defaults hash with the hasher itself; hashers can override as needed.
pub trait HashFinalize {
    /// Hashes a word given its bit planes as byte slices, each `storage_bytes`
    /// long. The default feeds them to this hasher and applies `finalize_hash`.
    #[inline(always)]
    fn word_hash(parts: &[&[u8]], storage_bytes: usize) -> u64
    where
        Self: BuildHasher + Default,
    {
        let mut hasher = Self::default().build_hasher();
        for part in parts {
            hasher.write(part);
        }
        Self::finalize_hash(hasher.finish(), storage_bytes)
    }

    /// Finalizes `Hasher::finish()` using the storage size in bytes per bit plane.
    /// Callers supply a compile-time size, allowing width-dependent branches to fold away.
    #[inline(always)]
    fn finalize_hash(raw: u64, _storage_bytes: usize) -> u64 {
        raw
    }

    /// Apply the map-index transform to a finalized structural digest.
    #[inline(always)]
    fn index_hash(raw: u64) -> u64
    where
        Self: BuildHasher + Default,
    {
        let mut hasher = Self::default().build_hasher();
        hasher.write_u64(raw);
        hasher.finish()
    }
}

impl HashFinalize for fxhash::FxBuildHasher {
    // FxHash can collapse distinct X/Z planes to the same digest.
    #[inline(always)]
    fn word_hash(parts: &[&[u8]], _storage_bytes: usize) -> u64 {
        let mut hasher = foldhash::fast::FixedState::default().build_hasher();
        for part in parts {
            hasher.write(part);
        }
        hasher.finish()
    }

    #[cfg(target_pointer_width = "64")]
    #[inline(always)]
    fn index_hash(raw: u64) -> u64 {
        raw.wrapping_mul(0x517c_c1b7_2722_0a95)
    }
}

/// Hashes the X/Z planes and applies the map-lookup transform.
/// Shared by scalar and column storage so equal words receive identical digests.
#[inline(always)]
pub(crate) fn structural_hash<A, H>(x: &A, z: &A, _nqubits: usize) -> u64
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    // Qubit count is omitted because each sum uses a common width.
    let structural = H::word_hash(
        &[bytemuck::bytes_of(x), bytemuck::bytes_of(z)],
        std::mem::size_of::<A>(),
    );
    // Direct bucket lookup requires the same transform as hashing a cached digest.
    H::index_hash(structural)
}

/// Writes the cached key hash as one `u64` for pass-through identity hashing.
impl<A, H> Hash for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline(always)]
    fn hash<S: Hasher>(&self, state: &mut S) {
        state.write_u64(self.key_hash());
    }
}

/// Returns the cached digest, refreshed eagerly by constructors and mutations.
/// Reading a key hash does not mutate the word or recompute its digest.
impl<A, H> Indexable for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline(always)]
    fn key_hash(&self) -> u64 {
        self.hash_cache
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
        assert_eq!(bh.hash_one(w), w.key_hash());
    }

    #[test]
    fn cache_is_stable_across_clone_and_mutation() {
        let w: PauliWord = "XYZI".into();
        let h0 = w.key_hash();
        let c = w;
        assert_eq!(c.key_hash(), h0, "clone copies the cached digest");

        let mut m = w;
        m.set_x_bit(3, true); // I -> X on qubit 3, a structural change
        assert_ne!(m.key_hash(), h0, "mutation refreshes the digest");
    }

    // The eager cache is a plain `u64`; stored keys have no interior mutability.
    #[test]
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

    // Every Pauli string on the top 7 of 64 qubits. Maps use the cached hash
    // directly, so its low bits must spread these words.
    #[test]
    fn top_qubit_words_spread() {
        use std::collections::HashSet;
        let buckets: HashSet<u64> = (0..1usize << 14)
            .map(|code| {
                let mut w: PauliWord<u64> = PauliWord::new(64);
                for j in 0..7 {
                    let pauli = (code >> (2 * j)) & 3;
                    w.set_x_bit(57 + j, pauli & 1 != 0);
                    w.set_z_bit(57 + j, pauli & 2 != 0);
                }
                IdentityBuildHasher.hash_one(w) & 0xFFFF
            })
            .collect();
        assert!(buckets.len() > 1 << 13, "{} buckets", buckets.len());
    }

    const RAW: u64 = 0xDEAD_BEEF_0000_0001;

    #[test]
    fn fxhash_index_fast_path_matches_hasher() {
        for raw in [0, 1, RAW, u64::MAX] {
            let mut hasher = fxhash::FxBuildHasher::default().build_hasher();
            hasher.write_u64(raw);
            assert_eq!(
                <fxhash::FxBuildHasher as HashFinalize>::index_hash(raw),
                hasher.finish()
            );
        }
    }
}
