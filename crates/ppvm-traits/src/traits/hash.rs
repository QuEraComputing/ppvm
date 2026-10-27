// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Hashers for internal maps and cached Pauli-word hashes.

use std::collections::HashMap;
use std::hash::{BuildHasher, Hasher};

/// Fast hasher for internal maps and word hashes, with a fixed seed.
/// Iteration is reproducible for the same build and insertion order.
pub type FoldBuildHasher = foldhash::fast::FixedState;

/// A [`HashMap`] using [`FoldBuildHasher`].
pub type FoldHashMap<K, V> = HashMap<K, V, FoldBuildHasher>;

/// Computes and finalizes cached Pauli-word hashes.
/// Custom hashers can use the defaults with an empty implementation.
pub trait HashFinalize {
    /// Hash a Pauli word given its bit arrays as byte slices (X, Z, and for
    /// lossy words the loss bits), each `storage_bytes` long.
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

    /// Adjust the raw digest computed by the default
    /// [`word_hash`](Self::word_hash). The default returns it unchanged.
    #[inline(always)]
    fn finalize_hash(raw: u64, _storage_bytes: usize) -> u64 {
        raw
    }
}

impl HashFinalize for fxhash::FxBuildHasher {
    // FxHash can collapse distinct X/Z planes to the same digest; use foldhash
    // for the word's cached hash while retaining FxHash for map lookup.
    #[inline(always)]
    fn word_hash(parts: &[&[u8]], _storage_bytes: usize) -> u64 {
        let mut hasher = FoldBuildHasher::default().build_hasher();
        for part in parts {
            hasher.write(part);
        }
        hasher.finish()
    }
}

#[cfg(all(
    feature = "gxhash",
    target_feature = "aes",
    any(target_arch = "x86_64", target_arch = "aarch64")
))]
impl HashFinalize for gxhash::GxBuildHasher {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    const KEYS: u64 = 1 << 14;

    /// Distinct values among the low 16 bits, which a hash map uses as the bucket.
    fn buckets(hashes: impl Iterator<Item = u64>) -> usize {
        hashes.map(|h| h & 0xFFFF).collect::<HashSet<_>>().len()
    }

    #[test]
    fn fold_build_hasher_spreads_high_bit_keys() {
        let fold = FoldBuildHasher::default();
        let usize_keys = (0..KEYS).map(|k| fold.hash_one((k as usize) << (usize::BITS - 14)));
        let u128_keys = (0..KEYS).map(|k| fold.hash_one(u128::from(k) << 114));
        assert!(buckets(usize_keys) > KEYS as usize / 2);
        assert!(buckets(u128_keys) > KEYS as usize / 2);
    }

    // Every Pauli string on the top 7 qubits, hashed the way a map with FxHash
    // sees a word: FxHash over the cached word hash.
    #[test]
    fn fxhash_word_hash_spreads_top_qubit_words() {
        let fx = fxhash::FxBuildHasher::default();
        for bytes in [8, 16, 64] {
            let n_qubits = 8 * bytes;
            let hashes = (0..KEYS).map(|code| {
                let (mut x, mut z) = (vec![0u8; bytes], vec![0u8; bytes]);
                for j in 0..7 {
                    let (q, pauli) = (n_qubits - 7 + j, (code >> (2 * j)) & 3);
                    x[q / 8] |= ((pauli & 1) as u8) << (q % 8);
                    z[q / 8] |= ((pauli >> 1) as u8) << (q % 8);
                }
                let cached = <fxhash::FxBuildHasher as HashFinalize>::word_hash(&[&x, &z], bytes);
                fx.hash_one(cached)
            });
            assert!(buckets(hashes) > KEYS as usize / 2, "{n_qubits} qubits");
        }
    }
}
