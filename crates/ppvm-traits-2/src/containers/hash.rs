// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Indexable keys and the identity pass-through hasher.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Indexable values" and
//! §"The pass-through storage contract". `Indexable` is *not* the universal key
//! bound (that is `Eq + Clone`); it is required only on the hash backends.

use std::hash::{BuildHasher, Hash, Hasher};

/// A key whose finalized structural digest is first class.
///
/// The digest is avalanche-quality in both the low bits (the hashbrown bucket)
/// and the top 7 (the control tag), so it can be consumed *directly* as the map
/// hash. Contracts:
///
///   * `Hash for Self` is exactly `state.write_u64(self.key_hash())`;
///   * structurally equal keys return equal digests; and
///   * `KeyColumn::hash_into` reproduces this value bit for bit.
///
/// This exposes the digest *value*, not the cache mechanics — there is no cache
/// type or invalidation hook in the contract.
///
/// Design: §"Indexable values".
pub trait Indexable: Clone + Eq + Hash {
    /// The finalized structural digest of this key.
    fn key_hash(&self) -> u64;
}

/// A pass-through `Hasher`: a key writes its already-finalized `key_hash()` as a
/// single `u64` and this hands it back verbatim, so the digest reaches
/// hashbrown untouched.
///
/// Design: §"The pass-through storage contract".
#[derive(Debug, Default, Clone)]
pub struct IdentityHasher(u64);

impl Hasher for IdentityHasher {
    #[inline]
    fn write_u64(&mut self, n: u64) {
        self.0 = n; // store the digest
    }

    fn write(&mut self, _: &[u8]) {
        unreachable!("Indexable keys write exactly one u64 (their key_hash())")
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.0 // hand it back verbatim
    }
}

/// `BuildHasher` for [`IdentityHasher`]; the storage aliases in
/// `ppvm-pauli-sum-2` bake this into their `HashMap` so `finish() ==
/// key.key_hash()`.
///
/// Design: §"The pass-through storage contract".
#[derive(Debug, Default, Clone)]
pub struct IdentityBuildHasher;

impl BuildHasher for IdentityBuildHasher {
    type Hasher = IdentityHasher;

    #[inline]
    fn build_hasher(&self) -> IdentityHasher {
        IdentityHasher::default()
    }
}
