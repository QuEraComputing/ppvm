// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Structural key digests and a pass-through hasher.
//! [`Indexable`] is required by hash backends, not by all container keys.

use std::hash::{BuildHasher, Hash, Hasher};

/// A key with an avalanche-quality structural digest; equal keys have equal digests.
/// `Hash` must write exactly one `u64`: `state.write_u64(self.key_hash())`.
/// Column hashing must reproduce that digest bit for bit.
pub trait Indexable: Clone + Eq + Hash {
    /// The finalized structural digest of this key.
    fn key_hash(&self) -> u64;
}

/// A pass-through hasher returning the single `u64` digest written by a key.
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

/// Builds [`IdentityHasher`] instances so map hashes equal the supplied key digests.
#[derive(Debug, Default, Clone)]
pub struct IdentityBuildHasher;

impl BuildHasher for IdentityBuildHasher {
    type Hasher = IdentityHasher;

    #[inline]
    fn build_hasher(&self) -> IdentityHasher {
        IdentityHasher::default()
    }
}
