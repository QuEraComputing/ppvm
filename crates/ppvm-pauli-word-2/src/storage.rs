// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The private representation parameters of a packed [`crate::PauliWord`]: the
//! backing-storage blob bound [`PauliStorage`] and the per-hasher, per-width
//! finalization fold [`HashFinalize`].
//!
//! Design: `word-data-structures.md` §"`PauliWord` packed representation" (`A`
//! and `H` are private representation parameters, on the same footing, neither
//! exposed through `Word`/`Indexable`) and `traits-2-configuration-and-hashing.md`
//! §"Concrete word hashing" (the finalization fold is a **private utility inside
//! this crate**, not part of the algebra-agnostic `Indexable` contract). Both
//! are ported from `ppvm-pauli-word`'s `PauliStorage` / `HashFinalize` to keep
//! the hot paths and the bucket distribution at parity.

use bitvec::view::BitViewSized;
use std::hash::{BuildHasher, Hash, Hasher};

/// Native-word default storage: `u64` on native targets and `usize` on wasm32.
///
/// `bitvec` implements `BitStore` for `u64` only on 64-bit pointer-width
/// targets; browser wasm is 32-bit and therefore uses its native `usize`.
#[cfg(not(target_arch = "wasm32"))]
pub type DefaultStorage = u64;
/// See the native definition above.
#[cfg(target_arch = "wasm32")]
pub type DefaultStorage = usize;

/// Backing storage for a [`crate::PauliWord`] — a fixed-size, `Copy`-able block
/// of bits (typically `u64`, `[u8; N]`, or `[u64; N]`).
///
/// The [`bytemuck::Pod`] bound guarantees plain-old-data with no padding and all
/// bit patterns valid, which lets the structural hash view the blob as a `&[u8]`
/// without `unsafe`. Design: `word-data-structures.md` §"`PauliWord` packed
/// representation" (the `A` parameter).
pub trait PauliStorage:
    BitViewSized + Clone + Copy + Hash + Eq + Send + Sync + std::fmt::Debug + bytemuck::Pod
{
}

impl<A> PauliStorage for A where
    A: BitViewSized + Clone + Copy + Hash + Eq + Send + Sync + std::fmt::Debug + bytemuck::Pod
{
}

/// Per-hasher, per-width finalization fold applied to a word's raw digest before
/// it becomes the finalized [`ppvm_traits_2::Indexable::key_hash`] value.
///
/// The cached digest is split by `hashbrown` into a bucket index (low bits) and
/// a control tag (top 7 bits); whether the raw digest is good enough for that
/// split is a property of the **hasher**, not of the Pauli word. A weak-but-fast
/// hasher on a short key must fold; a strong hasher folds nothing (the identity
/// default). Ported from `ppvm-pauli-word`'s `HashFinalize`; kept a crate-private
/// utility per `traits-2-configuration-and-hashing.md` §"Concrete word hashing".
pub trait HashFinalize {
    /// Finalize `raw` (from `Hasher::finish`) for a key whose backing storage is
    /// `storage_bytes` wide per bit-array. The width is a compile-time constant
    /// at every call site (`size_of` of the storage), so any branch on it is
    /// monomorphized away.
    #[inline(always)]
    fn finalize_hash(raw: u64, storage_bytes: usize) -> u64 {
        let _ = storage_bytes;
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
    /// `FxHasher` avalanches weakly for short inputs: a word fitting in a single
    /// `u64` per bit-array (`[u8; 8]` and narrower) goes through only a couple of
    /// multiply-rotate rounds, leaving the low bits — the ones `hashbrown` uses
    /// to choose a bucket — correlated. Folding the high half into the low half
    /// decorrelates them. Wider storage already distributes its low bits, so it
    /// passes through (folding there would couple the tag back into the bucket).
    #[inline(always)]
    fn finalize_hash(raw: u64, storage_bytes: usize) -> u64 {
        if storage_bytes <= std::mem::size_of::<u64>() {
            raw ^ (raw >> 32)
        } else {
            raw
        }
    }

    #[cfg(target_pointer_width = "64")]
    #[inline(always)]
    fn index_hash(raw: u64) -> u64 {
        raw.wrapping_mul(0x517c_c1b7_2722_0a95)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const RAW: u64 = 0xDEAD_BEEF_0000_0001;

    #[test]
    fn fxhash_folds_narrow_storage() {
        for width in [1, 2, 4, 8] {
            assert_eq!(
                <fxhash::FxBuildHasher as HashFinalize>::finalize_hash(RAW, width),
                RAW ^ (RAW >> 32),
            );
        }
    }

    #[test]
    fn fxhash_passes_wide_storage_through() {
        for width in [16, 32, 64] {
            assert_eq!(
                <fxhash::FxBuildHasher as HashFinalize>::finalize_hash(RAW, width),
                RAW,
            );
        }
    }

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
