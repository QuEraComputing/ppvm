// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bitvec::view::BitViewSized;
use std::hash::Hash;

/// Backing storage for a [`PauliWord`](crate::word::PauliWord) — a
/// fixed-size, `Copy`-able block of bits (typically `[u8; N]` or
/// `[u64; N]`).
///
/// The [`bytemuck::Pod`] bound guarantees the storage is a plain-old-data
/// type with no padding and all bit patterns valid, which lets
/// [`PauliWord::rehash`](crate::word::PauliWord) view it as a `&[u8]` for
/// fast byte-slice hashing *without* `unsafe`. Every concrete storage used
/// here (`[u8; N]`, `[u64; N]`) is already `Pod`.
pub trait PauliStorage:
    BitViewSized
    + Clone
    + Copy
    + Hash
    + Eq
    + PartialEq
    + Send
    + Sync
    + std::fmt::Debug
    + Ord
    + PartialOrd
    + bytemuck::Pod
{
    /// Finalize the raw `FxHasher` digest of a Pauli word's bytes into the
    /// value cached in [`PauliWord`](crate::word::PauliWord).
    ///
    /// This is the single place that adapts the hash to the storage width,
    /// and it is a compile-time decision: `size_of::<Self>()` is a constant
    /// per monomorphization, so the branch below folds away to either a bare
    /// `finish()` or `finish()` plus one shift/xor — no runtime dispatch.
    ///
    /// **Why width matters.** `FxHasher` consumes the bytes one `u64` word at
    /// a time, and its avalanche is weak for short inputs: a word backed by a
    /// single `u64` per bit-array (`[u8; 8]` and narrower) goes through only a
    /// couple of multiply-rotate rounds, leaving the *low* bits — exactly the
    /// ones `hashbrown` uses to choose a bucket — highly correlated. At high
    /// fill that clusters distinct words into a few oversized buckets and the
    /// `insert_unique` probe length explodes (measured ~7x at 64 qubits in
    /// `[u8; 8]`; see `examples/trotter_storage_cliff.rs`). Folding the high
    /// half into the low half (`h ^ (h >> 32)`) decorrelates them and removes
    /// the cliff.
    ///
    /// For wider storage (`[u8; 16]` and up) the extra `u64` rounds already
    /// distribute the low bits well. Folding there is not just unnecessary but
    /// mildly *harmful* — it mixes the top bits, which `hashbrown` reserves for
    /// its control-byte tag, back into the bucket index, coupling two values
    /// it wants independent (~4–6% slower in the same benchmark). So we only
    /// fold for storage that fits in a single `u64` word per bit-array and
    /// pass the digest through unchanged otherwise.
    #[inline(always)]
    fn finalize_hash(raw: u64) -> u64 {
        if std::mem::size_of::<Self>() <= std::mem::size_of::<u64>() {
            raw ^ (raw >> 32)
        } else {
            raw
        }
    }
}
impl<
    A: BitViewSized
        + Clone
        + Copy
        + Hash
        + Eq
        + PartialEq
        + Send
        + Sync
        + std::fmt::Debug
        + Ord
        + PartialOrd
        + bytemuck::Pod,
> PauliStorage for A
{
}
