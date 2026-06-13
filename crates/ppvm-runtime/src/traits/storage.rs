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
