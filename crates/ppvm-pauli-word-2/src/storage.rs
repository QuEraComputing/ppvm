// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bitvec::view::{BitView, BitViewSized};
use num::PrimInt;
use std::hash::Hash;

/// Default storage is `u64`, or `usize` on wasm32 where `bitvec` lacks `u64` support.
#[cfg(not(target_arch = "wasm32"))]
pub type DefaultStorage = u64;

/// wasm32 default storage; `bitvec` requires a supported native integer type.
#[cfg(target_arch = "wasm32")]
pub type DefaultStorage = usize;

/// Fixed-size storage for one Pauli bit plane: `u64`, `[u8; N]`, or `[u64; N]`.
/// `Pod` permits safe byte views for hashing; `PrimInt` enables whole-integer kernels.
pub trait PauliStorage:
    BitViewSized
    + BitView<Store: PrimInt>
    + Copy
    + Hash
    + Eq
    + Send
    + Sync
    + std::fmt::Debug
    + bytemuck::Pod
{
}

impl<A> PauliStorage for A where
    A: BitViewSized
        + BitView<Store: PrimInt>
        + Copy
        + Hash
        + Eq
        + Send
        + Sync
        + std::fmt::Debug
        + bytemuck::Pod
{
}
