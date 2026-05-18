// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

use num::PrimInt;
use std::hash::Hash;
use std::ops::{BitAnd, BitOrAssign, BitXor, Shl};

/// Bit-string index type used by
/// [`GeneralizedTableau`](crate::data::GeneralizedTableau) to address
/// individual branches of its sparse coefficient vector.
///
/// Blanket-implemented for every primitive (and `bnum`-style) integer
/// type that supports the required bit operations. Pick:
/// * `usize` for ≤ 64 qubits,
/// * `u128` for ≤ 128,
/// * `bnum::types::U256` / `U512` / `U1024` for the wide regime.
pub trait TableauIndex:
    PartialEq
    + Eq
    + Hash
    + Copy
    + From<u8>
    + Shl<usize, Output = Self>
    + BitOrAssign<Self>
    + BitAnd<Self, Output = Self>
    + BitXor<Output = Self>
    + PrimInt
{
}

impl<I> TableauIndex for I where
    I: PartialEq
        + Eq
        + Hash
        + Copy
        + From<u8>
        + Shl<usize, Output = Self>
        + BitOrAssign<Self>
        + BitAnd<I, Output = I>
        + BitXor<Output = I>
        + PrimInt
{
}
