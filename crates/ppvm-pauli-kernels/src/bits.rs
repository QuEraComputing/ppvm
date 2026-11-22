// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The value types the kernels run on.

use std::ops::{BitAnd, BitOr, BitXor, Not};

/// A value supporting the Boolean operations the kernels use.
///
/// `bool` holds one Pauli's bit; an integer holds one bit per position, each
/// position an independent Pauli. The operations act bitwise, so a rule on an
/// integer is the rule applied to each bit position.
pub trait Bits:
    Copy + BitAnd<Output = Self> + BitOr<Output = Self> + BitXor<Output = Self> + Not<Output = Self>
{
    /// The number of set bits.
    fn popcount(self) -> u32; // needs its own method since no rust trait guarantees this; can't use num::PrimInt because of bool
}

impl Bits for bool {
    #[inline(always)]
    fn popcount(self) -> u32 {
        self as u32
    }
}

macro_rules! impl_int {
    ($($t:ty),*) => {$(
        impl Bits for $t {
            #[inline(always)]
            fn popcount(self) -> u32 {
                self.count_ones()
            }
        }
    )*};
}

impl_int!(u8, u16, u32, u64, u128, usize);
