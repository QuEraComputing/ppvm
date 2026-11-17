// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::ops::{BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign, Not};

use rand::Rng;

pub(crate) const LANES: usize = 8;
/// Shots simulated together, one bit each.
pub(crate) const SHOTS: usize = 64 * LANES;

/// One bit per shot of a batch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Bits(pub [u64; LANES]);

impl Bits {
    pub const ZERO: Self = Self([0; LANES]);
    pub const ONES: Self = Self([!0; LANES]);

    pub fn random(rng: &mut impl Rng) -> Self {
        Self(std::array::from_fn(|_| rng.next_u64()))
    }

    pub fn get(&self, shot: usize) -> bool {
        self.0[shot / 64] >> (shot % 64) & 1 == 1
    }

    pub fn set(&mut self, shot: usize, value: bool) {
        let mask = 1 << (shot % 64);
        let word = &mut self.0[shot / 64];
        *word = if value { *word | mask } else { *word & !mask };
    }

    pub fn flip(&mut self, shot: usize) {
        self.0[shot / 64] ^= 1 << (shot % 64);
    }
}

macro_rules! bitwise {
    ($trait:ident, $method:ident, $assign_trait:ident, $assign_method:ident) => {
        impl $trait for Bits {
            type Output = Self;
            #[inline]
            fn $method(self, rhs: Self) -> Self {
                Self(std::array::from_fn(|i| self.0[i].$method(rhs.0[i])))
            }
        }

        impl $assign_trait for Bits {
            #[inline]
            fn $assign_method(&mut self, rhs: Self) {
                for (a, b) in self.0.iter_mut().zip(rhs.0) {
                    a.$assign_method(b);
                }
            }
        }
    };
}

bitwise!(BitAnd, bitand, BitAndAssign, bitand_assign);
bitwise!(BitOr, bitor, BitOrAssign, bitor_assign);
bitwise!(BitXor, bitxor, BitXorAssign, bitxor_assign);

impl Not for Bits {
    type Output = Self;
    #[inline]
    fn not(self) -> Self {
        Self(self.0.map(|w| !w))
    }
}
