// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod clifford;
mod data;
mod mul;
mod trace;

pub use data::PauliWord;

/// Bitwise-OR of `K` equal-length byte slices, popcounted as a single
/// total. Used by `PauliWord::weight`, `LossyPauliWord::weight`, and
/// `LossyPauliWord::loss_weight` to bypass `bitvec`'s `Domain` /
/// per-element abstractions, which block LLVM from inlining the native
/// popcount.
///
/// Slices are processed in u64 / u32 / u16 / u8 chunks; all branches
/// const-fold at each call site once `K` and slice length are known
/// after monomorphisation.
#[inline]
pub(crate) fn or_popcount<const K: usize>(parts: [&[u8]; K]) -> usize {
    let n = parts[0].len();
    let mut total: u32 = 0;
    let mut i = 0usize;

    while i + 8 <= n {
        let mut acc: u64 = 0;
        for p in parts.iter() {
            acc |= u64::from_ne_bytes(p[i..i + 8].try_into().unwrap());
        }
        total += acc.count_ones();
        i += 8;
    }
    if i + 4 <= n {
        let mut acc: u32 = 0;
        for p in parts.iter() {
            acc |= u32::from_ne_bytes(p[i..i + 4].try_into().unwrap());
        }
        total += acc.count_ones();
        i += 4;
    }
    if i + 2 <= n {
        let mut acc: u16 = 0;
        for p in parts.iter() {
            acc |= u16::from_ne_bytes(p[i..i + 2].try_into().unwrap());
        }
        total += acc.count_ones();
        i += 2;
    }
    if i < n {
        let mut acc: u8 = 0;
        for p in parts.iter() {
            acc |= p[i];
        }
        total += acc.count_ones();
    }

    total as usize
}
