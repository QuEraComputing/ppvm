// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The storage word and the only code that knows its width. Everything else
//! touches words through these items and `Bits`, so a wider word can be tried
//! by changing the alias.

/// Storage word. Swap to benchmark other widths.
pub(super) type BitWord = u64;
pub(super) const BITWORD_BITS: usize = BitWord::BITS as usize;
pub(super) const ZERO: BitWord = 0;

/// Bit `i` of a packed bit string.
#[inline(always)]
pub(super) fn get_bit(words: &[BitWord], i: usize) -> bool {
    (words[i / BITWORD_BITS] >> (i % BITWORD_BITS)) & 1 == 1
}

#[inline(always)]
pub(super) fn flip_bit(words: &mut [BitWord], i: usize) {
    words[i / BITWORD_BITS] ^= 1 << (i % BITWORD_BITS);
}

/// Word index and mask of bit `i`.
#[inline(always)]
pub(super) fn bit_position(i: usize) -> (usize, BitWord) {
    (i / BITWORD_BITS, 1 << (i % BITWORD_BITS))
}

/// All ones if `b`, else zero (for branch-free masked updates).
#[inline(always)]
pub(super) fn mask_if(b: bool) -> BitWord {
    ZERO.wrapping_sub(BitWord::from(b))
}
