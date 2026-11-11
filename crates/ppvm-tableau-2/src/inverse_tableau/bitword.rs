// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The storage word and the only code that knows its width. Everything else
//! touches words through these items and `Bits`, so a wider word can be tried
//! by changing the alias.

/// Storage word. Swap to benchmark other widths.
pub(super) type BitWord = u64;
pub(super) const BITWORD_BITS: usize = BitWord::BITS as usize;
pub(super) const ZERO: BitWord = 0;
/// The highest bit of a word.
pub(super) const TOP_BIT: BitWord = 1 << (BITWORD_BITS - 1);

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

/// Indices of the set bits of a packed bit string, ascending.
pub(super) fn set_bits(words: &[BitWord]) -> impl Iterator<Item = usize> + '_ {
    words.iter().enumerate().flat_map(|(w, &word)| {
        let mut rest = word;
        std::iter::from_fn(move || {
            (rest != ZERO).then(|| {
                let b = rest.trailing_zeros() as usize;
                rest &= rest - 1;
                w * BITWORD_BITS + b
            })
        })
    })
}

/// Transpose a 64 × 64 bit matrix in place: bit `j` of word `i` swaps with bit
/// `i` of word `j`. Width-specific (64-bit words).
#[inline(always)]
pub(super) fn transpose_64x64(a: &mut [BitWord; 64]) {
    const PASSES: [(BitWord, usize); 6] = [
        (0x5555_5555_5555_5555, 1),
        (0x3333_3333_3333_3333, 2),
        (0x0F0F_0F0F_0F0F_0F0F, 4),
        (0x00FF_00FF_00FF_00FF, 8),
        (0x0000_FFFF_0000_FFFF, 16),
        (0x0000_0000_FFFF_FFFF, 32),
    ];
    for (mask, shift) in PASSES {
        for k in 0..64 {
            if k & shift == 0 {
                let (x, y) = (a[k], a[k + shift]);
                a[k] = (x & mask) | ((y & mask) << shift);
                a[k + shift] = ((x & !mask) >> shift) | (y & !mask);
            }
        }
    }
}
