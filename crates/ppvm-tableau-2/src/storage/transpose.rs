// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! In-place transposition of a square `n × n` bit matrix stored as `n` majors
//! of `stride` `u64` words.
//!
//! This is the orientation change [`TableauData`](super::TableauData) performs
//! when a caller needs generator-contiguous rows instead of qubit-contiguous
//! columns. It is the direct analogue of Stim's `do_transpose_quadrants`
//! (`stim/stabilizers/tableau.inl`), and it is the reason the tableau is stored
//! as **four square quadrants** rather than two `2n × n` rectangles: a
//! rectangular bit-matrix transpose is not an in-place permutation, a square one
//! is.
//!
//! The unit of work is a `64 × 64` bit block, transposed by the standard
//! recursive shift-mask algorithm (Hacker's Delight §7-3). Blocks on the
//! diagonal are transposed in place; off-diagonal blocks `(bi, bj)` with
//! `bi < bj` are transposed and swapped with `(bj, bi)`.
//!
//! Ragged edges (`n` not a multiple of 64) are handled by gathering each block
//! into a zero-filled scratch and scattering it back, which keeps **padding
//! bits zero by construction** — the invariant the bulk equality and hashing
//! paths depend on.
//!
//! A block that is only populated in its top-left `e × e` corner — the whole
//! matrix when `n < 64`, the trailing diagonal block otherwise — is transposed
//! over a **span** of `e.next_power_of_two()` rows rather than all 64. See
//! [`transpose_block`] for why that is exact.

/// Bits per machine word, and the side length of one transpose block.
pub(crate) const BITS_PER_WORD: usize = 64;

/// The shift-mask round masks, indexed by `log₂ j`: `j` ones alternating with
/// `j` zeros.
const ROUND_MASKS: [u64; 6] = [
    0x5555_5555_5555_5555,
    0x3333_3333_3333_3333,
    0x0f0f_0f0f_0f0f_0f0f,
    0x00ff_00ff_00ff_00ff,
    0x0000_ffff_0000_ffff,
    0x0000_0000_ffff_ffff,
];

/// Transpose the top-left `span × span` corner of a bit block held as rows of
/// one `u64` each, where the span is `block.len()`, a power of two, and every
/// bit outside that corner is zero.
///
/// Row `r`, column `c` is bit `c` of `block[r]`. Afterwards row `r`, column `c`
/// is the original row `c`, column `r`.
///
/// Truncating to `span` is exact rather than approximate. A round `j` exchanges
/// the off-diagonal `j × j` sub-blocks of each `2j × 2j` block, so with `j ≥
/// span` both operands are outside the populated corner and the round is the
/// identity; and every round with `j < span` keeps its exchanges inside a
/// `2j ≤ span` block, so the corner is closed under the rounds that remain.
/// Rows and columns at index `≥ span` are therefore zero on exit as well as on
/// entry — the padding invariant this module owes its callers. (Rows between
/// the populated extent `e` and `span` are *not* zero throughout: the rounds
/// stage bits there and take them back before the last one returns.)
#[inline]
fn transpose_block(block: &mut [u64]) {
    let span = block.len();
    debug_assert!(span.is_power_of_two() && span <= BITS_PER_WORD);
    // Split so the full-block case reaches `shift_mask_rounds` with a literal
    // round count and stays fully unrolled, as it was before the span existed.
    if span == BITS_PER_WORD {
        shift_mask_rounds(block, 6);
    } else {
        shift_mask_rounds(block, span.trailing_zeros() as usize);
    }
}

/// The `log_span` shift-mask rounds, largest `j` first.
#[inline(always)]
fn shift_mask_rounds(block: &mut [u64], log_span: usize) {
    for round in (0..log_span).rev() {
        let (j, m) = (1usize << round, ROUND_MASKS[round]);
        let mut k = 0usize;
        while k < block.len() {
            // Exchange the off-diagonal `j × j` sub-blocks of the `2j × 2j`
            // block whose top-left corner is row `k`.
            let t = ((block[k] >> j) ^ block[k + j]) & m;
            block[k + j] ^= t;
            block[k] ^= t << j;
            k = (k + j + 1) & !j;
        }
    }
}

/// Read the block whose top-left corner is `(row0, col0)` into a zero-filled
/// scratch of `out.len()` rows.
///
/// Rows and columns past `n` read as zero, so the scratch always describes a
/// full span and the transpose never has to special-case the ragged edge.
#[inline]
fn gather(words: &[u64], stride: usize, n: usize, row0: usize, col0: usize, out: &mut [u64]) {
    let word = col0 / BITS_PER_WORD;
    let rows = (n - row0).min(out.len());
    out[rows..].fill(0);
    for (r, slot) in out[..rows].iter_mut().enumerate() {
        *slot = words[(row0 + r) * stride + word];
    }
}

/// Write a transposed block back at `(row0, col0)`.
///
/// Rows past `n` are dropped; the bits they would have carried are zero because
/// [`gather`] zero-filled them.
#[inline]
fn scatter(words: &mut [u64], stride: usize, n: usize, row0: usize, col0: usize, src: &[u64]) {
    let word = col0 / BITS_PER_WORD;
    let rows = (n - row0).min(src.len());
    for (r, &value) in src[..rows].iter().enumerate() {
        words[(row0 + r) * stride + word] = value;
    }
}

/// Transpose the square `n × n` bit matrix held in `words` in place.
///
/// `words` is `n` majors of `stride` words each, `stride >= n.div_ceil(64)`.
/// Bits at index `>= n` inside a major, and words past `n.div_ceil(64)`, are
/// required to be zero on entry and are left zero on exit.
///
/// Only the trailing diagonal block can be partially populated — an
/// off-diagonal block `(bi, bj)` has `bi < bj ≤ blocks - 1`, so its 64 rows are
/// all inside `n` — which is what lets a sub-64 matrix skip almost the whole
/// shift-mask ladder while a large one pays exactly what it did before.
pub(crate) fn transpose_square(words: &mut [u64], stride: usize, n: usize) {
    debug_assert!(stride >= n.div_ceil(BITS_PER_WORD));
    debug_assert!(words.len() >= n * stride);
    let blocks = n.div_ceil(BITS_PER_WORD);

    let mut a = [0u64; BITS_PER_WORD];
    let mut b = [0u64; BITS_PER_WORD];
    for bi in 0..blocks {
        let row0 = bi * BITS_PER_WORD;
        let span = (n - row0).min(BITS_PER_WORD).next_power_of_two();

        // Diagonal block: transpose in place.
        gather(words, stride, n, row0, row0, &mut a[..span]);
        transpose_block(&mut a[..span]);
        scatter(words, stride, n, row0, row0, &a[..span]);

        // Off-diagonal blocks: transpose both and swap.
        for bj in (bi + 1)..blocks {
            let col0 = bj * BITS_PER_WORD;
            gather(words, stride, n, row0, col0, &mut a);
            gather(words, stride, n, col0, row0, &mut b);
            transpose_block(&mut a);
            transpose_block(&mut b);
            scatter(words, stride, n, row0, col0, &b);
            scatter(words, stride, n, col0, row0, &a);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Bit `(r, c)` of a matrix held as `n` majors of `stride` words.
    fn get(words: &[u64], stride: usize, r: usize, c: usize) -> bool {
        words[r * stride + c / BITS_PER_WORD] >> (c % BITS_PER_WORD) & 1 == 1
    }

    fn set(words: &mut [u64], stride: usize, r: usize, c: usize, v: bool) {
        let w = &mut words[r * stride + c / BITS_PER_WORD];
        let mask = 1u64 << (c % BITS_PER_WORD);
        if v {
            *w |= mask;
        } else {
            *w &= !mask;
        }
    }

    /// A cheap reproducible bit pattern that is *not* symmetric, so a transpose
    /// that silently does nothing cannot pass.
    fn pattern(r: usize, c: usize) -> bool {
        (r * 7 + c * 13 + r * c * 3).is_multiple_of(5)
    }

    #[test]
    fn transpose_block_matches_naive() {
        let mut block = [0u64; 64];
        for (r, slot) in block.iter_mut().enumerate() {
            for c in 0..64 {
                if pattern(r, c) {
                    *slot |= 1u64 << c;
                }
            }
        }
        let original = block;
        transpose_block(&mut block);
        for r in 0..64 {
            for c in 0..64 {
                let got = block[r] >> c & 1 == 1;
                let want = original[c] >> r & 1 == 1;
                assert_eq!(got, want, "bit ({r}, {c})");
            }
        }
    }

    /// The truncated span: a block populated only in its top-left `e × e`
    /// corner transposes that corner and leaves everything outside it zero.
    #[test]
    fn transpose_block_span_matches_naive_and_keeps_padding_zero() {
        for e in 1usize..=64 {
            let span = e.next_power_of_two();
            let mut block = [0u64; 64];
            for (r, slot) in block[..e].iter_mut().enumerate() {
                for c in 0..e {
                    if pattern(r, c) {
                        *slot |= 1u64 << c;
                    }
                }
            }
            let original = block;
            transpose_block(&mut block[..span]);
            for (r, &row) in block.iter().enumerate() {
                for (c, &col) in original.iter().enumerate() {
                    let want = r < e && c < e && col >> r & 1 == 1;
                    assert_eq!(row >> c & 1 == 1, want, "e={e} bit ({r}, {c})");
                }
            }
        }
    }

    #[test]
    fn transpose_square_matches_naive_and_is_an_involution() {
        // Every width below one block, then either side of the block boundary,
        // including ragged ones.
        for &n in &[
            1usize, 2, 3, 4, 5, 6, 7, 8, 63, 64, 65, 100, 127, 128, 129, 200,
        ] {
            let stride = n.div_ceil(BITS_PER_WORD).next_multiple_of(4).max(4);
            let mut words = vec![0u64; n * stride];
            for r in 0..n {
                for c in 0..n {
                    set(&mut words, stride, r, c, pattern(r, c));
                }
            }
            let original = words.clone();

            transpose_square(&mut words, stride, n);
            for r in 0..n {
                for c in 0..n {
                    assert_eq!(
                        get(&words, stride, r, c),
                        get(&original, stride, c, r),
                        "n={n} bit ({r}, {c})"
                    );
                }
            }

            transpose_square(&mut words, stride, n);
            assert_eq!(words, original, "n={n}: transpose is not an involution");
        }
    }

    #[test]
    fn transpose_leaves_padding_zero() {
        for &n in &[1usize, 2, 3, 4, 5, 6, 7, 8, 63, 65, 100, 129] {
            let stride = n.div_ceil(BITS_PER_WORD).next_multiple_of(4).max(4);
            let mut words = vec![0u64; n * stride];
            for r in 0..n {
                for c in 0..n {
                    set(&mut words, stride, r, c, pattern(r, c));
                }
            }
            transpose_square(&mut words, stride, n);
            for r in 0..n {
                for c in n..stride * BITS_PER_WORD {
                    assert!(!get(&words, stride, r, c), "n={n}: padding bit ({r}, {c})");
                }
            }
        }
    }
}
