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
//! into a zero-filled `[u64; 64]` scratch and scattering it back, which keeps
//! **padding bits zero by construction** — the invariant the bulk equality and
//! hashing paths depend on.

/// Bits per machine word, and the side length of one transpose block.
pub(crate) const BITS_PER_WORD: usize = 64;

/// Transpose one `64 × 64` bit block held as 64 rows of one `u64` each.
///
/// Row `r`, column `c` is bit `c` of `block[r]`. Afterwards row `r`, column `c`
/// is the original row `c`, column `r`.
#[inline]
fn transpose_block(block: &mut [u64; BITS_PER_WORD]) {
    let mut j = 32usize;
    let mut m = 0x0000_0000_ffff_ffffu64;
    while j != 0 {
        let mut k = 0usize;
        while k < BITS_PER_WORD {
            // Exchange the off-diagonal `j × j` sub-blocks of the `2j × 2j`
            // block whose top-left corner is row `k`.
            let t = ((block[k] >> j) ^ block[k + j]) & m;
            block[k + j] ^= t;
            block[k] ^= t << j;
            k = (k + j + 1) & !j;
        }
        j >>= 1;
        m ^= m << j;
    }
}

/// Read the `64 × 64` block whose top-left corner is `(row0, col0)` into a
/// zero-filled scratch.
///
/// Rows and columns past `n` read as zero, so the scratch always describes a
/// full block and the transpose never has to special-case the ragged edge.
#[inline]
fn gather(
    words: &[u64],
    stride: usize,
    n: usize,
    row0: usize,
    col0: usize,
    out: &mut [u64; BITS_PER_WORD],
) {
    let word = col0 / BITS_PER_WORD;
    let rows = (n - row0).min(BITS_PER_WORD);
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
fn scatter(
    words: &mut [u64],
    stride: usize,
    n: usize,
    row0: usize,
    col0: usize,
    src: &[u64; BITS_PER_WORD],
) {
    let word = col0 / BITS_PER_WORD;
    let rows = (n - row0).min(BITS_PER_WORD);
    for (r, &value) in src[..rows].iter().enumerate() {
        words[(row0 + r) * stride + word] = value;
    }
}

/// Transpose the square `n × n` bit matrix held in `words` in place.
///
/// `words` is `n` majors of `stride` words each, `stride >= n.div_ceil(64)`.
/// Bits at index `>= n` inside a major, and words past `n.div_ceil(64)`, are
/// required to be zero on entry and are left zero on exit.
pub(crate) fn transpose_square(words: &mut [u64], stride: usize, n: usize) {
    debug_assert!(stride >= n.div_ceil(BITS_PER_WORD));
    debug_assert!(words.len() >= n * stride);
    let blocks = n.div_ceil(BITS_PER_WORD);

    let mut a = [0u64; BITS_PER_WORD];
    let mut b = [0u64; BITS_PER_WORD];
    for bi in 0..blocks {
        let row0 = bi * BITS_PER_WORD;

        // Diagonal block: transpose in place.
        gather(words, stride, n, row0, row0, &mut a);
        transpose_block(&mut a);
        scatter(words, stride, n, row0, row0, &a);

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

    #[test]
    fn transpose_square_matches_naive_and_is_an_involution() {
        // Sizes either side of the 64-bit block boundary, including ragged ones.
        for &n in &[1usize, 7, 63, 64, 65, 100, 127, 128, 129, 200] {
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
        for &n in &[7usize, 65, 100, 129] {
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
