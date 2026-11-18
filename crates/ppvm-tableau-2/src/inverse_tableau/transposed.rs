// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! A column-major copy of the tableau, for batches of changes of basis.
//!
//! Column `c` (virtual qubit `c`) of a row block is `stride` contiguous words,
//! one bit per row, so a change of basis runs the gate rules on 64 rows per
//! word (`B = u64`) instead of one bit per row. Building the copy and writing
//! it back each cost a transpose, so it only pays for a batch.

use ppvm_pauli_kernels::heisenberg;

use super::InverseTableau;
use super::bitword::{BitWord, LANES, ZERO, bit_position, flip_bit, get_bit, transpose_64x64};

pub(super) struct TransposedTableau {
    n: usize,
    stride: usize,
    /// `x[k][c * stride + w]`: X bits of column `c` for rows `64w..` of
    /// block `k` (0: the `ix` rows, 1: the `iz` rows).
    x: [Vec<BitWord>; 2],
    z: [Vec<BitWord>; 2],
    /// `signs[k]`: bit `i` is the sign of row `i` of block `k`.
    signs: [Vec<BitWord>; 2],
}

impl TransposedTableau {
    /// The column-major copy of `t`.
    pub(super) fn from_rows(t: &InverseTableau) -> Self {
        let (n, stride) = (t.n, t.stride);
        let empty = || [vec![ZERO; n * stride], vec![ZERO; n * stride]];
        let mut c = Self {
            n,
            stride,
            x: empty(),
            z: empty(),
            signs: [vec![ZERO; stride], vec![ZERO; stride]],
        };
        for k in 0..2 {
            for plane in 0..2 {
                let dst = if plane == 0 { &mut c.x[k] } else { &mut c.z[k] };
                for_each_block(stride, |rows, cols, lanes| {
                    let mut block = [[ZERO; LANES]; 64];
                    for (r, words) in block.iter_mut().enumerate() {
                        let i = rows + r;
                        if i < n {
                            let (x, z) = InverseTableau::offsets(n, stride, k, i);
                            let at = if plane == 0 { x } else { z } + cols / 64;
                            // A loop: `copy_from_slice` (or a zip) of < 4 words calls `memmove`.
                            #[allow(clippy::manual_memcpy)]
                            for l in 0..lanes {
                                words[l] = t.words[at + l];
                            }
                        }
                    }
                    transpose_64x64(&mut block);
                    for (j, words) in block.iter().enumerate() {
                        for (l, &word) in words[..lanes].iter().enumerate() {
                            let c = cols + 64 * l + j;
                            if c < n {
                                dst[c * stride + rows / 64] = word;
                            }
                        }
                    }
                });
            }
            for i in 0..n {
                if t.negative(k * n + i) {
                    flip_bit(&mut c.signs[k], i);
                }
            }
        }
        c
    }

    /// Write this copy back into the rows of `t`.
    pub(super) fn write_rows(&self, t: &mut InverseTableau) {
        let (n, stride) = (self.n, self.stride);
        for k in 0..2 {
            for plane in 0..2 {
                let src = if plane == 0 { &self.x[k] } else { &self.z[k] };
                for_each_block(stride, |rows, cols, lanes| {
                    let mut block = [[ZERO; LANES]; 64];
                    for (j, words) in block.iter_mut().enumerate() {
                        for (l, word) in words[..lanes].iter_mut().enumerate() {
                            let c = cols + 64 * l + j;
                            if c < n {
                                *word = src[c * stride + rows / 64];
                            }
                        }
                    }
                    transpose_64x64(&mut block);
                    for (r, words) in block.iter().enumerate() {
                        let i = rows + r;
                        if i < n {
                            let (x, z) = InverseTableau::offsets(n, stride, k, i);
                            let at = if plane == 0 { x } else { z } + cols / 64;
                            #[allow(clippy::manual_memcpy)]
                            for l in 0..lanes {
                                t.words[at + l] = words[l];
                            }
                        }
                    }
                });
            }
            for i in 0..n {
                if t.negative(k * n + i) != get_bit(&self.signs[k], i) {
                    t.flip_sign(k * n + i);
                }
            }
        }
    }

    /// X and Z bits of row `i` of block `k` at column `c`.
    pub(super) fn bits(&self, k: usize, i: usize, c: usize) -> [bool; 2] {
        let w = c * self.stride;
        [
            get_bit(&self.x[k][w..w + self.stride], i),
            get_bit(&self.z[k][w..w + self.stride], i),
        ]
    }

    /// Columns where row `i` of block `k` has an X bit, ascending, into `out`:
    /// word `i / 64` of each column, without touching the Z bits.
    pub(super) fn x_support(&self, k: usize, i: usize, out: &mut Vec<usize>) {
        let (w, mask) = bit_position(i);
        let words = self.x[k][w..].iter().step_by(self.stride);
        out.extend(
            words
                .enumerate()
                .filter(|&(_, &word)| word & mask != ZERO)
                .map(|(c, _)| c),
        );
    }

    pub(super) fn negative(&self, k: usize, i: usize) -> bool {
        get_bit(&self.signs[k], i)
    }

    /// [`InverseTableau::fan_out`] without CZ targets, on columns: the rules
    /// run on 64 rows per word.
    pub(super) fn fan_out(
        &mut self,
        pivot: usize,
        cx_targets: &[usize],
        last: impl Fn(BitWord, BitWord) -> ([BitWord; 2], BitWord),
    ) {
        let s = self.stride;
        for k in 0..2 {
            let (x, z) = (&mut self.x[k], &mut self.z[k]);
            for w in 0..s {
                let p = pivot * s + w;
                let (mut px, mut pz, mut flips) = (x[p], z[p], ZERO);
                for &t in cx_targets {
                    let i = t * s + w;
                    let ([nx, nz, tx, tz], f) = heisenberg::cnot([px, pz, x[i], z[i]]);
                    (x[i], z[i], px, pz) = (tx, tz, nx, nz);
                    flips ^= f;
                }
                let ([nx, nz], f) = last(px, pz);
                (x[p], z[p]) = (nx, nz);
                self.signs[k][w] ^= flips ^ f;
            }
        }
    }
}

/// Visit the 64 × 64 blocks of a `stride`-word square bit matrix in groups of up to
/// [`LANES`] side by side: first row, first column, and blocks in the group.
fn for_each_block(stride: usize, mut f: impl FnMut(usize, usize, usize)) {
    for rows in (0..stride).map(|w| 64 * w) {
        for w in (0..stride).step_by(LANES) {
            f(rows, 64 * w, LANES.min(stride - w));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::test_support::*;
    use super::*;

    /// Building the copy and writing it back restores the tableau exactly,
    /// and the copy's bits and signs are the rows' bits and signs.
    #[test]
    fn round_trip_and_bits_match() {
        // 300 qubits: one full group of `LANES` column blocks and a partial one.
        for (i, n) in SIZES.into_iter().chain([300]).enumerate() {
            let mut t = InverseTableau::new(n);
            random_circuit(n, 150, 0x7a + i as u64)
                .into_iter()
                .for_each(|g| apply(&mut t, g));
            let c = TransposedTableau::from_rows(&t);
            for q in 0..n {
                for col in 0..n {
                    assert_eq!(c.bits(0, q, col), [t.ix(q).x_bit(col), t.ix(q).z_bit(col)]);
                    assert_eq!(c.bits(1, q, col), [t.iz(q).x_bit(col), t.iz(q).z_bit(col)]);
                }
                assert_eq!(c.negative(1, q), t.iz(q).negative);
            }
            let mut back = InverseTableau::new(n);
            random_circuit(n, 20, 0x99)
                .into_iter()
                .for_each(|g| apply(&mut back, g));
            c.write_rows(&mut back);
            assert_eq!(back, t, "n = {n}");
        }
    }
}
