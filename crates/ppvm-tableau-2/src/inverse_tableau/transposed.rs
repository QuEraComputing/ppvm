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
use super::bitword::{BitWord, ZERO, flip_bit, get_bit, transpose_64x64};

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
                for_each_block(n, stride, |rows, cols| {
                    let mut block = [ZERO; 64];
                    for (r, word) in block.iter_mut().enumerate() {
                        let i = rows + r;
                        if i < n {
                            let (x, z) = InverseTableau::offsets(n, stride, k, i);
                            *word = t.words[if plane == 0 { x } else { z } + cols / 64];
                        }
                    }
                    transpose_64x64(&mut block);
                    for (j, &word) in block.iter().enumerate() {
                        if cols + j < n {
                            dst[(cols + j) * stride + rows / 64] = word;
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
                for_each_block(n, stride, |rows, cols| {
                    let mut block = [ZERO; 64];
                    for (j, word) in block.iter_mut().enumerate() {
                        if cols + j < n {
                            *word = src[(cols + j) * stride + rows / 64];
                        }
                    }
                    transpose_64x64(&mut block);
                    for (r, &word) in block.iter().enumerate() {
                        let i = rows + r;
                        if i < n {
                            let (x, z) = InverseTableau::offsets(n, stride, k, i);
                            t.words[if plane == 0 { x } else { z } + cols / 64] = word;
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

    pub(super) fn negative(&self, k: usize, i: usize) -> bool {
        get_bit(&self.signs[k], i)
    }

    pub(super) fn num_qubits(&self) -> usize {
        self.n
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

/// Visit the 64 × 64 blocks of an `n × n` bit matrix by their first row and
/// first column.
fn for_each_block(n: usize, stride: usize, mut f: impl FnMut(usize, usize)) {
    for rows in (0..stride).map(|w| 64 * w).filter(|&r| r < n) {
        for cols in (0..stride).map(|w| 64 * w).filter(|&c| c < n) {
            f(rows, cols);
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
        for (i, n) in SIZES.into_iter().enumerate() {
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
