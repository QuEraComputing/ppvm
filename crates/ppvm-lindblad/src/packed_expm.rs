// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Uniform-coupling generator rows stored as signed parent positions.

use crate::direct_expm::{ActionScratch, DirectOp};
use quspin_types::{LinearOperator, QuSpinError};
use rayon::prelude::*;

struct Range {
    first: usize,
    offsets: Vec<u32>,
    entries: Vec<i32>,
    diagonal: Vec<f64>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compact_index::CompactBasisIndex;
    use crate::{ExpmBackend, JumpInput, LindbladSpec, word_from_codes};
    use num::Complex;

    #[test]
    fn signed_rows_match_projection_with_diagonal_noise_and_negative_couplings() {
        let h = vec![
            ("XXII".into(), 0.25),
            ("IYYI".into(), -0.25),
            ("ZIIX".into(), 0.25),
        ];
        let jumps = vec![JumpInput {
            lincomb: vec![("IZII".into(), Complex::new(1.0, 0.0))],
            rate: 0.3,
        }];
        let spec = LindbladSpec::<2>::new(4, &h, &jumps)
            .unwrap()
            .with_expm_backend(ExpmBackend::Packed);
        let basis: Vec<_> = (0..256)
            .filter(|i| i % 3 != 0)
            .map(|i| {
                word_from_codes::<2>(
                    &(0..4)
                        .map(|q| ((i >> (2 * q)) & 3) as u8)
                        .collect::<Vec<_>>(),
                )
                .unwrap()
            })
            .collect();
        let direct = DirectOp {
            spec: &spec,
            basis: &basis,
            index: CompactBasisIndex::new(&basis),
        };
        let (packed, metadata) = PackedOp::build(&direct, 0.5);
        assert_eq!(metadata, direct.per_col());
        let input: Vec<_> = (0..basis.len())
            .map(|i| if i % 7 == 0 { 0.0 } else { (i as f64).cos() })
            .collect();
        let mut expected = vec![0.0; basis.len()];
        for (row, col, value) in spec.generator(&basis) {
            if input[col] != 0.0 {
                expected[row] += value * input[col];
            }
        }
        for threads in [1, 2, 6] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let mut actual = vec![0.0; basis.len()];
            pool.install(|| packed.dot(true, &input, &mut actual))
                .unwrap();
            assert_eq!(actual, expected);
            pool.install(|| packed.dot(false, &input, &mut actual))
                .unwrap();
            assert_eq!(actual, expected.iter().map(|v| v + v).collect::<Vec<_>>());
        }
    }
}

pub(crate) struct PackedOp {
    ranges: Vec<Range>,
    shift: u32,
    dimension: usize,
    scale: f64,
}

impl PackedOp {
    pub(crate) fn build<const C: usize>(
        direct: &DirectOp<'_, C>,
        scale: f64,
    ) -> (Self, Vec<(f64, f64)>) {
        let dimension = direct.basis.len();
        assert!(dimension < i32::MAX as usize);
        let shift = dimension
            .div_ceil(256)
            .max(64)
            .next_power_of_two()
            .trailing_zeros();
        let (ranges, per_range): (Vec<_>, Vec<_>) = direct
            .basis
            .par_chunks(1 << shift)
            .enumerate()
            .map(|(block, words)| {
                let first = block << shift;
                let mut scratch = ActionScratch::new();
                let mut offsets = Vec::with_capacity(words.len() + 1);
                offsets.push(0u32);
                let mut per_col = Vec::with_capacity(words.len());
                for row in first..first + words.len() {
                    per_col.push(direct.gather_row(row, &mut scratch, true));
                    let count = u32::try_from(scratch.row.len()).expect("row exceeds u32 entries");
                    offsets.push(
                        offsets
                            .last()
                            .unwrap()
                            .checked_add(count)
                            .expect("range exceeds u32 entries"),
                    );
                }
                let mut entries = vec![0; *offsets.last().unwrap() as usize];
                for row in 0..words.len() {
                    direct.gather_row(first + row, &mut scratch, false);
                    let start = offsets[row] as usize;
                    for (slot, &(column, value)) in entries[start..offsets[row + 1] as usize]
                        .iter_mut()
                        .zip(&scratch.row)
                    {
                        *slot = if column as usize == first + row {
                            0
                        } else {
                            assert!(value.abs() == scale);
                            let parent = column as i32 + 1;
                            if value < 0.0 { -parent } else { parent }
                        };
                    }
                }
                let diagonal = if per_col.iter().any(|(_, value)| *value != 0.0) {
                    per_col.iter().map(|(_, value)| *value).collect()
                } else {
                    Vec::new()
                };
                (
                    Range {
                        first,
                        offsets,
                        entries,
                        diagonal,
                    },
                    per_col,
                )
            })
            .unzip();
        let per_col = per_range.into_iter().flatten().collect();
        (
            Self {
                ranges,
                shift,
                dimension,
                scale,
            },
            per_col,
        )
    }
}

impl LinearOperator<f64> for PackedOp {
    fn dim(&self) -> usize {
        self.dimension
    }
    fn parallel_hint(&self) -> bool {
        false
    }
    fn trace(&self) -> f64 {
        unreachable!("trace supplied to exponential")
    }
    fn onenorm(&self, _: f64) -> f64 {
        unreachable!("norm supplied to exponential")
    }

    fn dot(&self, overwrite: bool, input: &[f64], output: &mut [f64]) -> Result<(), QuSpinError> {
        output
            .par_chunks_mut(1 << self.shift)
            .zip(self.ranges.par_iter())
            .for_each(|(output, range)| {
                for (row, result) in output.iter_mut().enumerate() {
                    let mut sum = 0.0;
                    for &entry in
                        &range.entries[range.offsets[row] as usize..range.offsets[row + 1] as usize]
                    {
                        let (column, value) = if entry == 0 {
                            (range.first + row, range.diagonal[row])
                        } else {
                            (
                                entry.unsigned_abs() as usize - 1,
                                if entry < 0 { -self.scale } else { self.scale },
                            )
                        };
                        if input[column] != 0.0 {
                            sum += value * input[column];
                        }
                    }
                    if overwrite {
                        *result = sum;
                    } else {
                        *result += sum;
                    }
                }
            });
        Ok(())
    }

    fn dot_transpose(&self, _: bool, _: &[f64], _: &mut [f64]) -> Result<(), QuSpinError> {
        unreachable!("transpose unused by single-vector exponential")
    }
    fn dot_many(
        &self,
        _: bool,
        _: ndarray::ArrayView2<'_, f64>,
        _: ndarray::ArrayViewMut2<'_, f64>,
    ) -> Result<(), QuSpinError> {
        unreachable!("multiple vectors unused by single-vector exponential")
    }
    fn dot_chunk(&self, _: bool, _: &[f64], _: &mut [f64], _: usize) -> Result<(), QuSpinError> {
        unreachable!("chunk products unused by single-vector exponential")
    }
    fn dot_transpose_chunk(
        &self,
        _: &[f64],
        _: &[<f64 as quspin_types::ExpmComputation>::Atomic],
        _: std::ops::Range<usize>,
    ) -> Result<(), QuSpinError> {
        unreachable!("transpose chunks unused by single-vector exponential")
    }
}
