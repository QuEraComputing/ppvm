// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Uniform-coupling generator rows stored as signed parent positions.

use crate::direct_expm::{ActionScratch, DirectOp};
use quspin_types::{LinearOperator, QuSpinError};
use rayon::prelude::*;

struct Range {
    first: usize,
    offsets: Vec<u32>,
    diagonal: Vec<f64>,
    entries: Vec<i32>,
}

pub(crate) enum Norm {
    Columns(Vec<(f64, f64)>),
    Uniform(f64),
}

pub(crate) struct PackedOp {
    ranges: Vec<Range>,
    shift: u32,
    dimension: usize,
    scale: f64,
}

impl PackedOp {
    pub(crate) fn allocation_bytes(&self) -> usize {
        self.ranges.capacity() * std::mem::size_of::<Range>()
            + self
                .ranges
                .iter()
                .map(|r| {
                    r.offsets.capacity() * 4 + r.diagonal.capacity() * 8 + r.entries.capacity() * 4
                })
                .sum::<usize>()
    }

    pub(crate) fn build<const C: usize>(
        direct: &DirectOp<'_, C>,
        scale: f64,
        previous: Option<Self>,
    ) -> (Self, Norm) {
        let dimension = direct.basis.len();
        assert!(dimension < i32::MAX as usize);
        let shift = dimension
            .div_ceil(256)
            .max(64)
            .next_power_of_two()
            .trailing_zeros();
        let rows = 1 << shift;
        let noisy = direct
            .spec
            .pauli_plan
            .as_ref()
            .unwrap()
            .has_diagonal_noise();
        let mut ranges = previous.map(|p| p.ranges).unwrap_or_default();
        ranges.resize_with(dimension.div_ceil(rows), || Range {
            first: 0,
            offsets: Vec::new(),
            diagonal: Vec::new(),
            entries: Vec::new(),
        });
        let metadata: Vec<_> = ranges
            .par_iter_mut()
            .enumerate()
            .map_init(ActionScratch::new, |scratch, (block, range)| {
                let first = block * rows;
                let length = rows.min(dimension - first);
                range.first = first;
                range.entries.clear();
                range.offsets.clear();
                range.offsets.reserve_exact(length + 1);
                range.offsets.push(0);
                let mut per_col = if noisy {
                    Vec::with_capacity(length)
                } else {
                    Vec::new()
                };
                let mut bound = 0.0_f64;
                for row in first..first + length {
                    let value = direct.gather_row(row, scratch, noisy);
                    bound = bound.max(value.0);
                    if noisy {
                        per_col.push(value);
                    }
                    if range.entries.len() + scratch.row.len() > range.entries.capacity() {
                        let capacity =
                            (range.entries.len() + scratch.row.len()).next_power_of_two();
                        range.entries.reserve_exact(capacity - range.entries.len());
                    }
                    for &(column, value) in &scratch.row {
                        range.entries.push(if column as usize == row {
                            0
                        } else {
                            assert!(value.abs() == scale);
                            let parent = column as i32 + 1;
                            if value < 0.0 { -parent } else { parent }
                        });
                    }
                    range.offsets.push(
                        u32::try_from(range.entries.len()).expect("range exceeds u32 entries"),
                    );
                }
                range.diagonal.clear();
                if per_col.iter().any(|(_, v)| *v != 0.0) {
                    range.diagonal.extend(per_col.iter().map(|(_, v)| *v));
                }
                (per_col, bound)
            })
            .collect();
        let mut per_col = if noisy {
            Vec::with_capacity(dimension)
        } else {
            Vec::new()
        };
        let mut norm = 0.0_f64;
        for (data, bound) in metadata {
            per_col.extend(data);
            norm = norm.max(bound);
        }
        (
            Self {
                ranges,
                shift,
                dimension,
                scale,
            },
            if noisy {
                Norm::Columns(per_col)
            } else {
                Norm::Uniform(norm)
            },
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compact_index::CompactBasisIndex;
    use crate::{ExpmBackend, JumpInput, LindbladSpec, word_from_codes};
    use num::Complex;

    #[test]
    fn uniform_bound_matches_canonical_metadata_without_noise() {
        let h = vec![
            ("XXII".into(), 0.1),
            ("IYYI".into(), -0.1),
            ("ZIIX".into(), 0.1),
        ];
        let spec = LindbladSpec::<2>::new(4, &h, &[]).unwrap();
        let basis: Vec<_> = (0..256)
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
        let (_, Norm::Uniform(bound)) = PackedOp::build(&direct, 0.2, None) else {
            panic!("Hamiltonian has zero diagonal");
        };
        assert_eq!(
            bound,
            direct
                .per_col()
                .iter()
                .map(|(raw, _)| *raw)
                .fold(0.0_f64, f64::max)
        );
    }

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
        let (packed, Norm::Columns(metadata)) = PackedOp::build(&direct, 0.5, None) else {
            panic!("noise requires column bounds");
        };
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
