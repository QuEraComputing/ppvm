// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Uniform-coupling generator rows stored as signed parent positions.

use crate::direct_expm::{ActionScratch, DirectOp};
use quspin_types::{LinearOperator, QuSpinError};
use rayon::prelude::*;

struct Range {
    first: usize,
    offsets: Vec<u32>,
    start: usize,
    ends: Vec<u32>,
    diagonal: Vec<f64>,
    entries: Vec<i32>,
}

pub(crate) enum Norm {
    Columns(Vec<(f64, f64)>),
    Uniform(f64),
}

pub(crate) struct PackedOp {
    ranges: Vec<Range>,
    entries: Vec<i32>,
    staging: Vec<Vec<i32>>,
    shift: u32,
    dimension: usize,
    scale: f64,
}

impl PackedOp {
    pub(crate) fn allocation_bytes(&self) -> usize {
        self.staging.iter().map(|v| v.capacity() * 4).sum::<usize>()
            + self.entries.capacity() * 4
            + self.ranges.capacity() * std::mem::size_of::<Range>()
            + self
                .ranges
                .iter()
                .map(|r| {
                    (r.offsets.capacity() + r.ends.capacity()) * 4
                        + r.diagonal.capacity() * 8
                        + r.entries.capacity() * 4
                })
                .sum::<usize>()
    }

    pub(crate) fn build<const C: usize>(
        direct: &DirectOp<'_, C>,
        scale: f64,
        previous: Option<Self>,
    ) -> (Self, Norm) {
        if direct.spec.expm_backend == crate::ExpmBackend::PackedSegmented {
            return Self::build_segmented(direct, scale, previous);
        }
        if direct.spec.expm_backend == crate::ExpmBackend::PackedBatched {
            return Self::build_batched(direct, scale, previous);
        }
        let dimension = direct.basis.len();
        assert!(dimension < i32::MAX as usize);
        let shift = dimension
            .div_ceil(256)
            .max(64)
            .next_power_of_two()
            .trailing_zeros();
        let upper = direct.spec.expm_backend == crate::ExpmBackend::PackedUpper;
        let noisy = direct
            .spec
            .pauli_plan
            .as_ref()
            .unwrap()
            .has_diagonal_noise();
        let built: Vec<_> = direct
            .basis
            .par_chunks(1 << shift)
            .enumerate()
            .map(|(block, words)| {
                let first = block << shift;
                let mut scratch = ActionScratch::new();
                let mut offsets = Vec::with_capacity(words.len() + 1);
                offsets.push(0u32);
                let mut per_col = if noisy {
                    Vec::with_capacity(words.len())
                } else {
                    Vec::new()
                };
                let mut norm = 0.0_f64;
                for row in first..first + words.len() {
                    let (metadata, count) = direct.count_row(row, &mut scratch, noisy, !upper);
                    if noisy {
                        per_col.push(metadata);
                    }
                    norm = norm.max(metadata.0);
                    let count = u32::try_from(count).expect("row exceeds u32 entries");
                    offsets.push(
                        offsets
                            .last()
                            .unwrap()
                            .checked_add(count)
                            .expect("range exceeds u32 entries"),
                    );
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
                        start: 0,
                        entries: Vec::new(),
                        ends: if upper {
                            vec![0; words.len()]
                        } else {
                            Vec::new()
                        },
                        diagonal,
                    },
                    per_col,
                    norm,
                )
            })
            .collect();
        let mut ranges = Vec::with_capacity(built.len());
        let mut per_col = if noisy {
            Vec::with_capacity(dimension)
        } else {
            Vec::new()
        };
        let mut norm = 0.0_f64;
        for (range, metadata, bound) in built {
            ranges.push(range);
            per_col.extend(metadata);
            norm = norm.max(bound);
        }
        let mut count = 0;
        for range in &mut ranges {
            range.start = count;
            count += *range.offsets.last().unwrap() as usize;
        }
        let mut entries = previous.map(|packed| packed.entries).unwrap_or_default();
        if count > entries.capacity() {
            entries.reserve_exact(count - entries.len());
        }
        entries.resize(count, 0);
        let mut remaining = entries.as_mut_slice();
        let chunks: Vec<_> = ranges
            .iter()
            .map(|range| {
                let (chunk, tail) = std::mem::take(&mut remaining)
                    .split_at_mut(*range.offsets.last().unwrap() as usize);
                remaining = tail;
                chunk
            })
            .collect();
        chunks
            .into_par_iter()
            .zip(ranges.par_iter_mut())
            .for_each_init(ActionScratch::new, |scratch, (entries, range)| {
                for row in 0..range.offsets.len() - 1 {
                    direct.gather_row(range.first + row, scratch, false);
                    let start = range.offsets[row] as usize;
                    let end = start + scratch.row.len();
                    assert!(end <= range.offsets[row + 1] as usize);
                    if upper {
                        range.ends[row] = end as u32;
                    }
                    for (slot, &(column, value)) in entries[start..end].iter_mut().zip(&scratch.row)
                    {
                        *slot = if column as usize == range.first + row {
                            0
                        } else {
                            assert!(value.abs() == scale);
                            let parent = column as i32 + 1;
                            if value < 0.0 { -parent } else { parent }
                        };
                    }
                }
            });
        (
            Self {
                ranges,
                entries,
                staging: Vec::new(),
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
    fn build_batched<const C: usize>(
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
        let (mut entries, mut staging) =
            previous.map(|p| (p.entries, p.staging)).unwrap_or_default();
        entries.clear();
        staging.resize_with(16, Vec::new);
        let mut ranges = Vec::with_capacity(dimension.div_ceil(rows));
        let mut per_col = if noisy {
            Vec::with_capacity(dimension)
        } else {
            Vec::new()
        };
        let mut norm = 0.0_f64;
        for (batch, words) in direct.basis.chunks(rows * 16).enumerate() {
            let built: Vec<_> = words
                .par_chunks(rows)
                .zip(staging.par_iter_mut())
                .enumerate()
                .map_init(ActionScratch::new, |scratch, (block, (words, stage))| {
                    let first = (batch * 16 + block) * rows;
                    stage.clear();
                    let mut offsets = Vec::with_capacity(words.len() + 1);
                    offsets.push(0u32);
                    let mut metadata = if noisy {
                        Vec::with_capacity(words.len())
                    } else {
                        Vec::new()
                    };
                    let mut bound = 0.0_f64;
                    for row in first..first + words.len() {
                        let value = direct.gather_row(row, scratch, noisy);
                        bound = bound.max(value.0);
                        if noisy {
                            metadata.push(value);
                        }
                        for &(column, value) in &scratch.row {
                            stage.push(if column as usize == row {
                                0
                            } else {
                                assert!(value.abs() == scale);
                                let parent = column as i32 + 1;
                                if value < 0.0 { -parent } else { parent }
                            });
                        }
                        offsets
                            .push(u32::try_from(stage.len()).expect("range exceeds u32 entries"));
                    }
                    let diagonal = if metadata.iter().any(|(_, v)| *v != 0.0) {
                        metadata.iter().map(|(_, v)| *v).collect()
                    } else {
                        Vec::new()
                    };
                    (
                        Range {
                            first,
                            offsets,
                            start: 0,
                            entries: Vec::new(),
                            ends: Vec::new(),
                            diagonal,
                        },
                        metadata,
                        bound,
                    )
                })
                .collect();
            let required = entries.len()
                + staging
                    .iter()
                    .take(built.len())
                    .map(Vec::len)
                    .sum::<usize>();
            if required > entries.capacity() {
                entries.reserve_exact((required - entries.len()).max(entries.capacity() / 4));
            }
            for ((mut range, metadata, bound), stage) in built.into_iter().zip(&staging) {
                range.start = entries.len();
                entries.extend_from_slice(stage);
                ranges.push(range);
                per_col.extend(metadata);
                norm = norm.max(bound);
            }
        }
        for stage in &mut staging {
            stage.clear();
        }
        (
            Self {
                ranges,
                entries,
                staging,
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
    fn build_segmented<const C: usize>(
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
            start: 0,
            offsets: Vec::new(),
            ends: Vec::new(),
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
                range.start = 0;
                range.ends.clear();
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
                entries: Vec::new(),
                staging: Vec::new(),
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
                let entries = if self.entries.is_empty() {
                    &range.entries[..]
                } else {
                    &self.entries[range.start..]
                };
                for (row, result) in output.iter_mut().enumerate() {
                    let mut sum = 0.0;
                    for &entry in &entries[range.offsets[row] as usize..if range.ends.is_empty() {
                        range.offsets[row + 1]
                    } else {
                        range.ends[row]
                    }
                        as usize]
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
        for backend in [
            ExpmBackend::Packed,
            ExpmBackend::PackedUpper,
            ExpmBackend::PackedBatched,
            ExpmBackend::PackedSegmented,
        ] {
            let spec = LindbladSpec::<2>::new(4, &h, &jumps)
                .unwrap()
                .with_expm_backend(backend);
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
}
