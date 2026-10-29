// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Coefficient-free Hamiltonian connections referencing the original model terms.

use crate::compact_index::CompactBasisIndex;
use crate::pauli_plan::{PauliActionPlan, UnitAction};
use crate::word::{Word, word_hash};
use quspin_types::{LinearOperator, QuSpinError};
use rayon::prelude::*;

#[derive(Clone, Copy)]
struct Connection {
    parent: u32,
    unit: i32,
}

const _: () = assert!(std::mem::size_of::<Connection>() == 8);

struct Range {
    offsets: Vec<u32>,
    connections: Vec<Connection>,
    diagonal: Vec<f64>,
}

pub(crate) struct UnitCache {
    ranges: Vec<Range>,
    shift: u32,
    dimension: usize,
}

impl UnitCache {
    pub(crate) fn allocation_bytes(&self) -> usize {
        self.ranges.capacity() * std::mem::size_of::<Range>()
            + self
                .ranges
                .iter()
                .map(|r| {
                    r.offsets.capacity() * 4
                        + r.connections.capacity() * 8
                        + r.diagonal.capacity() * 8
                })
                .sum::<usize>()
    }

    pub(crate) fn build<const C: usize>(
        plan: &PauliActionPlan<C>,
        basis: &[Word<C>],
        previous: Option<Self>,
    ) -> (Self, (f64, f64)) {
        let dimension = basis.len();
        let shift = dimension
            .div_ceil(256)
            .max(64)
            .next_power_of_two()
            .trailing_zeros();
        let rows = 1 << shift;
        let index = CompactBasisIndex::new(basis);
        let noisy = plan.has_diagonal_noise();
        let mut ranges = previous.map(|cache| cache.ranges).unwrap_or_default();
        ranges.resize_with(dimension.div_ceil(rows), || Range {
            offsets: Vec::new(),
            connections: Vec::new(),
            diagonal: Vec::new(),
        });
        let metadata: Vec<_> = ranges
            .par_iter_mut()
            .enumerate()
            .map_init(
                || {
                    (
                        Vec::<u64>::new(),
                        Vec::<UnitAction<C>>::new(),
                        Vec::<Connection>::new(),
                    )
                },
                |(active, actions, connections), (block, range)| {
                    let first = block * rows;
                    let length = rows.min(dimension - first);
                    range.connections.clear();
                    range.offsets.clear();
                    range.offsets.reserve_exact(length + 1);
                    range.offsets.push(0);
                    range.diagonal.clear();
                    if noisy {
                        range.diagonal.resize(length, 0.0);
                    }
                    let mut columns = if noisy {
                        Vec::with_capacity(length)
                    } else {
                        Vec::new()
                    };
                    let mut norm = 0.0_f64;
                    for (row, word) in basis[first..first + length].iter().enumerate() {
                        plan.fill_units(word, active, actions);
                        actions.sort_unstable_by(|a, b| {
                            word_hash(&a.word)
                                .cmp(&word_hash(&b.word))
                                .then_with(|| a.word.cmp(&b.word))
                        });
                        connections.clear();
                        let mut raw = 0.0;
                        let mut diagonal = 0.0;
                        for action in actions.iter() {
                            raw += action.value.abs();
                            if action.unit == 0 {
                                diagonal = action.value;
                            }
                            if let Some(parent) = index.get(&action.word) {
                                connections.push(Connection {
                                    parent,
                                    unit: action.unit,
                                });
                            }
                        }
                        norm = norm.max(raw);
                        if noisy {
                            range.diagonal[row] = diagonal;
                            columns.push((raw, diagonal));
                        }
                        // One row sums parents in the same order as the general sparse product.
                        connections.sort_unstable_by_key(|connection| connection.parent);
                        let required = range.connections.len() + connections.len();
                        if required > range.connections.capacity() {
                            range.connections.reserve_exact(
                                required.next_power_of_two() - range.connections.len(),
                            );
                        }
                        range.connections.extend_from_slice(connections);
                        range.offsets.push(
                            u32::try_from(range.connections.len())
                                .expect("unit row range exceeds u32 entries"),
                        );
                    }
                    (columns, norm)
                },
            )
            .collect();
        let mut columns = Vec::new();
        let mut norm = 0.0_f64;
        for (values, bound) in metadata {
            columns.extend(values);
            norm = norm.max(bound);
        }
        let bounds = if noisy {
            crate::mf_expm::real_bounds(&columns)
        } else {
            (0.0, norm)
        };
        (
            Self {
                ranges,
                shift,
                dimension,
            },
            bounds,
        )
    }
}

pub(crate) struct UnitOp<'a, const C: usize> {
    pub(crate) cache: &'a UnitCache,
    pub(crate) plan: &'a PauliActionPlan<C>,
}

impl<const C: usize> LinearOperator<f64> for UnitOp<'_, C> {
    fn dim(&self) -> usize {
        self.cache.dimension
    }
    fn parallel_hint(&self) -> bool {
        false
    }
    fn dot(&self, overwrite: bool, input: &[f64], output: &mut [f64]) -> Result<(), QuSpinError> {
        output
            .par_chunks_mut(1 << self.cache.shift)
            .zip(self.cache.ranges.par_iter())
            .for_each(|(output, range)| {
                for (row, result) in output.iter_mut().enumerate() {
                    let mut sum = 0.0;
                    for connection in &range.connections
                        [range.offsets[row] as usize..range.offsets[row + 1] as usize]
                    {
                        let value = input[connection.parent as usize];
                        if value != 0.0 {
                            let coefficient = if connection.unit == 0 {
                                range.diagonal[row]
                            } else {
                                self.plan.row_coefficient(connection.unit)
                            };
                            sum += coefficient * value;
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
    fn trace(&self) -> f64 {
        unreachable!("trace supplied to exponential")
    }
    fn onenorm(&self, _: f64) -> f64 {
        unreachable!("norm supplied to exponential")
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
    use crate::{JumpInput, LindbladSpec, PcStepConfig, word_from_codes};
    use num::Complex;

    #[test]
    fn unit_products_match_general_generator_with_unequal_couplings_and_noise() {
        let h = vec![
            ("XXII".into(), 0.3),
            ("IYYI".into(), -0.71),
            ("ZIIX".into(), 0.17),
        ];
        let basis: Vec<_> = (0..256usize)
            .filter(|code| code % 3 != 0)
            .map(|code| {
                word_from_codes::<2>(
                    &(0..4)
                        .map(|q| ((code >> (2 * q)) & 3) as u8)
                        .collect::<Vec<_>>(),
                )
                .unwrap()
            })
            .collect();
        let input: Vec<_> = (0..basis.len())
            .map(|i| if i % 7 == 0 { 0.0 } else { (i as f64).sin() })
            .collect();
        for rate in [0.0, 0.23] {
            let jumps = vec![JumpInput {
                lincomb: vec![("IZII".into(), Complex::new(1.0, 0.0))],
                rate,
            }];
            let spec = LindbladSpec::<2>::new(4, &h, &jumps).unwrap();
            let plan = spec.pauli_plan.as_ref().unwrap();
            let (cache, _) = UnitCache::build(plan, &basis, None);
            let op = UnitOp {
                cache: &cache,
                plan,
            };
            let mut expected = vec![0.0; basis.len()];
            for (row, column, value) in spec.generator(&basis) {
                if input[column] != 0.0 {
                    expected[row] += value * input[column];
                }
            }
            for threads in [1, 2, 6] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let mut actual = vec![0.0; basis.len()];
                pool.install(|| op.dot(true, &input, &mut actual)).unwrap();
                assert_eq!(actual, expected);
                pool.install(|| op.dot(false, &input, &mut actual)).unwrap();
                assert_eq!(
                    actual,
                    expected
                        .iter()
                        .map(|value| value + value)
                        .collect::<Vec<_>>()
                );
            }
        }
    }

    #[test]
    fn adding_a_general_dissipator_replaces_a_retained_unit_cache() {
        let h = vec![("XXI".into(), 0.3), ("IYY".into(), -0.17)];
        let mut reference = LindbladSpec::<2>::new(3, &h, &[]).unwrap();
        let mut optimized = LindbladSpec::<2>::new(3, &h, &[]).unwrap();
        let seed = word_from_codes::<2>(&[2, 0, 0]).unwrap();
        let (mut a, mut b) = (vec![seed], vec![seed]);
        let (mut ca, mut cb) = (vec![1.0], vec![1.0]);
        let cfg = PcStepConfig {
            max_basis: 8,
            admit_basis: Some(24),
            num_threads: Some(2),
            ..Default::default()
        };
        for step in 0..8 {
            if step == 3 {
                let operators = vec![vec![
                    ("XII".into(), Complex::new(0.5, 0.0)),
                    ("YII".into(), Complex::new(0.0, 0.5)),
                ]];
                let k = vec![vec![Complex::new(0.2, 0.0)]];
                reference.add_kossakowski(&operators, &k).unwrap();
                optimized.add_kossakowski(&operators, &k).unwrap();
            }
            reference
                .pc_step_reference(&mut a, &mut ca, 0.1, &[seed], &cfg)
                .unwrap();
            optimized
                .pc_step(&mut b, &mut cb, 0.1, &[seed], &cfg)
                .unwrap();
            assert_eq!(a, b);
            assert_eq!(ca, cb);
        }
    }
}
