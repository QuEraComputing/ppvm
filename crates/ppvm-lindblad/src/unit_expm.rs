// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Coefficient-free Hamiltonian connections referencing the original model terms.
//!
//! A unit plan exists only when every dissipator is a Hermitian Pauli jump, so all
//! off-diagonal generator entries come from the Hamiltonian, with any noise on the
//! diagonal. In the orthonormal Pauli basis a Hamiltonian commutator is real
//! antisymmetric, term by term: row `i` reaches parent `j` through term `k` exactly when
//! row `j` reaches `i` through `k` with the opposite sign. Projection onto the working
//! basis keeps this, so the cache stores each connected pair once and the product
//! applies it in both directions.

use crate::compact_index::CompactBasisIndex;
use crate::pauli_plan::{PauliActionPlan, UnitAction};
use crate::word::Word;
use quspin_types::{LinearOperator, QuSpinError};
use rayon::prelude::*;

/// Signed Hamiltonian term code stored with each pair: 16 bits when the model has at
/// most `i16::MAX` terms, otherwise 32 bits.
trait Code: Copy + Default + Send + Sync {
    fn pack(unit: i32) -> Self;
    fn unit(self) -> i32;
}

impl Code for i16 {
    fn pack(unit: i32) -> Self {
        unit as i16
    }
    fn unit(self) -> i32 {
        self as i32
    }
}

impl Code for i32 {
    fn pack(unit: i32) -> Self {
        unit
    }
    fn unit(self) -> i32 {
        self
    }
}

/// One stored pair: block-local row and parent positions and the signed term code
/// for the row's direction. The parent's block is given by the tile holding the entry.
#[derive(Clone, Copy, Default)]
struct Entry<U> {
    row: u16,
    parent: u16,
    unit: U,
}

const _: () = assert!(std::mem::size_of::<Entry<i16>>() == 6);
const _: () = assert!(std::mem::size_of::<Entry<i32>>() == 8);

/// Largest block-size exponent that keeps block-local positions in 16 bits.
const MAX_SHIFT: u32 = 16;

/// Rows of one block and the pairs it owns, grouped by the other row's block.
struct Block<U> {
    /// `entries[tiles[b]..tiles[b + 1]]` pair rows of this block with rows of block `b`.
    tiles: Vec<u32>,
    entries: Vec<Entry<U>>,
    diagonal: Vec<f64>,
}

impl<U> Default for Block<U> {
    fn default() -> Self {
        Self {
            tiles: Vec::new(),
            entries: Vec::new(),
            diagonal: Vec::new(),
        }
    }
}

impl<U> Block<U> {
    fn allocation_bytes(&self) -> usize {
        self.tiles.capacity() * 4
            + self.entries.capacity() * std::mem::size_of::<Entry<U>>()
            + self.diagonal.capacity() * 8
    }
}

enum Blocks {
    Narrow(Vec<Block<i16>>),
    Wide(Vec<Block<i32>>),
}

pub(crate) struct UnitCache {
    blocks: Blocks,
    shift: u32,
    dimension: usize,
}

/// The block storing pairs between blocks `a` and `b`. Alternating the owner by parity
/// gives every block about half of its pairs, which balances both product passes.
fn owner(a: usize, b: usize) -> usize {
    if (a + b).is_multiple_of(2) {
        a.max(b)
    } else {
        a.min(b)
    }
}

/// Running reduction for the exponential's shift and norm bound without per-row storage.
///
/// [`crate::mf_expm::real_bounds`] takes `mu` as the mean diagonal and the bound as the
/// largest `raw - |d| + |d - mu|` over rows. Since `|d - mu| = max(d - mu, mu - d)`, that
/// maximum equals `max(max(raw - |d| + d) - mu, max(raw - |d| - d) + mu)`, which needs only
/// two running maxima and the diagonal sum.
#[derive(Clone, Copy)]
struct Bounds {
    diagonal_sum: f64,
    above: f64,
    below: f64,
    raw: f64,
}

impl Default for Bounds {
    fn default() -> Self {
        Self {
            diagonal_sum: 0.0,
            above: f64::NEG_INFINITY,
            below: f64::NEG_INFINITY,
            raw: 0.0,
        }
    }
}

impl Bounds {
    fn add(&mut self, raw: f64, diagonal: f64) {
        let off = raw - diagonal.abs();
        self.diagonal_sum += diagonal;
        self.above = self.above.max(off + diagonal);
        self.below = self.below.max(off - diagonal);
        self.raw = self.raw.max(raw);
    }

    fn merge(self, other: Self) -> Self {
        Self {
            diagonal_sum: self.diagonal_sum + other.diagonal_sum,
            above: self.above.max(other.above),
            below: self.below.max(other.below),
            raw: self.raw.max(other.raw),
        }
    }

    fn finish(self, dimension: usize, noisy: bool) -> (f64, f64) {
        if !noisy {
            return (0.0, self.raw);
        }
        let mu = self.diagonal_sum / dimension as f64;
        (mu, (self.above - mu).max(self.below + mu).max(0.0))
    }
}

impl UnitCache {
    pub(crate) fn allocation_bytes(&self) -> usize {
        match &self.blocks {
            Blocks::Narrow(b) => {
                b.capacity() * std::mem::size_of::<Block<i16>>()
                    + b.iter().map(Block::allocation_bytes).sum::<usize>()
            }
            Blocks::Wide(b) => {
                b.capacity() * std::mem::size_of::<Block<i32>>()
                    + b.iter().map(Block::allocation_bytes).sum::<usize>()
            }
        }
    }

    #[cfg(test)]
    fn block_count(&self) -> usize {
        match &self.blocks {
            Blocks::Narrow(b) => b.len(),
            Blocks::Wide(b) => b.len(),
        }
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
            .trailing_zeros()
            .min(MAX_SHIFT);
        let index = CompactBasisIndex::new(basis);
        let previous = previous.map(|cache| cache.blocks);
        let (blocks, bounds) = if plan.term_count() <= i16::MAX as usize {
            let blocks = match previous {
                Some(Blocks::Narrow(blocks)) => blocks,
                _ => Vec::new(),
            };
            let (blocks, bounds) = build_blocks(plan, basis, &index, shift, blocks);
            (Blocks::Narrow(blocks), bounds)
        } else {
            let blocks = match previous {
                Some(Blocks::Wide(blocks)) => blocks,
                _ => Vec::new(),
            };
            let (blocks, bounds) = build_blocks(plan, basis, &index, shift, blocks);
            (Blocks::Wide(blocks), bounds)
        };
        (
            Self {
                blocks,
                shift,
                dimension,
            },
            bounds,
        )
    }
}

fn build_blocks<const C: usize, U: Code>(
    plan: &PauliActionPlan<C>,
    basis: &[Word<C>],
    index: &CompactBasisIndex<'_, C>,
    shift: u32,
    mut blocks: Vec<Block<U>>,
) -> (Vec<Block<U>>, (f64, f64)) {
    let dimension = basis.len();
    let rows = 1 << shift;
    let count = dimension.div_ceil(rows);
    let noisy = plan.has_diagonal_noise();
    blocks.resize_with(count, Block::default);
    let metadata: Vec<_> = blocks
        .par_iter_mut()
        .enumerate()
        .map_init(
            || {
                (
                    Vec::<u64>::new(),
                    Vec::<UnitAction<C>>::new(),
                    Vec::<(u32, Entry<U>)>::new(),
                    Vec::<u32>::new(),
                )
            },
            |(active, actions, staged, cursor), (block, storage)| {
                let first = block * rows;
                let length = rows.min(dimension - first);
                staged.clear();
                storage.diagonal.clear();
                if noisy {
                    storage.diagonal.resize(length, 0.0);
                }
                let mut bounds = Bounds::default();
                for (row, word) in basis[first..first + length].iter().enumerate() {
                    plan.fill_units(word, active, actions);
                    // Hash order keeps successive index probes close together.
                    actions.sort_unstable_by(|a, b| {
                        a.hash.cmp(&b.hash).then_with(|| a.word.cmp(&b.word))
                    });
                    let mut raw = 0.0;
                    let mut diagonal = 0.0;
                    for action in actions.iter() {
                        raw += action.value.abs();
                        if action.unit == 0 {
                            diagonal = action.value;
                            continue;
                        }
                        let Some(parent) = index.get_hashed(&action.word, action.hash) else {
                            continue;
                        };
                        let parent = parent as usize;
                        let other = parent >> shift;
                        let stored = if other == block {
                            parent < first + row
                        } else {
                            owner(block, other) == block
                        };
                        if stored {
                            staged.push((
                                other as u32,
                                Entry {
                                    row: row as u16,
                                    parent: (parent - (other << shift)) as u16,
                                    unit: U::pack(action.unit),
                                },
                            ));
                        }
                    }
                    bounds.add(raw, diagonal);
                    if noisy {
                        storage.diagonal[row] = diagonal;
                    }
                }
                // Group entries by tile, keeping row order within each tile.
                storage.tiles.clear();
                storage.tiles.resize(count + 1, 0);
                for &(tile, _) in staged.iter() {
                    storage.tiles[tile as usize + 1] += 1;
                }
                for tile in 0..count {
                    storage.tiles[tile + 1] += storage.tiles[tile];
                }
                // Entry storage keeps its capacity between builds and grows only when a
                // build needs more, with 1/16 headroom, so capacity tracks the block's
                // entry count and the allocation stays stable across exponentials.
                storage.entries.clear();
                if storage.entries.capacity() < staged.len() {
                    storage.entries = Vec::new();
                    storage
                        .entries
                        .reserve_exact(staged.len() + staged.len() / 16);
                }
                storage.entries.resize(staged.len(), Entry::default());
                cursor.clear();
                cursor.extend_from_slice(&storage.tiles[..count]);
                for &(tile, entry) in staged.iter() {
                    let slot = &mut cursor[tile as usize];
                    storage.entries[*slot as usize] = entry;
                    *slot += 1;
                }
                bounds
            },
        )
        .collect();
    let bounds = metadata
        .into_iter()
        .fold(Bounds::default(), Bounds::merge)
        .finish(dimension, noisy);
    (blocks, bounds)
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
        let shift = self.cache.shift;
        match &self.cache.blocks {
            Blocks::Narrow(blocks) => product(blocks, shift, self.plan, overwrite, input, output),
            Blocks::Wide(blocks) => product(blocks, shift, self.plan, overwrite, input, output),
        }
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

fn product<const C: usize, U: Code>(
    blocks: &[Block<U>],
    shift: u32,
    plan: &PauliActionPlan<C>,
    overwrite: bool,
    input: &[f64],
    output: &mut [f64],
) {
    let rows = 1usize << shift;
    let slice = |block: usize| &input[block * rows..((block + 1) * rows).min(input.len())];
    // Rows gather from the pairs their block owns, in the stored direction.
    output
        .par_chunks_mut(rows)
        .zip(blocks.par_iter())
        .enumerate()
        .for_each(|(block, (output, storage))| {
            let own = slice(block);
            if storage.diagonal.is_empty() {
                if overwrite {
                    output.fill(0.0);
                }
            } else {
                for ((result, &d), &x) in output.iter_mut().zip(&storage.diagonal).zip(own) {
                    let value = if x != 0.0 { d * x } else { 0.0 };
                    *result = if overwrite { value } else { *result + value };
                }
            }
            for other in 0..blocks.len() {
                let x = slice(other);
                let tile = storage.tiles[other] as usize..storage.tiles[other + 1] as usize;
                for entry in &storage.entries[tile] {
                    let value = x[entry.parent as usize];
                    if value != 0.0 {
                        output[entry.row as usize] +=
                            plan.row_coefficient(entry.unit.unit()) * value;
                    }
                }
            }
        });
    // Rows receive the mirrored direction of pairs owned by other blocks (and of
    // their own block's pairs), with the opposite sign.
    output
        .par_chunks_mut(rows)
        .enumerate()
        .for_each(|(block, output)| {
            for (other, storage) in blocks.iter().enumerate() {
                let x = slice(other);
                let tile = storage.tiles[block] as usize..storage.tiles[block + 1] as usize;
                for entry in &storage.entries[tile] {
                    let value = x[entry.row as usize];
                    if value != 0.0 {
                        output[entry.parent as usize] -=
                            plan.row_coefficient(entry.unit.unit()) * value;
                    }
                }
            }
        });
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
                // Paired storage sums each row in a different order than the general
                // generator, so products agree to rounding rather than bit for bit.
                let close = |actual: &[f64], factor: f64| {
                    for (a, e) in actual.iter().zip(&expected) {
                        assert!(
                            (a - factor * e).abs() <= 1e-14 * (1.0 + e.abs()),
                            "{a} vs {e}"
                        );
                    }
                };
                let mut actual = vec![0.0; basis.len()];
                pool.install(|| op.dot(true, &input, &mut actual)).unwrap();
                close(&actual, 1.0);
                pool.install(|| op.dot(false, &input, &mut actual)).unwrap();
                close(&actual, 2.0);
            }
        }
    }

    #[test]
    fn streaming_bounds_match_per_row_bounds() {
        let rows: Vec<(f64, f64)> = (0..997)
            .map(|i| {
                let d = -((i * 37 % 101) as f64) * 0.013;
                (d.abs() + ((i * 53 % 89) as f64) * 0.07, d)
            })
            .collect();
        let expected = crate::mf_expm::real_bounds(&rows);
        let merged = rows
            .chunks(64)
            .map(|chunk| {
                let mut b = Bounds::default();
                for &(raw, d) in chunk {
                    b.add(raw, d);
                }
                b
            })
            .fold(Bounds::default(), Bounds::merge)
            .finish(rows.len(), true);
        assert!((merged.0 - expected.0).abs() <= 1e-15);
        assert!((merged.1 - expected.1).abs() <= 1e-13 * expected.1);
    }

    #[test]
    fn wide_term_codes_give_identical_products() {
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
        let input: Vec<_> = (0..basis.len()).map(|i| (i as f64).sin()).collect();
        let jumps = vec![JumpInput {
            lincomb: vec![("IZII".into(), Complex::new(1.0, 0.0))],
            rate: 0.23,
        }];
        let spec = LindbladSpec::<2>::new(4, &h, &jumps).unwrap();
        let plan = spec.pauli_plan.as_ref().unwrap();
        let (narrow, bounds) = UnitCache::build(plan, &basis, None);
        assert!(matches!(narrow.blocks, Blocks::Narrow(_)));
        let index = CompactBasisIndex::new(&basis);
        let (blocks, wide_bounds) =
            build_blocks::<2, i32>(plan, &basis, &index, narrow.shift, Vec::new());
        assert_eq!(bounds, wide_bounds);
        let wide = UnitCache {
            blocks: Blocks::Wide(blocks),
            shift: narrow.shift,
            dimension: narrow.dimension,
        };
        let mut a = vec![0.0; basis.len()];
        let mut b = vec![0.0; basis.len()];
        UnitOp {
            cache: &narrow,
            plan,
        }
        .dot(true, &input, &mut a)
        .unwrap();
        UnitOp { cache: &wide, plan }
            .dot(true, &input, &mut b)
            .unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn paired_products_cover_many_blocks_and_are_thread_count_invariant() {
        // Five qubits with a three-body term and unequal couplings; the basis spans
        // more than four blocks, so pairs of both owner parities and within-block
        // pairs are exercised.
        let h = vec![
            ("XXIII".into(), 0.31),
            ("IYYII".into(), -0.52),
            ("IIZXI".into(), 0.27),
            ("XIYIZ".into(), 0.11),
            ("IIIYY".into(), 0.9),
        ];
        let basis: Vec<_> = (0..1024usize)
            .filter(|code| code % 5 != 1)
            .map(|code| {
                word_from_codes::<2>(
                    &(0..5)
                        .map(|q| ((code >> (2 * q)) & 3) as u8)
                        .collect::<Vec<_>>(),
                )
                .unwrap()
            })
            .collect();
        let input: Vec<_> = (0..basis.len())
            .map(|i| {
                if i % 11 == 0 {
                    0.0
                } else {
                    (0.7 * i as f64).cos()
                }
            })
            .collect();
        for rate in [0.0, 0.4] {
            let jumps = vec![JumpInput {
                lincomb: vec![("IIZII".into(), Complex::new(1.0, 0.0))],
                rate,
            }];
            let spec = LindbladSpec::<2>::new(5, &h, &jumps).unwrap();
            let plan = spec.pauli_plan.as_ref().unwrap();
            let (cache, _) = UnitCache::build(plan, &basis, None);
            assert!(cache.block_count() > 4);
            let op = UnitOp {
                cache: &cache,
                plan,
            };
            let mut expected = vec![0.0; basis.len()];
            for (row, column, value) in spec.generator(&basis) {
                expected[row] += value * input[column];
            }
            let mut first: Option<Vec<f64>> = None;
            for threads in [1, 3, 6] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let mut actual = vec![0.0; basis.len()];
                pool.install(|| op.dot(true, &input, &mut actual)).unwrap();
                for (a, e) in actual.iter().zip(&expected) {
                    assert!((a - e).abs() <= 1e-13 * (1.0 + e.abs()), "{a} vs {e}");
                }
                match &first {
                    Some(reference) => assert_eq!(&actual, reference),
                    None => first = Some(actual),
                }
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
