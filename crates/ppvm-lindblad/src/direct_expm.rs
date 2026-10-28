// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Row-gather action for a Hamiltonian plus diagonal Pauli dissipators.

use crate::compact_index::CompactBasisIndex;
use crate::{LindbladSpec, Word};
use fxhash::FxHashMap;
use num::Complex;
use quspin_types::{LinearOperator, QuSpinError};
use rayon::prelude::*;

pub(crate) struct DirectOp<'a, const C: usize> {
    pub(crate) spec: &'a LindbladSpec<C>,
    pub(crate) basis: &'a [Word<C>],
    pub(crate) index: CompactBasisIndex<'a, C>,
}

pub(crate) struct ActionScratch<const C: usize> {
    support: Vec<u32>,
    candidates: Vec<u32>,
    local: FxHashMap<Word<C>, Complex<f64>>,
    active: Vec<u64>,
    pub(crate) terms: Vec<(Word<C>, f64)>,
    pub(crate) row: Vec<(u32, f64)>,
}

impl<const C: usize> ActionScratch<C> {
    pub(crate) fn new() -> Self {
        Self {
            support: Vec::new(),
            candidates: Vec::new(),
            local: FxHashMap::default(),
            active: Vec::new(),
            terms: Vec::new(),
            row: Vec::new(),
        }
    }
    pub(crate) fn fill(&mut self, spec: &LindbladSpec<C>, word: &Word<C>, canonical: bool) {
        if let Some(plan) = &spec.pauli_plan {
            plan.fill(word, &mut self.active, &mut self.terms);
            if canonical {
                self.terms.sort_unstable_by(|(a, _), (b, _)| {
                    crate::word::word_hash(a)
                        .cmp(&crate::word::word_hash(b))
                        .then_with(|| a.cmp(b))
                });
            }
            return;
        }

        self.terms = spec.compute_action_terms(
            word,
            &mut self.support,
            &mut self.candidates,
            &mut self.local,
        );
    }
}

impl<const C: usize> DirectOp<'_, C> {
    pub(crate) fn count_row(
        &self,
        position: usize,
        scratch: &mut ActionScratch<C>,
        canonical: bool,
        lookup: bool,
    ) -> ((f64, f64), usize) {
        let word = &self.basis[position];
        scratch.fill(self.spec, word, canonical);
        let mut raw = 0.0;
        let mut diagonal = 0.0;
        let mut count = 0;
        for &(parent, value) in &scratch.terms {
            raw += value.abs();
            if parent == *word {
                diagonal = value;
            }
            count += usize::from(!lookup || self.index.get(&parent).is_some());
        }
        ((raw, diagonal), count)
    }

    pub(crate) fn gather_row(
        &self,
        position: usize,
        scratch: &mut ActionScratch<C>,
        canonical: bool,
    ) -> (f64, f64) {
        let word = &self.basis[position];
        scratch.fill(self.spec, word, canonical);
        scratch.row.clear();
        let mut raw = 0.0;
        let mut diagonal = 0.0;
        for &(parent, value) in &scratch.terms {
            raw += value.abs();
            if parent == *word {
                diagonal = value;
            }
            if let Some(column) = self.index.get(&parent) {
                scratch
                    .row
                    .push((column, if parent == *word { value } else { -value }));
            }
        }
        scratch.row.sort_unstable_by_key(|entry| entry.0);
        (raw, diagonal)
    }

    pub(crate) fn per_col(&self) -> Vec<(f64, f64)> {
        self.basis
            .par_iter()
            .map_init(ActionScratch::new, |scratch, word| {
                scratch.fill(self.spec, word, true);
                let mut raw = 0.0;
                let mut diag = 0.0;
                for &(output, value) in &scratch.terms {
                    raw += value.abs();
                    if output == *word {
                        diag = value;
                    }
                }
                (raw, diag)
            })
            .collect()
    }
}

impl<const C: usize> LinearOperator<f64> for DirectOp<'_, C> {
    fn dim(&self) -> usize {
        self.basis.len()
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
        output.par_iter_mut().enumerate().for_each_init(
            ActionScratch::new,
            |scratch, (position, result)| {
                let word = &self.basis[position];
                scratch.fill(self.spec, word, false);
                scratch.row.clear();
                for &(parent, value) in &scratch.terms {
                    if let Some(column) = self.index.get(&parent) {
                        // Hamiltonian rows are negative columns; Pauli dissipators are diagonal.
                        scratch
                            .row
                            .push((column, if parent == *word { value } else { -value }));
                    }
                }
                // Match the cached generator's accumulation order across thread counts.
                scratch.row.sort_unstable_by_key(|entry| entry.0);
                let mut sum = 0.0;
                for &(column, value) in &scratch.row {
                    if input[column as usize] != 0.0 {
                        sum += value * input[column as usize];
                    }
                }
                if overwrite {
                    *result = sum;
                } else {
                    *result += sum;
                }
            },
        );
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
    use crate::{ExpmBackend, JumpInput, PcStepConfig, word_from_codes};

    #[test]
    fn row_gather_matches_projected_signed_generator_across_threads() {
        let h = vec![
            ("XXII".into(), 0.3),
            ("IYYI".into(), -0.7),
            ("ZIIX".into(), 0.2),
            ("XXII".into(), -0.1),
        ];
        let jumps = vec![JumpInput {
            lincomb: vec![("IZII".into(), Complex::new(1.0, 0.0))],
            rate: 0.4,
        }];
        let spec = LindbladSpec::<2>::new(4, &h, &jumps).unwrap();
        let basis: Vec<_> = (0..256usize)
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
        let input: Vec<f64> = (0..basis.len())
            .map(|i| if i % 7 == 0 { 0.0 } else { (i as f64).sin() })
            .collect();
        let mut expected = vec![0.0; basis.len()];
        for (row, col, value) in spec.generator(&basis) {
            if input[col] != 0.0 {
                expected[row] += value * input[col];
            }
        }
        let op = DirectOp {
            spec: &spec,
            basis: &basis,
            index: CompactBasisIndex::new(&basis),
        };
        for threads in [1, 2, 6] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let mut actual = vec![0.0; basis.len()];
            pool.install(|| op.dot(true, &input, &mut actual)).unwrap();
            assert_eq!(actual, expected, "thread count {threads}");
            pool.install(|| op.dot(false, &input, &mut actual)).unwrap();
            assert_eq!(actual, expected.iter().map(|v| v + v).collect::<Vec<_>>());
        }
    }

    #[test]
    fn direct_exponential_preserves_adaptive_steps() {
        let h = vec![
            ("XXI".into(), 0.3),
            ("IYY".into(), 0.7),
            ("ZII".into(), -0.2),
        ];
        let cached = LindbladSpec::<2>::new(3, &h, &[]).unwrap();
        let direct = LindbladSpec::<2>::new(3, &h, &[])
            .unwrap()
            .with_expm_backend(ExpmBackend::MatrixFree)
            .with_leakage_backend(crate::LeakageBackend::Indexed);
        let seed = word_from_codes::<2>(&[3, 0, 0]).unwrap();
        let mut a = vec![seed];
        let mut b = a.clone();
        let mut ca = vec![1.0];
        let mut cb = ca.clone();
        let mut reused_basis = a.clone();
        let mut reused_coeffs = ca.clone();
        let mut workspace = crate::PcStepWorkspace::default();
        let tail = LindbladSpec::<2>::new(3, &h, &[])
            .unwrap()
            .with_expm_backend(ExpmBackend::MatrixFree)
            .with_leakage_backend(crate::LeakageBackend::Tail);
        let mut tail_basis = a.clone();
        let mut tail_coeffs = ca.clone();
        let mut tail_workspace = crate::PcStepWorkspace::default();
        let cfg = PcStepConfig {
            max_basis: 8,
            admit_basis: Some(24),
            num_threads: Some(2),
            ..Default::default()
        };
        for _ in 0..10 {
            cached.pc_step(&mut a, &mut ca, 0.1, &[seed], &cfg).unwrap();
            direct.pc_step(&mut b, &mut cb, 0.1, &[seed], &cfg).unwrap();
            direct
                .pc_step_with_workspace(
                    &mut reused_basis,
                    &mut reused_coeffs,
                    0.1,
                    &[seed],
                    &cfg,
                    &mut workspace,
                )
                .unwrap();
            tail.pc_step_with_workspace(
                &mut tail_basis,
                &mut tail_coeffs,
                0.1,
                &[seed],
                &cfg,
                &mut tail_workspace,
            )
            .unwrap();
            assert_eq!(a, tail_basis);
            assert_eq!(ca, tail_coeffs);
            assert_eq!(a, reused_basis);
            assert_eq!(ca, reused_coeffs);
            assert_eq!(a, b);
            assert_eq!(ca, cb);
        }
    }

    #[test]
    fn packed_workspace_preserves_predictor_corrector_and_noise() {
        let h = vec![
            ("XXI".into(), 0.3),
            ("IYY".into(), -0.3),
            ("ZII".into(), 0.3),
        ];
        let jumps = vec![JumpInput {
            lincomb: vec![("IZI".into(), Complex::new(1.0, 0.0))],
            rate: 0.2,
        }];
        for backend in [
            ExpmBackend::Packed,
            ExpmBackend::PackedUpper,
            ExpmBackend::PackedBatched,
        ] {
            let cached = LindbladSpec::<2>::new(3, &h, &jumps).unwrap();
            let optimized = LindbladSpec::<2>::new(3, &h, &jumps)
                .unwrap()
                .with_expm_backend(backend)
                .with_leakage_backend(crate::LeakageBackend::Tail);
            let seed = word_from_codes::<2>(&[3, 0, 0]).unwrap();
            let mut a = vec![seed];
            let mut b = a.clone();
            let mut ca = vec![1.0];
            let mut cb = ca.clone();
            let mut workspace = crate::PcStepWorkspace::default();
            for step in 0..12 {
                let cfg = PcStepConfig {
                    max_basis: if step == 0 { usize::MAX } else { 8 + step % 3 },
                    admit_basis: if step == 0 { None } else { Some(24 + step % 4) },
                    num_threads: Some(2),
                    ..Default::default()
                };
                cached.pc_step(&mut a, &mut ca, 0.1, &[seed], &cfg).unwrap();
                optimized
                    .pc_step_with_workspace(&mut b, &mut cb, 0.1, &[seed], &cfg, &mut workspace)
                    .unwrap();
                assert_eq!(a, b);
                assert_eq!(ca, cb);
            }
        }
    }

    #[test]
    fn compact_lookup_checks_words_and_wide_registers() {
        let basis: Vec<_> = (0..192)
            .map(|site| {
                let mut codes = vec![0; 128];
                codes[site / 3] = (site % 3 + 1) as u8;
                word_from_codes::<2>(&codes).unwrap()
            })
            .collect();
        let index = CompactBasisIndex::new(&basis);
        for (i, word) in basis.iter().enumerate() {
            assert_eq!(index.get(word), Some(i as u32));
        }
        assert_eq!(index.get(&word_from_codes::<2>(&[2; 128]).unwrap()), None);
        let wide = vec![word_from_codes::<8>(&vec![1; 512]).unwrap()];
        assert_eq!(CompactBasisIndex::new(&wide).get(&wide[0]), Some(0));
    }
}
