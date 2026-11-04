// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Row-gather action for a Hamiltonian plus diagonal Pauli dissipators.

use crate::compact_index::CompactBasisIndex;
use crate::{LindbladSpec, Word};
use fxhash::FxHashMap;
use num::Complex;

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
    /// Output words with their [`crate::word::word_hash`] and generator values.
    pub(crate) terms: Vec<(u64, Word<C>, f64)>,
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
            plan.fill_keyed(word, &mut self.active, &mut self.terms);
            if canonical {
                self.terms
                    .sort_unstable_by(|(ha, a, _), (hb, b, _)| ha.cmp(hb).then_with(|| a.cmp(b)));
            }
            return;
        }

        let terms = spec.compute_action_terms(
            word,
            &mut self.support,
            &mut self.candidates,
            &mut self.local,
        );
        self.terms.clear();
        self.terms.extend(
            terms
                .into_iter()
                .map(|(word, value)| (crate::word::word_hash(&word), word, value)),
        );
    }
}

impl<const C: usize> DirectOp<'_, C> {
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
        for &(hash, parent, value) in &scratch.terms {
            raw += value.abs();
            if parent == *word {
                diagonal = value;
            }
            if let Some(column) = self.index.get_hashed(&parent, hash) {
                scratch
                    .row
                    .push((column, if parent == *word { value } else { -value }));
            }
        }
        scratch.row.sort_unstable_by_key(|entry| entry.0);
        (raw, diagonal)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{JumpInput, PcStepConfig, word_from_codes};
    use quspin_types::LinearOperator;

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
        let (mat, _) = crate::mf_expm::build_row_cache(&spec, &basis, None);
        let op = crate::mf_expm::CsrOp { mat: &mat };
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
    fn default_exponential_preserves_nonuniform_adaptive_steps() {
        let h = vec![
            ("XXI".into(), 0.3),
            ("IYY".into(), 0.7),
            ("ZII".into(), -0.2),
        ];
        let cached = LindbladSpec::<2>::new(3, &h, &[]).unwrap();
        let direct = LindbladSpec::<2>::new(3, &h, &[]).unwrap();
        let seed = word_from_codes::<2>(&[3, 0, 0]).unwrap();
        let mut a = vec![seed];
        let mut b = a.clone();
        let mut ca = vec![1.0];
        let mut cb = ca.clone();
        let mut reused_basis = a.clone();
        let mut reused_coeffs = ca.clone();
        let mut workspace = crate::PcStepWorkspace::default();
        let tail = LindbladSpec::<2>::new(3, &h, &[]).unwrap();
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
            cached
                .pc_step_reference(&mut a, &mut ca, 0.1, &[seed], &cfg)
                .unwrap();
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
    fn default_workspace_preserves_predictor_corrector_and_noise() {
        let h = vec![
            ("XXI".into(), 0.3),
            ("IYY".into(), -0.3),
            ("ZII".into(), 0.3),
        ];
        let jumps = vec![JumpInput {
            lincomb: vec![("IZI".into(), Complex::new(1.0, 0.0))],
            rate: 0.2,
        }];

        let cached = LindbladSpec::<2>::new(3, &h, &jumps).unwrap();
        let optimized = LindbladSpec::<2>::new(3, &h, &jumps).unwrap();
        let seed = word_from_codes::<2>(&[3, 0, 0]).unwrap();
        let mut a = vec![seed];
        let mut b = a.clone();
        let mut ca = vec![1.0];
        let mut cb = ca.clone();
        b.reserve_exact(64 - b.len());
        cb.reserve_exact(64 - cb.len());
        let basis_storage = b.as_ptr();
        let coefficient_storage = cb.as_ptr();
        for step in 0..12 {
            let cfg = PcStepConfig {
                max_basis: if step == 0 { usize::MAX } else { 8 + step % 3 },
                admit_basis: if step == 0 { None } else { Some(24 + step % 4) },
                num_threads: Some(2),
                ..Default::default()
            };
            cached
                .pc_step_reference(&mut a, &mut ca, 0.1, &[seed], &cfg)
                .unwrap();
            optimized
                .pc_step(&mut b, &mut cb, 0.1, &[seed], &cfg)
                .unwrap();
            assert_eq!(b.as_ptr(), basis_storage);
            assert_eq!(cb.as_ptr(), coefficient_storage);
            assert_eq!(a, b);
            assert_eq!(ca, cb);
        }
    }

    #[test]
    fn default_workspace_matches_reference_beyond_one_hundred_twenty_eight_qubits() {
        let n = 129;
        let term = |symbol: char, sites: &[usize]| {
            let mut symbols = vec!['I'; n];
            for &site in sites {
                symbols[site] = symbol;
            }
            symbols.into_iter().collect::<String>()
        };
        let h = vec![
            (term('X', &[0, 128]), 0.25),
            (term('Y', &[64, 128]), -0.25),
            (term('Z', &[0]), 0.25),
        ];
        let jumps = vec![JumpInput {
            lincomb: vec![(term('Z', &[64]), Complex::new(1.0, 0.0))],
            rate: 0.2,
        }];
        let reference = LindbladSpec::<4>::new(n, &h, &jumps).unwrap();
        let optimized = LindbladSpec::<4>::new(n, &h, &jumps).unwrap();
        let mut codes = vec![0; n];
        codes[128] = 2;
        let seed = word_from_codes::<4>(&codes).unwrap();
        let mut a = vec![seed];
        let mut b = a.clone();
        let mut ca = vec![1.0];
        let mut cb = ca.clone();
        let mut workspace = crate::PcStepWorkspace::default();
        let cfg = PcStepConfig {
            max_basis: 16,
            admit_basis: Some(48),
            num_threads: Some(2),
            ..Default::default()
        };
        for _ in 0..8 {
            reference
                .pc_step_reference(&mut a, &mut ca, 0.1, &[seed], &cfg)
                .unwrap();
            optimized
                .pc_step_with_workspace(&mut b, &mut cb, 0.1, &[seed], &cfg, &mut workspace)
                .unwrap();
            let mut expected: Vec<_> = a.iter().copied().zip(ca.iter().copied()).collect();
            let mut actual: Vec<_> = b.iter().copied().zip(cb.iter().copied()).collect();
            expected.sort_unstable_by_key(|entry| entry.0);
            actual.sort_unstable_by_key(|entry| entry.0);
            assert_eq!(expected.len(), actual.len());
            for ((word_a, value_a), (word_b, value_b)) in expected.iter().zip(&actual) {
                assert_eq!(word_a, word_b);
                assert!((value_a - value_b).abs() < 2e-13);
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
