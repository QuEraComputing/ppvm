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
    use crate::{JumpInput, word_from_codes};
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
