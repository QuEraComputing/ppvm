// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::BuildHasher;

use ppvm_traits_2::{Clifford, Pauli, Reset};
use rand::RngExt;

use super::fingerprint::fingerprint;
use super::{Branch, GeneralizedTableauMixture};
use crate::{Bitstring, GeneralizedTableau, MeasureScratch, RowStorage};

impl<A, I, H> GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    pub(crate) fn for_each_z_branch(
        &mut self,
        qubit: usize,
        mut visit: impl FnMut(&mut GeneralizedTableau<A, I, H>, Option<bool>, f64),
    ) {
        let original_len = self.entries.len();
        let mut branches: Vec<Branch<A, I, H>> = Vec::with_capacity(original_len);
        let mut scratch = MeasureScratch::new();
        let mut other_scratch = MeasureScratch::new();

        for index in 0..original_len {
            let (tab, probability) = &mut self.entries[index];
            if tab.is_lost[qubit] {
                visit(tab, None, *probability);
                continue;
            }
            let (phase, stab, destab) = tab.compute_decomposition(qubit, Pauli::Z);
            let tab_seed = self.rng.random();

            if stab == I::zero() {
                let entries = tab.coefficients.take();
                let mut other = tab.fork(Some(tab_seed));
                let overlap =
                    GeneralizedTableau::<A, I, H>::compute_overlap_case_b(&entries, phase, destab);
                let p_one = 0.5 - 0.5 * overlap;
                let likely = p_one > 0.5;
                let p_likely = if likely { p_one } else { 1.0 - p_one };
                let p_other = 1.0 - p_likely;
                let other_probability = *probability * p_other;
                if p_other > self.sum_cutoff {
                    other.project_case_b(&entries, !likely, phase, destab);
                    visit(&mut other, Some(!likely), other_probability);
                    let fp = fingerprint(&other);
                    branches.push((other, other_probability, fp));
                }
                tab.project_case_b(&entries, likely, phase, destab);
                *probability *= p_likely;
                visit(tab, Some(likely), *probability);
            } else {
                scratch.coeff_map.clear();
                for (value, bitstring) in tab.coefficients.take() {
                    scratch.coeff_map.insert(bitstring, value);
                }
                other_scratch.coeff_map.clone_from(&scratch.coeff_map);
                let mut other = tab.fork(Some(tab_seed));
                let mask = tab.odd_phase_destabilizer_mask();
                let overlap = GeneralizedTableau::<A, I, H>::compute_overlap_case_a(
                    &scratch.coeff_map,
                    phase,
                    destab,
                    stab,
                    mask,
                );
                let p_one = 0.5 - 0.5 * overlap;
                let likely = p_one > 0.5;
                let p_likely = if likely { p_one } else { 1.0 - p_one };
                let p_other = 1.0 - p_likely;
                let other_probability = *probability * p_other;
                if p_other > self.sum_cutoff {
                    other.project_case_a(!likely, &mut other_scratch, phase, stab, destab, qubit);
                    visit(&mut other, Some(!likely), other_probability);
                    let fp = fingerprint(&other);
                    branches.push((other, other_probability, fp));
                }
                tab.project_case_a(likely, &mut scratch, phase, stab, destab, qubit);
                *probability *= p_likely;
                visit(tab, Some(likely), *probability);
            }
        }

        self.mark_dirty();
        if self.insert_branches(branches) {
            self.normalize_probabilities();
        }
        self.truncate();
    }

    /// Analytic Z-basis measurement probabilities `(zero, one, lost)`.
    pub fn measure(&mut self, qubit: usize) -> (f64, f64, f64) {
        let mut probabilities = (Vec::new(), Vec::new(), Vec::new());
        self.for_each_z_branch(qubit, |_, outcome, probability| match outcome {
            Some(false) => probabilities.0.push(probability),
            Some(true) => probabilities.1.push(probability),
            None => probabilities.2.push(probability),
        });
        (
            probabilities.0.into_iter().sum(),
            probabilities.1.into_iter().sum(),
            probabilities.2.into_iter().sum(),
        )
    }
}

impl<A, I, H> Reset for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn reset(&mut self, qubit: usize) {
        self.for_each_z_branch(qubit, |tab, outcome, _| {
            if outcome == Some(true) {
                tab.x(qubit);
            }
        });
    }
}
