use ppvm_runtime::config::Config;
use ppvm_tableau::data::GeneralizedTableau;
use ppvm_tableau::measure::MeasureScratch;
use ppvm_tableau::sparsevec::SparseVector;
use ppvm_tableau::tableau_index::TableauIndex;

use crate::prelude::*;
use crate::storage::{EntryStore, phase_loss_hash, word_fingerprint};
use bitvec::view::BitView;
use num::PrimInt;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, ToPrimitive, Zero};
use rand::RngExt;
use std::fmt::Debug;

impl<T, I, C, S> GeneralizedTableauSum<T, I, C, S>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>
        + PartialOrd
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Debug + Send + Sync,
    S: EntryStore<T, I, C>,
{
    /// Routine for mid-circuit measurements (for obtaining results use sampling)
    /// branch into different outcomes and return probabilities for outcomes (zero, one, lost)
    pub fn measure(&mut self, addr0: usize) -> (T::Coeff, T::Coeff, T::Coeff) {
        let mut p_zero = Vec::<T::Coeff>::new();
        let mut p_one = Vec::<T::Coeff>::new();
        let mut p_lost = Vec::<T::Coeff>::new();
        let n_entries = self.entries.len();

        let mut branches =
            Vec::<(GeneralizedTableau<T, I, C>, T::Coeff, u64, u64)>::with_capacity(n_entries);

        let mut scratch = MeasureScratch::<I, T::Coeff>::new();
        let mut scratch_other_outcome = scratch.clone();

        self.entries
            .for_each_mut_with_keys(|tab, p_sum, word_fp, phase_loss_fp| {
                if tab.is_lost[addr0] {
                    // NOTE: deterministically outputs lost, no branching
                    p_lost.push(p_sum.clone());
                    return;
                }

                let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
                    tab.compute_decomposition(addr0, Pauli::Z);

                let tab_seed = self.rng.random::<u64>();

                if stab_anticomm_bits == I::zero() {
                    // case b
                    let entries: Vec<(Complex<T::Coeff>, I)> =
                        std::mem::replace(&mut tab.coefficients, C::new())
                            .into_iter()
                            .collect();

                    // NOTE: fork AFTER draining coefficients, so we only copy an
                    // empty coefficients vec
                    let mut tab_other_outcome = tab.fork(Some(tab_seed));

                    // Pass 1: compute overlap (read-only, real-only accumulation)
                    // Since conj(c)*c = |c|^2 (always real), the phase factor contribution
                    // to z_overlap.re is: phase 0 → +|c|^2, phase 2 → −|c|^2,
                    // phase 1,3 → 0 (imaginary × real = imaginary, doesn't contribute to .re)
                    let z_overlap_re = GeneralizedTableau::<T, I, C>::compute_overlap_case_b(
                        &entries,
                        phase_decomp,
                        destab_anticomm_bits,
                    );

                    let prob_1 = 0.5 - 0.5 * z_overlap_re;
                    let prob_0 = 1.0 - prob_1;
                    p_one.push(p_sum.clone() * prob_1);
                    p_zero.push(p_sum.clone() * prob_0);

                    // make the existing term the more likely outcome
                    let likely_outcome = prob_1 > 0.5;
                    let (p_likely, p_other) = if likely_outcome {
                        (prob_1, prob_0)
                    } else {
                        (prob_0, prob_1)
                    };

                    debug_assert!(
                        phase_decomp == 0 || phase_decomp == 2,
                        "Measurement result cannot be imaginary!"
                    );

                    // project
                    // NOTE: avoid projecting into zero amplitude state
                    // intentionally stricter than normal truncate
                    // dropping terms < sum_cutoff even if they'd merge with another one
                    if Into::<T::Coeff>::into(p_other) > self.sum_cutoff {
                        tab_other_outcome.project_case_b(
                            &entries,
                            !likely_outcome,
                            phase_decomp,
                            destab_anticomm_bits,
                        );
                        branches.push((
                            tab_other_outcome,
                            p_sum.clone() * p_other,
                            word_fp,
                            phase_loss_fp,
                        ));
                    }

                    // update existing entry
                    tab.project_case_b(
                        &entries,
                        likely_outcome,
                        phase_decomp,
                        destab_anticomm_bits,
                    );
                    *p_sum *= p_likely;
                } else {
                    // case a
                    scratch.coeff_map.clear();
                    scratch.coeff_map.reserve(tab.coefficients.len());
                    {
                        let coeff_map = &mut scratch.coeff_map;
                        tab.coefficients.retain(|(v, i)| {
                            coeff_map.insert(*i, *v);
                            false // drain — keeps allocation
                        });
                    }

                    scratch_other_outcome
                        .coeff_map
                        .clone_from(&scratch.coeff_map);

                    // NOTE: fork AFTER draining coefficients, so we only copy an
                    // empty coefficients vec
                    let mut tab_other_outcome = tab.fork(Some(tab_seed));

                    // Compute z_overlap.re directly (the imaginary part is always ~0).
                    // The mask is a pure function of destabilizer phases — cache it across
                    // measurements until `update_tableau_according_to_outcome` invalidates it.
                    let odd_phase_mask = *scratch
                        .odd_phase_mask
                        .get_or_insert_with(|| tab.odd_phase_destabilizer_mask());
                    let z_overlap_re = GeneralizedTableau::<T, I, C>::compute_overlap_case_a(
                        &scratch.coeff_map,
                        phase_decomp,
                        destab_anticomm_bits,
                        stab_anticomm_bits,
                        odd_phase_mask,
                    );

                    let prob_1 = 0.5 - 0.5 * z_overlap_re;
                    let prob_0 = 1.0 - prob_1;

                    // make the existing term the more likely outcome
                    let likely_outcome = prob_1 > 0.5;
                    let (p_likely, p_other) = if likely_outcome {
                        (prob_1, prob_0)
                    } else {
                        (prob_0, prob_1)
                    };

                    // project
                    if Into::<T::Coeff>::into(p_other) > self.sum_cutoff {
                        tab_other_outcome.project_case_a(
                            !likely_outcome,
                            &mut scratch_other_outcome,
                            phase_decomp,
                            stab_anticomm_bits,
                            destab_anticomm_bits,
                            addr0,
                        );
                        let word_fp_other = word_fingerprint(&tab_other_outcome);
                        let phase_loss_other = phase_loss_hash(&tab_other_outcome);
                        branches.push((
                            tab_other_outcome,
                            p_sum.clone() * p_other,
                            word_fp_other,
                            phase_loss_other,
                        ));
                    }

                    // update exisiting entry in-place
                    tab.project_case_a(
                        likely_outcome,
                        &mut scratch,
                        phase_decomp,
                        stab_anticomm_bits,
                        destab_anticomm_bits,
                        addr0,
                    );
                    *p_sum *= p_likely;
                }
            });

        self.entries.mark_keys_dirty();
        let needs_normalize = self
            .entries
            .insert_or_merge_batch(branches, &self.sum_cutoff);
        if needs_normalize {
            self.normalize_probabilities();
        }
        self.truncate();

        let p_0 = p_zero
            .iter()
            .fold(T::Coeff::zero(), |acc, p| acc + p.clone());
        let p_1 = p_one
            .iter()
            .fold(T::Coeff::zero(), |acc, p| acc + p.clone());
        let p_l = p_lost
            .iter()
            .fold(T::Coeff::zero(), |acc, p| acc + p.clone());
        (p_0, p_1, p_l)
    }
}
