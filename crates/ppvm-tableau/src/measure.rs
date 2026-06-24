// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::data::{compute_phase_with_mask_static, symplectic_inner};
use crate::prelude::*;
use bitvec::view::BitView;
use fxhash::FxHashMap as HashMap;
use num::PrimInt;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, ToPrimitive, Zero};
use rand::RngExt;
use std::fmt::Debug;

impl<T: Config> Measure for Tableau<T>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
{
    /// Measure qubit `addr0` in Z basis
    fn measure(&mut self, addr0: usize) -> bool {
        let q = self.find_z_anticommuting_stabilizer(addr0);
        match q {
            Some(q_idx) => {
                // Case a: random measurement outcome
                // At least one stabilizer anticommutes with Z_addr0

                // Generate random measurement outcome (50/50)
                let outcome = self.rng.random::<bool>();

                self.update_tableau_according_to_outcome(addr0, q_idx, outcome);

                outcome
            }
            None => {
                // Case b: deterministic measurement outcome

                self.get_deterministic_outcome(addr0)
            }
        }
    }
}

const COMPLEX_PHASE_CONVERSION: [Complex64; 4] = [
    Complex64::new(1.0, 0.0),  // +1
    Complex64::new(0.0, 1.0),  // +i
    Complex64::new(-1.0, 0.0), // -1
    Complex64::new(0.0, -1.0), // -i
];

/// Per-measurement scratch buffers, reused across qubits within a single
/// `measure_all` invocation — and, when threaded through
/// [`measure_all_with_scratch`](GeneralizedTableau::measure_all_with_scratch),
/// across many shots of a sampler too. Reusing one scratch keeps the case-a
/// HashMap and the b-entries Vec out of the per-shot allocator churn.
///
/// - `odd_phase_mask` is lazily computed and cached until the destabilizers
///   change (i.e. until a case-a measurement runs `update_tableau_according_to_outcome`).
/// - `coeff_map` is the case-a HashMap holding `(idx → amplitude)` between
///   the overlap, partition, and merge passes.
/// - `b_entries` is the case-a partition's "k-bit = 1" scratch Vec.
///
/// Construct one per active sampling thread; the type is not meant to be
/// shared across threads concurrently.
#[derive(Clone)]
pub struct MeasureScratch<I, R> {
    pub odd_phase_mask: Option<I>,
    pub coeff_map: HashMap<I, Complex<R>>,
    pub b_entries: Vec<(I, Complex<R>)>,
}

impl<I, R> MeasureScratch<I, R> {
    pub fn new() -> Self {
        Self {
            odd_phase_mask: None,
            coeff_map: HashMap::default(),
            b_entries: Vec::new(),
        }
    }
}

impl<I, R> Default for MeasureScratch<I, R> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> LossyMeasure
    for GeneralizedTableau<T, I, C>
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
        + PartialOrd<f64>,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Debug,
{
    fn measure(&mut self, addr0: usize) -> Option<bool> {
        if self.is_lost[addr0] {
            self.measurement_record.push(None);
            return None;
        }

        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            self.compute_decomposition(addr0, Pauli::Z);

        // Standalone callers don't get cross-call cache benefits; `measure_all`
        // threads through a long-lived scratch.
        let mut scratch = MeasureScratch::new();
        self.measure_with_scratch(
            addr0,
            &mut scratch,
            phase_decomp,
            stab_anticomm_bits,
            destab_anticomm_bits,
        )
    }

    /// Override the trait default (a per-target `measure` loop, which allocates
    /// a fresh `MeasureScratch` on every call) with a single scratch reused
    /// across the whole batch, amortizing the case-a HashMap / `b_entries`
    /// allocations and the cached odd-phase-destabilizer mask. Outcomes, the
    /// measurement record, and the RNG-draw order are identical to measuring
    /// each target individually — only the internal allocation pattern changes.
    fn measure_many(&mut self, targets: &[usize]) -> Vec<Option<bool>> {
        let mut scratch = MeasureScratch::new();
        self.measure_many_with_scratch(targets, &mut scratch)
    }
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> GeneralizedTableau<T, I, C>
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
        + PartialOrd<f64>,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Debug,
{
    pub(crate) fn measure_with_scratch(
        &mut self,
        addr0: usize,
        scratch: &mut MeasureScratch<I, T::Coeff>,
        phase_decomp: u8,
        stab_anticomm_bits: I,
        destab_anticomm_bits: I,
    ) -> Option<bool> {
        if stab_anticomm_bits == I::zero() {
            // Case b (fast path): Z is already a stabilizer. Overlap + filter in place.

            // Compute overlap directly on self.coefficients without draining.
            let mut z_overlap_re = 0.0f64;
            for &(coeff, idx) in self.coefficients.iter() {
                let phase = (phase_decomp
                    + (2 * symplectic_inner(destab_anticomm_bits, idx) as u8) % 4)
                    % 4;
                if !phase.is_multiple_of(2) {
                    continue;
                }
                let norm_sq = coeff.norm_sqr().to_f64().unwrap_or(0.0);
                if phase == 0 {
                    z_overlap_re += norm_sq;
                } else {
                    z_overlap_re -= norm_sq;
                }
            }

            let prob_1 = 0.5 - 0.5 * z_overlap_re;
            let outcome = self.tableau.rng.random::<f64>() < prob_1;

            debug_assert!(
                phase_decomp == 0 || phase_decomp == 2,
                "Measurement result cannot be imaginary!"
            );

            let old_len = self.coefficients.len();
            let z_sign = phase_decomp == 2;
            self.coefficients.retain(|(_, alpha)| {
                let parity = symplectic_inner(*alpha, destab_anticomm_bits) % 2 != 0;
                (parity ^ outcome) == z_sign
            });
            if self.coefficients.len() < old_len {
                self.coefficients.normalize();
            }

            // Case-b doesn't mutate destabilizers, so the cached mask remains valid.
            self.measurement_record.push(Some(outcome));
            Some(outcome)
        } else {
            // Case a: Z is not a stabilizer — sort-merge instead of HashMap.
            let mut by_idx: Vec<(I, Complex<T::Coeff>)> =
                std::mem::replace(&mut self.coefficients, C::new())
                    .into_iter()
                    .map(|(c, i)| (i, c))
                    .collect();
            {
                let mut sorted = true;
                let mut prev_k: Option<I> = None;
                for &(k, _) in &by_idx {
                    if let Some(p) = prev_k
                        && k < p
                    {
                        sorted = false;
                        break;
                    }
                    prev_k = Some(k);
                }
                if !sorted {
                    by_idx.sort_unstable_by_key(|a| a.0);
                }
            }

            let odd_phase_mask = *scratch
                .odd_phase_mask
                .get_or_insert_with(|| self.odd_phase_destabilizer_mask());

            // OVERLAP: 2-way merge of by_idx and shifted (each entry XOR'd by stab_anticomm_bits).
            // At equal key k: by_idx has (k, a_k) and shifted has (k, a_{k^s}), matching
            // the HashMap overlap's coeff / coeff_branch pair — counted once per key from
            // each side, exactly as the HashMap iterates every (idx, branch_index) pair.
            let mut shifted: Vec<(I, Complex<T::Coeff>)> = by_idx
                .iter()
                .map(|&(i, c)| (i ^ stab_anticomm_bits, c))
                .collect();
            shifted.sort_unstable_by_key(|a| a.0);

            let mut z_overlap_re = 0.0f64;
            {
                let mut ii = 0usize;
                let mut jj = 0usize;
                while ii < by_idx.len() && jj < shifted.len() {
                    match by_idx[ii].0.cmp(&shifted[jj].0) {
                        std::cmp::Ordering::Less => {
                            ii += 1;
                        }
                        std::cmp::Ordering::Greater => {
                            jj += 1;
                        }
                        std::cmp::Ordering::Equal => {
                            let (idx, a) = by_idx[ii];
                            let (_, b) = shifted[jj];
                            let phase = (phase_decomp
                                + compute_phase_with_mask_static(
                                    destab_anticomm_bits,
                                    idx,
                                    stab_anticomm_bits,
                                    odd_phase_mask,
                                ))
                                % 4;
                            let a_re = a.re.to_f64().unwrap_or(0.0);
                            let a_im = a.im.to_f64().unwrap_or(0.0);
                            let b_re = b.re.to_f64().unwrap_or(0.0);
                            let b_im = b.im.to_f64().unwrap_or(0.0);
                            let re_w = a_re * b_re + a_im * b_im;
                            let im_w = a_re * b_im - a_im * b_re;
                            match phase {
                                0 => z_overlap_re += re_w,
                                1 => z_overlap_re += im_w,
                                2 => z_overlap_re -= re_w,
                                3 => z_overlap_re -= im_w,
                                _ => unreachable!(),
                            }
                            ii += 1;
                            jj += 1;
                        }
                    }
                }
            }

            let prob_1 = 0.5 - 0.5 * z_overlap_re;
            let outcome = self.tableau.rng.random::<f64>() < prob_1;

            // PROJECTION: partition A (k-bit=0) and B (k-bit=1), transform B, merge.
            let q_idx = stab_anticomm_bits.trailing_zeros() as usize;
            let k = I::one() << q_idx;
            let alpha = if outcome {
                (phase_decomp + 2) % 4
            } else {
                phase_decomp
            };

            let mut a: Vec<(I, Complex<T::Coeff>)> = Vec::new();
            let mut bt: Vec<(I, Complex<T::Coeff>)> = Vec::new();
            for (idx, coeff) in by_idx {
                if (idx & k) == I::zero() {
                    a.push((idx, coeff));
                } else {
                    let symp = symplectic_inner(idx, destab_anticomm_bits);
                    let phase_idx =
                        ((alpha as i32 + if symp % 2 == 1 { 2 } else { 0 }) % 4) as usize;
                    let q: Complex<T::Coeff> = COMPLEX_PHASE_CONVERSION[phase_idx].into();
                    bt.push((idx ^ stab_anticomm_bits, q * coeff));
                }
            }
            // `a` is already sorted (subset of sorted by_idx); bt needs sorting.
            bt.sort_unstable_by_key(|e| e.0);

            // 2-way merge summing equal keys → sorted merged output.
            let mut merged: Vec<(I, Complex<T::Coeff>)> = Vec::with_capacity(a.len() + bt.len());
            {
                let mut i = 0usize;
                let mut j = 0usize;
                while i < a.len() && j < bt.len() {
                    match a[i].0.cmp(&bt[j].0) {
                        std::cmp::Ordering::Less => {
                            merged.push(a[i]);
                            i += 1;
                        }
                        std::cmp::Ordering::Greater => {
                            merged.push(bt[j]);
                            j += 1;
                        }
                        std::cmp::Ordering::Equal => {
                            let mut sv = a[i].1;
                            sv += bt[j].1;
                            merged.push((a[i].0, sv));
                            i += 1;
                            j += 1;
                        }
                    }
                }
                while i < a.len() {
                    merged.push(a[i]);
                    i += 1;
                }
                while j < bt.len() {
                    merged.push(bt[j]);
                    j += 1;
                }
            }

            let norm_sqr = merged
                .iter()
                .fold(T::Coeff::zero(), |acc, (_, c)| acc + c.norm_sqr());
            let cutoff_sq = self.coefficient_threshold.clone() * self.coefficient_threshold.clone();
            let threshold = cutoff_sq.to_f64().unwrap_or(0.0) * norm_sqr.to_f64().unwrap_or(0.0);
            self.coefficients.reserve(merged.len());
            for (idx, coeff) in merged {
                if coeff.norm_sqr() > threshold {
                    self.coefficients.unsafe_insert(idx, coeff);
                }
            }

            self.coefficients.normalize();
            self.tableau
                .update_tableau_according_to_outcome(addr0, q_idx, outcome);
            scratch.odd_phase_mask = None;
            self.measurement_record.push(Some(outcome));
            Some(outcome)
        }
    }

    pub fn project_case_a(
        &mut self,
        outcome: bool,
        scratch: &mut MeasureScratch<I, T::Coeff>,
        phase_decomp: u8,
        stab_anticomm_bits: I,
        destab_anticomm_bits: I,
        addr0: usize,
    ) {
        // Case a: Z is not a stabilizer — need HashMap for cross-index lookups.
        // Drain self.coefficients into scratch.coeff_map via `retain` so the
        // Vec's capacity survives and we can refill it at the end without a
        // fresh allocation.

        let q_idx = stab_anticomm_bits.trailing_zeros() as usize;

        let one = I::one();
        let zero = I::zero();
        let k = one << q_idx;

        let alpha = if outcome {
            (phase_decomp + 2) % 4
        } else {
            phase_decomp
        };

        // Partition into A (k-bit=0) and B (k-bit=1) via retain, then merge.
        // Split the borrow so `retain` can mutate coeff_map while the closure
        // pushes into b_entries.
        scratch.b_entries.clear();
        let MeasureScratch {
            coeff_map,
            b_entries,
            ..
        } = scratch;
        b_entries.reserve(coeff_map.len() / 2 + 1);
        coeff_map.retain(|idx, coeff| {
            if (*idx & k) != zero {
                b_entries.push((*idx, *coeff));
                false // remove B entry
            } else {
                true // keep A entry
            }
        });
        // Merge B entries into their A partners with phase adjustment.
        Self::merge_b_into_a(
            coeff_map,
            b_entries,
            alpha,
            destab_anticomm_bits,
            stab_anticomm_bits,
        );

        // Keep entries where |c|/norm > threshold.
        let norm_sqr = coeff_map
            .values()
            .fold(T::Coeff::zero(), |acc, c: &Complex<T::Coeff>| {
                acc + c.norm_sqr()
            });

        let cutoff_sq = self.coefficient_threshold.clone() * self.coefficient_threshold.clone();
        let threshold = cutoff_sq.to_f64().unwrap_or(0.0) * norm_sqr.to_f64().unwrap_or(0.0);
        // self.coefficients is already empty here (drained via retain above);
        // reserve is mostly a no-op since the prior capacity is still there.
        self.coefficients.reserve(coeff_map.len());
        for (idx, coeff) in coeff_map.drain() {
            if coeff.norm_sqr() > threshold {
                self.coefficients.unsafe_insert(idx, coeff);
            }
        }

        self.coefficients.normalize();

        self.tableau
            .update_tableau_according_to_outcome(addr0, q_idx, outcome);
        // Destabilizer phases just changed, invalidate the cached mask.
        scratch.odd_phase_mask = None;
    }

    /// project state in case b (Z is a stabilizer) according to sampled outcome
    pub fn project_case_b(
        &mut self,
        entries: &[(Complex<T::Coeff>, I)],
        outcome: bool,
        phase_decomp: u8,
        destab_anticomm_bits: I,
    ) {
        let old_len = entries.len();

        let z_sign = phase_decomp == 2;

        // Pass 2: filter directly into self.coefficients (no retain needed)
        self.coefficients.reserve(entries.len());
        for &(coeff, alpha) in entries {
            let parity = symplectic_inner(alpha, destab_anticomm_bits) % 2 != 0;
            if (parity ^ outcome) == z_sign {
                self.coefficients.unsafe_insert(alpha, coeff);
            }
        }

        if self.coefficients.len() < old_len {
            self.coefficients.normalize();
        }
    }
}

impl<T, I, C> GeneralizedTableau<T, I, C>
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
        + PartialOrd<f64>,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Debug,
{
    /// Measure qubit `addr0` in Z basis with optional readout noise.
    ///
    /// Behaves like [`measure`](LossyMeasure::measure), then with probability
    /// `flip_prob` flips the *recorded* bit. The qubit's quantum state stays
    /// consistent with the true outcome — only the returned value flips.
    /// `flip_prob = 0.0` is equivalent to `measure`.
    ///
    /// If the qubit is lost, returns `None` regardless of `flip_prob`.
    pub fn measure_noisy(&mut self, addr0: usize, flip_prob: f64) -> Option<bool> {
        debug_assert!(
            (0.0..=1.0).contains(&flip_prob),
            "flip_prob must be in [0, 1], got {flip_prob}"
        );
        // `measure` already pushed the (un-flipped) outcome onto the record.
        // Overwrite that last entry with the post-noise value so exactly one
        // push occurs per logical measurement and the record matches what we
        // return.
        let outcome = self.measure(addr0)?;
        let noisy = self.flip_with_prob(outcome, flip_prob);
        self.overwrite_last_measurement_record(Some(noisy));
        Some(noisy)
    }

    /// Sample a Bernoulli(`p`) outcome using the tableau's internal RNG.
    /// Used by Stim measurement-noise dispatch in `ppvm-stim`.
    pub fn bernoulli(&mut self, p: f64) -> bool {
        debug_assert!((0.0..=1.0).contains(&p), "p must be in [0, 1], got {p}");
        self.tableau.rng.random::<f64>() < p
    }

    /// Flip `bit` with probability `p`. Used by Stim MR/MPad readout-noise
    /// dispatch in `ppvm-stim`. Returns `bit` unchanged when `p <= 0.0`.
    pub fn flip_with_prob(&mut self, bit: bool, p: f64) -> bool {
        debug_assert!((0.0..=1.0).contains(&p), "p must be in [0, 1], got {p}");
        if p > 0.0 && self.bernoulli(p) {
            !bit
        } else {
            bit
        }
    }
}

/// Measurement overlap helper functions, with optional rayon parallelism.
impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> GeneralizedTableau<T, I, C>
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
        + PartialOrd<f64>,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Debug,
{
    /// `⟨Z⟩` on qubit `addr0`, computed non-destructively (the state is not
    /// collapsed). Reuses the measurement overlap machinery; cost scales with
    /// the number of coefficients (and n²).
    pub fn z_expectation(&self, addr0: usize) -> f64 {
        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            self.compute_decomposition(addr0, Pauli::Z);

        if stab_anticomm_bits == I::zero() {
            // Case b: Z is a stabilizer — self-pairing overlap.
            let entries: Vec<(Complex<T::Coeff>, I)> = self.coefficients.iter().copied().collect();
            Self::compute_overlap_case_b(&entries, phase_decomp, destab_anticomm_bits)
        } else {
            // Case a: cross-index pairing — clone coefficients into a map (read-only,
            // so unlike `measure` we don't drain `self.coefficients`).
            let coeff_map: HashMap<I, Complex<T::Coeff>> =
                self.coefficients.iter().map(|&(c, i)| (i, c)).collect();
            let odd_phase_mask = self.odd_phase_destabilizer_mask();
            Self::compute_overlap_case_a(
                &coeff_map,
                phase_decomp,
                destab_anticomm_bits,
                stab_anticomm_bits,
                odd_phase_mask,
            )
        }
    }

    /// Case_b overlap: self-pairing (branch_index = idx), so overlap = ±|c|^2.
    /// Only even phases contribute to the real part.
    pub fn compute_overlap_case_b(
        entries: &[(Complex<T::Coeff>, I)],
        phase_decomp: u8,
        destab_anticomm_bits: I,
    ) -> f64 {
        let mut z_overlap_re = 0.0f64;
        for &(coeff, idx) in entries {
            let phase =
                (phase_decomp + (2 * symplectic_inner(destab_anticomm_bits, idx) as u8) % 4) % 4;
            if !phase.is_multiple_of(2) {
                continue;
            }
            let norm_sq = coeff.norm_sqr().to_f64().unwrap_or(0.0);
            if phase == 0 {
                z_overlap_re += norm_sq;
            } else {
                z_overlap_re -= norm_sq;
            }
        }
        z_overlap_re
    }

    /// Case_a overlap: cross-index pairing via HashMap lookup.
    /// Accumulates only the real part of z_overlap.
    pub fn compute_overlap_case_a(
        coeff_map: &HashMap<I, Complex<T::Coeff>>,
        phase_decomp: u8,
        destab_anticomm_bits: I,
        stab_anticomm_bits: I,
        odd_phase_mask: I,
    ) -> f64 {
        let mut z_overlap_re = 0.0f64;
        for (&idx, coeff) in coeff_map {
            let branch_index = idx ^ stab_anticomm_bits;
            let phase = (phase_decomp
                + compute_phase_with_mask_static(
                    destab_anticomm_bits,
                    idx,
                    stab_anticomm_bits,
                    odd_phase_mask,
                ))
                % 4;
            let Some(coeff_branch) = coeff_map.get(&branch_index).copied() else {
                continue;
            };
            let a_re = coeff.re.to_f64().unwrap_or(0.0);
            let a_im = coeff.im.to_f64().unwrap_or(0.0);
            let b_re = coeff_branch.re.to_f64().unwrap_or(0.0);
            let b_im = coeff_branch.im.to_f64().unwrap_or(0.0);
            let re_w = a_re * b_re + a_im * b_im;
            let im_w = a_re * b_im - a_im * b_re;
            match phase {
                0 => z_overlap_re += re_w,
                1 => z_overlap_re += im_w,
                2 => z_overlap_re -= re_w,
                3 => z_overlap_re -= im_w,
                _ => unreachable!(),
            }
        }
        z_overlap_re
    }

    /// Merge B entries (k-bit=1) into their A counterparts in coeff_map.
    /// With rayon: parallel phase computation, sequential HashMap accumulation.
    fn merge_b_into_a(
        coeff_map: &mut HashMap<I, Complex<T::Coeff>>,
        b_entries: &[(I, Complex<T::Coeff>)],
        alpha: u8,
        destab_anticomm_bits: I,
        stab_anticomm_bits: I,
    ) {
        for &(idx, coeff) in b_entries {
            let symp_inner = symplectic_inner(idx, destab_anticomm_bits);
            let phase_idx = ((alpha as i32 + if symp_inner % 2 == 1 { 2 } else { 0 }) % 4) as usize;
            let q: Complex<T::Coeff> = COMPLEX_PHASE_CONVERSION[phase_idx].into();
            *coeff_map
                .entry(idx ^ stab_anticomm_bits)
                .or_insert(Complex::zero()) += q * coeff;
        }
    }
}
