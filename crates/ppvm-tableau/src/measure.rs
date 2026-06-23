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
            // Case b (fast path): Z is already a stabilizer.
            // branch_index = idx ^ 0 = idx, so the overlap is self-pairing:
            //   z_overlap = sum_alpha phase_factor(alpha).conj() * |c_alpha|^2
            // No HashMap needed — work directly on the Vec.

            // Collect into a Vec so we can iterate twice: once for overlap,
            // once for filtering. The collect from Vec IntoIter→Vec is a no-op.
            let entries: Vec<(Complex<T::Coeff>, I)> =
                std::mem::replace(&mut self.coefficients, C::new())
                    .into_iter()
                    .collect();

            // Pass 1: compute overlap (read-only, real-only accumulation)
            // Since conj(c)*c = |c|^2 (always real), the phase factor contribution
            // to z_overlap.re is: phase 0 → +|c|^2, phase 2 → −|c|^2,
            // phase 1,3 → 0 (imaginary × real = imaginary, doesn't contribute to .re)
            let z_overlap_re =
                Self::compute_overlap_case_b(&entries, phase_decomp, destab_anticomm_bits);

            let prob_1 = 0.5 - 0.5 * z_overlap_re;
            let outcome = self.tableau.rng.random::<f64>() < prob_1;

            debug_assert!(
                phase_decomp == 0 || phase_decomp == 2,
                "Measurement result cannot be imaginary!"
            );

            self.project_case_b(&entries, outcome, phase_decomp, destab_anticomm_bits);

            // Case-b doesn't mutate destabilizers, so the cached mask remains valid.
            self.measurement_record.push(Some(outcome));
            Some(outcome)
        } else {
            // Case a: Z is not a stabilizer — need HashMap for cross-index lookups.
            // Drain self.coefficients into scratch.coeff_map via `retain` so the
            // Vec's capacity survives and we can refill it at the end without a
            // fresh allocation.
            scratch.coeff_map.clear();
            scratch.coeff_map.reserve(self.coefficients.len());
            {
                let coeff_map = &mut scratch.coeff_map;
                self.coefficients.retain(|(v, i)| {
                    coeff_map.insert(*i, *v);
                    false // drain — keeps allocation
                });
            }

            // Compute z_overlap.re directly (the imaginary part is always ~0).
            // The mask is a pure function of destabilizer phases — cache it across
            // measurements until `update_tableau_according_to_outcome` invalidates it.
            let odd_phase_mask = *scratch
                .odd_phase_mask
                .get_or_insert_with(|| self.odd_phase_destabilizer_mask());
            let z_overlap_re = Self::compute_overlap_case_a(
                &scratch.coeff_map,
                phase_decomp,
                destab_anticomm_bits,
                stab_anticomm_bits,
                odd_phase_mask,
            );

            let prob_1 = 0.5 - 0.5 * z_overlap_re;
            let outcome = self.tableau.rng.random::<f64>() < prob_1;
            self.project_case_a(
                outcome,
                scratch,
                phase_decomp,
                stab_anticomm_bits,
                destab_anticomm_bits,
                addr0,
            );
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
