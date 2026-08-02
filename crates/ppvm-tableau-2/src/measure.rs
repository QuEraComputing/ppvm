// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Z-basis measurement: the pure-Clifford frame algorithm, the
//! coefficient-aware `O(n²)` algorithm, the reusable [`MeasureScratch`], and the
//! batched `measure_all` / `measure_many` entry points.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Behavioral traits"
//! ([`Measure`] shares the `Option<bool>` result type across both tableau types
//! but *not* the algorithm) and §"Pauli algebra traits" (measurement is built on
//! the [`StabilizerFrame`](ppvm_traits_2::StabilizerFrame) primitives).
//!
//! The dichotomy the pivot search rests on — the outcome is deterministic
//! exactly when the measured Pauli commutes with every stabilizer — is
//! machine-checked in `lean/PPVM/Tableau/Frame.lean` (`measurement_dichotomy`,
//! `measure_deterministic_iff_xfree`), and the case-a frame projection is proved
//! to keep the `2n` rows a symplectic basis there too
//! (`isSymplecticFrame_projectFrame`).
//!
//! What the case-a *arithmetic* computes is machine-checked in
//! `lean/PPVM/Tableau/Projection.lean`, which models the frame-conjugated `Z_q`
//! as `M|k⟩ = i^{φ k}|k ⊕ s⟩`:
//!
//! - `rustTerm_eq` — the overlap merge's four-way ℤ/4 dispatch
//!   (`0 => +re_w, 1 => +im_w, 2 => -re_w, 3 => -im_w`) is exactly
//!   `Re(conj(i^φ · a) · b)`. Note it is **not** `Re(i^φ · conj a · b)`: the odd
//!   branches carry the conjugated phase.
//! - `shiftOp_involutive` / `shiftOp_selfAdjoint` — `M² = I`, `M† = M`.
//! - `overlap_eq_inner` — `z_overlap_re` is `Re⟨c, M c⟩`.
//! - `proj_add` / `proj_idem` — `P₀ + P₁ = I`, `P_b² = P_b` for
//!   `P_b = (I + (−1)^b M)/2`.
//! - `probOne_eq` — `prob_1 = 0.5 − 0.5 · z_overlap_re` **is** the Born
//!   probability `⟨c, P₁ c⟩`.
//! - `projectRaw_eq_two_proj` — the keep-`A`/transform-`B`/merge below is
//!   `2 · P_b` on the surviving half; the unconditional `normalize()` supplies
//!   the missing factor.
//!
//! Scope note (see that file): the projection's phase `alpha + 2·⟨idx, destab⟩`
//! omits the odd-phase-destabilizer term the overlap folds in through
//! `compute_phase_with_mask_static`, because the two are read in different
//! frames (pre- and post-projection). Relating them needs a Hilbert-space model
//! of the frame that the Lean development does not have; old and this crate
//! agree verbatim, so it is a specification gap, not a port divergence.
//!
//! Ported from `ppvm-tableau/src/{measure,measure_all}.rs`.

use fxhash::FxHashMap as HashMap;
use num::complex::Complex64;
use ppvm_traits_2::{Clifford, Measure, Pauli, Reset};
use rand::RngExt;

use crate::data::{
    Bitstring, COMPLEX_PHASE_CONVERSION, GeneralizedTableau, RowStorage, Tableau,
    compute_phase_with_mask_static, symplectic_inner,
};

/// The pure Clifford measurement procedure.
///
/// # Behaviour note
///
/// The old crate split this as `Measure -> bool` (bare tableau) and
/// `LossyMeasure -> Option<bool>` (generalized tableau). The `-2` design removes
/// that split: both return `Option<bool>`, with `None` reserved for a lost qubit
/// (design §"Behavioral traits"). A bare frame carries no loss bits, so this
/// impl returns `Some(_)` unconditionally — the same information, one type.
///
/// Keeps no measurement record and never normalizes. Case a consumes exactly one
/// `random::<bool>()`; case b consumes **no** randomness.
impl<A: RowStorage, H> Measure for Tableau<A, H> {
    fn measure(&mut self, qubit: usize) -> Option<bool> {
        match self.find_z_anticommuting_stabilizer(qubit) {
            Some(q_idx) => {
                // Case a: at least one stabilizer anticommutes with `Z_qubit`,
                // so the outcome is a fair coin.
                let outcome = self.rng.random::<bool>();
                self.update_tableau_according_to_outcome(qubit, q_idx, outcome);
                Some(outcome)
            }
            // Case b: deterministic outcome, no RNG draw.
            None => Some(self.get_deterministic_outcome(qubit)),
        }
    }
}

/// Per-measurement scratch buffers, reused across qubits within a single
/// `measure_all` / `measure_many` — and, when threaded through
/// [`GeneralizedTableau::measure_all_with_scratch`], across many shots of a
/// sampler too.
///
/// - `odd_phase_mask` is lazily computed and cached until the destabilizer
///   phases change, i.e. until a **case-a** projection runs
///   `update_tableau_according_to_outcome`. A case-b measurement does not touch
///   the frame, so the cache survives it.
/// - `coeff_map` / `b_entries` back the map-based
///   [`GeneralizedTableau::project_case_a`].
/// - `by_idx` / `shifted` / `a` / `bt` / `merged` back the default sort-merge
///   case-a path. The old crate allocated these five fresh on **every** case-a
///   measurement even though the scratch existed to hold them; owning them here
///   is allocation-only — the merge order and arithmetic are unchanged, so the
///   result is byte-identical.
///
/// Construct one per active sampling thread; not meant to be shared
/// concurrently.
#[derive(Clone, Debug)]
pub struct MeasureScratch<I> {
    /// Cached odd-phase destabilizer mask; `None` means "recompute".
    pub odd_phase_mask: Option<I>,
    /// Case-a map working set (`project_case_a`).
    pub coeff_map: HashMap<I, Complex64>,
    /// Case-a partition scratch: the `k`-bit = 1 entries.
    pub b_entries: Vec<(I, Complex64)>,
    by_idx: Vec<(I, Complex64)>,
    shifted: Vec<(I, Complex64)>,
    a: Vec<(I, Complex64)>,
    bt: Vec<(I, Complex64)>,
    merged: Vec<(I, Complex64)>,
}

impl<I> MeasureScratch<I> {
    /// A fresh, empty scratch.
    pub fn new() -> Self {
        Self {
            odd_phase_mask: None,
            coeff_map: HashMap::default(),
            b_entries: Vec::new(),
            by_idx: Vec::new(),
            shifted: Vec::new(),
            a: Vec::new(),
            bt: Vec::new(),
            merged: Vec::new(),
        }
    }
}

impl<I> Default for MeasureScratch<I> {
    fn default() -> Self {
        Self::new()
    }
}

impl<A: RowStorage, I: Bitstring, H> Measure for GeneralizedTableau<A, I, H> {
    /// Measure `qubit` in the Z basis.
    ///
    /// A lost qubit pushes `None` onto the measurement record and returns `None`
    /// **without** touching the state. Otherwise exactly one `Some(outcome)` is
    /// pushed.
    ///
    /// The working buffers come from the tableau's own scratch rather than from
    /// a fresh [`MeasureScratch`] per call, so a `for q in .. { tab.measure(q) }`
    /// loop — the shape the MSD readout and the CNOT-chain measurement sweep use
    /// — does not re-allocate the five case-a sort-merge `Vec`s on every qubit.
    /// This is allocation only: `with_scratch` resets the cached odd-phase mask,
    /// so every call sees exactly the `None` a freshly constructed scratch would
    /// have given it, and the merge order and arithmetic are untouched.
    fn measure(&mut self, qubit: usize) -> Option<bool> {
        if self.is_lost[qubit] {
            self.measurement_record.push(None);
            return None;
        }

        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            self.compute_decomposition(qubit, Pauli::Z);

        self.with_scratch(|s, scratch| {
            s.measure_with_scratch(
                qubit,
                scratch,
                phase_decomp,
                stab_anticomm_bits,
                destab_anticomm_bits,
            )
        })
    }

    /// Override the trait default (a per-target `measure` loop) with one scratch
    /// held across the whole batch, so the sort-merge buffers are sized once.
    /// Outcomes, the measurement record, and the RNG-draw **order** are identical
    /// to measuring each target individually; only the internal allocation
    /// pattern changes.
    fn measure_many(&mut self, targets: &[usize]) -> Vec<Option<bool>> {
        self.with_scratch(|s, scratch| s.measure_many_with_scratch(targets, scratch))
    }
}

impl<A: RowStorage, I: Bitstring, H> GeneralizedTableau<A, I, H> {
    /// Run `f` with the tableau's own measurement working set.
    ///
    /// The scratch is moved out for the duration (so `f` gets `&mut self` too)
    /// and put back afterwards, keeping the buffer capacities alive for the next
    /// call. It is boxed, so "moved out" is one pointer, not a copy of the whole
    /// multi-`Vec` struct, and a tableau that has never been measured allocates
    /// nothing. Only capacity survives — the cached odd-phase mask is cleared on
    /// the way in, so every entry point sees exactly the `None` that old's
    /// freshly constructed `MeasureScratch` gave it.
    ///
    /// That clear is belt-and-braces rather than load-bearing: the mask is
    /// populated only inside the case-a branch and every path that populates it
    /// (`measure_with_scratch`'s case a, `project_case_a`) resets it to `None`
    /// before returning, so it is already `None` on entry. It is written
    /// explicitly because *without* it the safety of carrying the scratch across
    /// calls would depend on that non-local invariant — and a mask surviving a
    /// Clifford that flipped a destabilizer phase would silently mis-phase every
    /// subsequent case-a coefficient.
    #[inline]
    fn with_scratch<R>(&mut self, f: impl FnOnce(&mut Self, &mut MeasureScratch<I>) -> R) -> R {
        let mut scratch = self
            .scratch
            .take()
            .unwrap_or_else(|| Box::new(MeasureScratch::new()));
        scratch.odd_phase_mask = None;
        let out = f(self, &mut scratch);
        self.scratch = Some(scratch);
        out
    }

    /// Measure every qubit `0..n` in ascending order, reusing one scratch.
    pub fn measure_all(&mut self) -> Vec<Option<bool>> {
        self.with_scratch(|s, scratch| s.measure_all_with_scratch(scratch))
    }

    /// [`Self::measure_all`] with a caller-supplied scratch, reused across the
    /// `n` per-qubit measurements (and, if the caller chooses, across shots).
    ///
    /// This is the entry point a sampler should use: initialize one scratch
    /// alongside the sampler and thread it through every shot to amortize the
    /// case-a working set.
    pub fn measure_all_with_scratch(
        &mut self,
        scratch: &mut MeasureScratch<I>,
    ) -> Vec<Option<bool>> {
        (0..self.n_qubits())
            .map(|idx| self.measure_one_with_scratch(idx, scratch))
            .collect()
    }

    /// Measure the given qubit `indices` **in the caller's order**, reusing a
    /// caller-supplied scratch — the explicit-index analogue of
    /// [`Self::measure_all_with_scratch`].
    pub fn measure_many_with_scratch(
        &mut self,
        indices: &[usize],
        scratch: &mut MeasureScratch<I>,
    ) -> Vec<Option<bool>> {
        indices
            .iter()
            .map(|&idx| self.measure_one_with_scratch(idx, scratch))
            .collect()
    }

    /// One qubit, reusing `scratch`. A lost qubit pushes `None` and returns
    /// `None`, exactly as the standalone `measure` path does — which is what
    /// makes every batched path observationally identical to a per-qubit loop.
    fn measure_one_with_scratch(
        &mut self,
        idx: usize,
        scratch: &mut MeasureScratch<I>,
    ) -> Option<bool> {
        if self.is_lost[idx] {
            self.measurement_record.push(None);
            return None;
        }
        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            self.compute_decomposition(idx, Pauli::Z);
        self.measure_with_scratch(
            idx,
            scratch,
            phase_decomp,
            stab_anticomm_bits,
            destab_anticomm_bits,
        )
    }

    /// The coefficient-aware measurement kernel.
    ///
    /// The `stab_anticomm_bits == 0` dichotomy is an **explicit** branch, not a
    /// special case of the general path:
    ///
    /// * **Case b** (`Z` is already a stabilizer): the Z-overlap is computed in
    ///   place over `coefficients` — no drain, no map, no sort — the projection
    ///   is an in-place `retain` on the parity predicate with **no magnitude
    ///   filter at all**, and `normalize()` runs *only* if the support shrank.
    ///   The frame is untouched, so the cached odd-phase mask stays valid.
    /// * **Case a**: sort-merge overlap, then partition/transform/merge, then a
    ///   cutoff that is **relative to the merged norm**
    ///   (`|c|² > threshold² · ‖v‖²`), then an unconditional `normalize()`, then
    ///   the frame projection (which invalidates the cached mask).
    ///
    /// The case-b arithmetic is machine-checked in
    /// `lean/PPVM/Tableau/Projection.lean` as the `s = 0` specialization of the
    /// Born projector: `selfInverse_zero_phase_even` forces every per-index
    /// phase into `{0, 2}` (the `debug_assert!` below, and the
    /// `!phase.is_multiple_of(2)` guard in the overlap loop is then vacuous);
    /// `proj_zero_apply` shows the projector keeps survivors with factor **1**
    /// — not the factor `2` of case a (`projectRaw_eq_two_proj`) — which is why
    /// case b applies no magnitude filter; `proj_zero_eq_self` shows a case-b
    /// projection that drops nothing is the identity, hence exactly
    /// norm-preserving, which is why `normalize()` runs only when the support
    /// shrank; and `proj_zero_eq_caseB_retain` shows the `retain` predicate
    /// `(parity ⊕ outcome) == (phase_decomp == 2)` is exactly the indicator of
    /// the surviving set.
    pub(crate) fn measure_with_scratch(
        &mut self,
        addr0: usize,
        scratch: &mut MeasureScratch<I>,
        phase_decomp: u8,
        stab_anticomm_bits: I,
        destab_anticomm_bits: I,
    ) -> Option<bool> {
        if stab_anticomm_bits == I::zero() {
            // ── Case b (fast path) ───────────────────────────────────────
            let mut z_overlap_re = 0.0f64;
            for &(coeff, idx) in self.coefficients.iter() {
                let phase = (phase_decomp
                    + (2 * symplectic_inner(destab_anticomm_bits, idx) as u8) % 4)
                    % 4;
                if !phase.is_multiple_of(2) {
                    continue;
                }
                let norm_sq = coeff.norm_sqr();
                if phase == 0 {
                    z_overlap_re += norm_sq;
                } else {
                    z_overlap_re -= norm_sq;
                }
            }

            let prob_1 = 0.5 - 0.5 * z_overlap_re;
            let outcome = self.tableau.rng.random::<f64>() < prob_1;

            // A theorem, not a hope: `lean/PPVM/Tableau/BranchPhase.lean`'s
            // `frameInvolution_zero_iff` shows that with `stab_anticomm == 0`
            // the frame identity (`M² = I`) collapses to exactly
            // `phase_decomp ∈ {0, 2}`. The assert stays as a cheap witness.
            debug_assert!(
                phase_decomp == 0 || phase_decomp == 2,
                "Measurement result cannot be imaginary!"
            );

            let old_len = self.coefficients.len();
            let z_sign = phase_decomp == 2;
            self.coefficients.retain_entries(|(_, alpha)| {
                let parity = symplectic_inner(*alpha, destab_anticomm_bits) % 2 != 0;
                (parity ^ outcome) == z_sign
            });
            if self.coefficients.len() < old_len {
                self.coefficients.normalize();
            }

            // Case b doesn't mutate destabilizers, so the cached mask stays valid.
            self.measurement_record.push(Some(outcome));
            Some(outcome)
        } else {
            // ── Case a ───────────────────────────────────────────────────
            let mut by_idx = std::mem::take(&mut scratch.by_idx);
            by_idx.clear();
            by_idx.extend(self.coefficients.take().into_iter().map(|(c, i)| (i, c)));
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

            // OVERLAP: 2-way merge of `by_idx` and `shifted` (each key XOR'd by
            // `stab_anticomm_bits`). At an equal key `k`, `by_idx` has `(k, a_k)`
            // and `shifted` has `(k, a_{k^s})` — the same `coeff / coeff_branch`
            // pair a map-based overlap would visit, counted once per key.
            let mut shifted = std::mem::take(&mut scratch.shifted);
            shifted.clear();
            shifted.extend(by_idx.iter().map(|&(i, c)| (i ^ stab_anticomm_bits, c)));
            shifted.sort_unstable_by_key(|a| a.0);

            let mut z_overlap_re = 0.0f64;
            {
                let mut ii = 0usize;
                let mut jj = 0usize;
                while ii < by_idx.len() && jj < shifted.len() {
                    match by_idx[ii].0.cmp(&shifted[jj].0) {
                        std::cmp::Ordering::Less => ii += 1,
                        std::cmp::Ordering::Greater => jj += 1,
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
                            let re_w = a.re * b.re + a.im * b.im;
                            let im_w = a.re * b.im - a.im * b.re;
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

            // PROJECTION: partition A (k-bit = 0) and B (k-bit = 1), transform
            // B, merge.
            let q_idx = stab_anticomm_bits.trailing_zeros() as usize;
            let k = I::one() << q_idx;
            let alpha = if outcome {
                (phase_decomp + 2) % 4
            } else {
                phase_decomp
            };

            let mut a = std::mem::take(&mut scratch.a);
            let mut bt = std::mem::take(&mut scratch.bt);
            a.clear();
            bt.clear();
            for (idx, coeff) in by_idx.drain(..) {
                if (idx & k) == I::zero() {
                    a.push((idx, coeff));
                } else {
                    let symp = symplectic_inner(idx, destab_anticomm_bits);
                    let phase_idx =
                        ((alpha as i32 + if symp % 2 == 1 { 2 } else { 0 }) % 4) as usize;
                    let q = COMPLEX_PHASE_CONVERSION[phase_idx];
                    bt.push((idx ^ stab_anticomm_bits, q * coeff));
                }
            }
            // `a` is already sorted (a subset of the sorted `by_idx`); `bt` is not.
            bt.sort_unstable_by_key(|e| e.0);

            // 2-way merge summing equal keys → sorted merged output.
            let mut merged = std::mem::take(&mut scratch.merged);
            merged.clear();
            merged.reserve(a.len() + bt.len());
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

            // The case-a cutoff is RELATIVE to the merged norm — unlike the
            // gates' absolute rule and unlike case b, which has no filter.
            let norm_sqr = merged.iter().fold(0.0f64, |acc, (_, c)| acc + c.norm_sqr());
            let cutoff_sq = self.coefficient_threshold * self.coefficient_threshold;
            let threshold = cutoff_sq * norm_sqr;
            self.coefficients.reserve(merged.len());
            for &(idx, coeff) in merged.iter() {
                if coeff.norm_sqr() > threshold {
                    self.coefficients.unsafe_insert(idx, coeff);
                }
            }

            self.coefficients.normalize();
            self.tableau
                .update_tableau_according_to_outcome(addr0, q_idx, outcome);
            // Destabilizer phases just changed; invalidate the cached mask.
            scratch.odd_phase_mask = None;
            self.measurement_record.push(Some(outcome));

            // Hand the (now-cleared, still-allocated) buffers back.
            merged.clear();
            scratch.by_idx = by_idx;
            scratch.shifted = shifted;
            scratch.a = a;
            scratch.bt = bt;
            scratch.merged = merged;

            Some(outcome)
        }
    }

    /// Project the state onto a **given** case-a outcome using the scratch's
    /// map working set. Retained for callers (the stim executor) that sample the
    /// outcome themselves.
    pub fn project_case_a(
        &mut self,
        outcome: bool,
        scratch: &mut MeasureScratch<I>,
        phase_decomp: u8,
        stab_anticomm_bits: I,
        destab_anticomm_bits: I,
        addr0: usize,
    ) {
        let q_idx = stab_anticomm_bits.trailing_zeros() as usize;

        let one = I::one();
        let zero = I::zero();
        let k = one << q_idx;

        let alpha = if outcome {
            (phase_decomp + 2) % 4
        } else {
            phase_decomp
        };

        // Partition into A (k-bit = 0) and B (k-bit = 1) via `retain`, then
        // merge. Split the borrow so `retain` can mutate `coeff_map` while the
        // closure pushes into `b_entries`.
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
                false // remove the B entry
            } else {
                true // keep the A entry
            }
        });
        Self::merge_b_into_a(
            coeff_map,
            b_entries,
            alpha,
            destab_anticomm_bits,
            stab_anticomm_bits,
        );

        let norm_sqr = coeff_map.values().fold(0.0f64, |acc, c| acc + c.norm_sqr());

        let cutoff_sq = self.coefficient_threshold * self.coefficient_threshold;
        let threshold = cutoff_sq * norm_sqr;
        self.coefficients.reserve(coeff_map.len());
        for (idx, coeff) in coeff_map.drain() {
            if coeff.norm_sqr() > threshold {
                self.coefficients.unsafe_insert(idx, coeff);
            }
        }

        self.coefficients.normalize();

        self.tableau
            .update_tableau_according_to_outcome(addr0, q_idx, outcome);
        scratch.odd_phase_mask = None;
    }

    /// Project the state onto a given case-b outcome (`Z` is a stabilizer).
    ///
    /// No magnitude filter: only the parity predicate. `normalize()` runs only
    /// if the support shrank.
    pub fn project_case_b(
        &mut self,
        entries: &[(Complex64, I)],
        outcome: bool,
        phase_decomp: u8,
        destab_anticomm_bits: I,
    ) {
        let old_len = entries.len();
        let z_sign = phase_decomp == 2;

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

    /// Case-b overlap: self-pairing (`branch_index == idx`), so the overlap is
    /// `±|c|²` and only even phases contribute to the real part.
    pub fn compute_overlap_case_b(
        entries: &[(Complex64, I)],
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
            let norm_sq = coeff.norm_sqr();
            if phase == 0 {
                z_overlap_re += norm_sq;
            } else {
                z_overlap_re -= norm_sq;
            }
        }
        z_overlap_re
    }

    /// Case-a overlap: cross-index pairing via map lookup; accumulates only the
    /// real part.
    pub fn compute_overlap_case_a(
        coeff_map: &HashMap<I, Complex64>,
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
            let re_w = coeff.re * coeff_branch.re + coeff.im * coeff_branch.im;
            let im_w = coeff.re * coeff_branch.im - coeff.im * coeff_branch.re;
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

    /// Merge the B entries (`k`-bit = 1) into their A partners with the phase
    /// adjustment.
    fn merge_b_into_a(
        coeff_map: &mut HashMap<I, Complex64>,
        b_entries: &[(I, Complex64)],
        alpha: u8,
        destab_anticomm_bits: I,
        stab_anticomm_bits: I,
    ) {
        for &(idx, coeff) in b_entries {
            let symp_inner = symplectic_inner(idx, destab_anticomm_bits);
            let phase_idx = ((alpha as i32 + if symp_inner % 2 == 1 { 2 } else { 0 }) % 4) as usize;
            let q = COMPLEX_PHASE_CONVERSION[phase_idx];
            *coeff_map
                .entry(idx ^ stab_anticomm_bits)
                .or_insert(Complex64::new(0.0, 0.0)) += q * coeff;
        }
    }

    /// Measure `qubit` in the Z basis with readout noise.
    ///
    /// Behaves like [`Measure::measure`], then flips the *recorded* bit with
    /// probability `flip_prob`. The quantum state stays consistent with the
    /// **true** outcome — only the returned/recorded value flips — and exactly
    /// one record entry is pushed per logical measurement (the internal
    /// `measure` pushes it, this overwrites it). Returns `None` for a lost
    /// qubit, regardless of `flip_prob`.
    pub fn measure_noisy(&mut self, qubit: usize, flip_prob: f64) -> Option<bool> {
        debug_assert!(
            (0.0..=1.0).contains(&flip_prob),
            "flip_prob must be in [0, 1], got {flip_prob}"
        );
        let outcome = self.measure(qubit)?;
        let noisy = self.flip_with_prob(outcome, flip_prob);
        self.overwrite_last_measurement_record(Some(noisy));
        Some(noisy)
    }

    /// Sample a `Bernoulli(p)` outcome from the tableau's internal RNG.
    pub fn bernoulli(&mut self, p: f64) -> bool {
        debug_assert!((0.0..=1.0).contains(&p), "p must be in [0, 1], got {p}");
        self.tableau.rng.random::<f64>() < p
    }

    /// Flip `bit` with probability `p`. Returns `bit` unchanged — and draws
    /// **no** randomness — when `p <= 0.0`.
    pub fn flip_with_prob(&mut self, bit: bool, p: f64) -> bool {
        debug_assert!((0.0..=1.0).contains(&p), "p must be in [0, 1], got {p}");
        if p > 0.0 && self.bernoulli(p) {
            !bit
        } else {
            bit
        }
    }
}

// ─── Reset ────────────────────────────────────────────────────────────────

impl<A: RowStorage, H> Reset for Tableau<A, H> {
    /// Measure, then `X` if the outcome was `1`. No record exists on a bare
    /// frame.
    fn reset(&mut self, qubit: usize) {
        if let Some(true) = Measure::measure(self, qubit) {
            Clifford::x(self, qubit);
        }
    }
}

impl<A: RowStorage, I: Bitstring, H> Reset for GeneralizedTableau<A, I, H> {
    /// Measure, drop the record entry the internal `measure` pushed, then `X` if
    /// the outcome was `1`.
    ///
    /// A reset is not a measurement in stim's model, so the operation is
    /// **measurement-record-neutral**.
    fn reset(&mut self, qubit: usize) {
        let m = Measure::measure(self, qubit);
        self.measurement_record.pop();
        if let Some(true) = m {
            Clifford::x(self, qubit);
        }
    }
}
