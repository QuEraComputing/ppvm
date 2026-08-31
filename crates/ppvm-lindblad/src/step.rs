// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Predictor-corrector adaptive step `O ← exp(dt·L*) O`.

use crate::spec::LindbladSpec;
use crate::word::Word;
use crate::{Error, PcStepConfig, mf_expm};
use fxhash::FxHashSet;
use std::time::Instant;

/// Per-phase timing breakdown (microseconds) returned by
/// [`LindbladSpec::pc_step_timed`].
#[derive(Default, Clone, Copy, Debug)]
pub struct PcStepTimings {
    pub leakage1_us: u64,
    pub expand1_us: u64,
    pub expm1_us: u64,
    pub leakage2_us: u64,
    pub expand2_us: u64,
    pub expm2_us: u64,
}

impl PcStepTimings {
    pub fn total_us(&self) -> u64 {
        self.leakage1_us
            + self.expand1_us
            + self.expm1_us
            + self.leakage2_us
            + self.expand2_us
            + self.expm2_us
    }
}

/// Compact `basis` / `coeffs` in place: drop entries whose absolute
/// coefficient is below `drop_tol` unless the word appears in `protected`.
/// No-op when `drop_tol ≤ 0`.
fn prune_basis(basis: &mut Vec<Word>, coeffs: &mut Vec<f64>, drop_tol: f64, protected: &[Word]) {
    if drop_tol <= 0.0 {
        return;
    }
    debug_assert_eq!(basis.len(), coeffs.len());
    let protected_set: FxHashSet<&Word> = protected.iter().collect();
    let mut write = 0;
    for read in 0..basis.len() {
        if coeffs[read].abs() >= drop_tol || protected_set.contains(&basis[read]) {
            if write != read {
                basis.swap(write, read);
                coeffs.swap(write, read);
            }
            write += 1;
        }
    }
    basis.truncate(write);
    coeffs.truncate(write);
}

/// Global max-basis cap (PauliStrings.jl-style top-M trim): keep only the
/// `max_basis` largest-|coeff| terms (protected strings always kept),
/// dropping the rest. Rank-based total-basis bound; dual of `drop_tol`.
/// A `max_basis` large enough to cover the whole basis is a no-op.
fn cap_basis(basis: &mut Vec<Word>, coeffs: &mut Vec<f64>, max_basis: usize, protected: &[Word]) {
    if basis.len() <= max_basis {
        return;
    }
    let protected_set: FxHashSet<&Word> = protected.iter().collect();
    let n_prot = basis.iter().filter(|w| protected_set.contains(w)).count();
    let slots = max_basis.saturating_sub(n_prot);
    let mut mags: Vec<f64> = basis
        .iter()
        .zip(coeffs.iter())
        .filter(|(w, _)| !protected_set.contains(w))
        .map(|(_, c)| c.abs())
        .collect();
    let cutoff = if slots == 0 {
        f64::INFINITY
    } else if slots >= mags.len() {
        return;
    } else {
        let k = slots - 1;
        mags.select_nth_unstable_by(k, |a, b| {
            b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal)
        });
        mags[k]
    };
    let mut write = 0;
    for read in 0..basis.len() {
        if protected_set.contains(&basis[read]) || coeffs[read].abs() >= cutoff {
            if write != read {
                basis.swap(write, read);
                coeffs.swap(write, read);
            }
            write += 1;
        }
    }
    basis.truncate(write);
    coeffs.truncate(write);
}

/// Add the largest leakage strings to the basis, up to the available room
/// `room = max_basis − basis.len()` — so the in-step basis (hence the
/// expm/leakage peak memory) never exceeds `max_basis`. New strings get
/// coefficient 0; the surrounding expm fills them. No magnitude filter: the
/// top-`room` by `|leakage|` are added (a large `max_basis` adds them all).
fn add_leakage_capped(
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<f64>,
    mut leak: Vec<(Word, f64)>,
    max_basis: usize,
) {
    let room = max_basis.saturating_sub(basis.len());
    if leak.len() > room {
        if room > 0 {
            leak.select_nth_unstable_by(room - 1, |a, b| {
                b.1.abs()
                    .partial_cmp(&a.1.abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }
        leak.truncate(room);
    }
    for (w, _) in leak {
        basis.push(w);
        coeffs.push(0.0);
    }
}

impl LindbladSpec {
    /// One predictor-corrector step `O ← exp(dt·L*) O` in the adaptive
    /// real-coefficient Pauli basis: first-hop leakage admission, predictor
    /// exponential, second-hop admission from the predicted state, corrector
    /// exponential from the saved pre-step state, then truncation (prune +
    /// rank cap) per [`PcStepConfig`]. Exact in `dt` within the working
    /// basis — the only error is basis truncation.
    ///
    /// `protected` words are never dropped. All tuning knobs live in `cfg`.
    pub fn pc_step(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word],
        cfg: &PcStepConfig,
    ) -> Result<(), Error> {
        self.run_in_pool(cfg, |this| {
            this.pc_step_inner(basis, coeffs, dt, protected, cfg)
                .map(|_| ())
        })
    }

    /// Same as [`Self::pc_step`] but also returns a per-phase timing
    /// breakdown (microseconds), for profiling parallel scaling and hot
    /// spots.
    pub fn pc_step_timed(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word],
        cfg: &PcStepConfig,
    ) -> Result<PcStepTimings, Error> {
        self.run_in_pool(cfg, |this| {
            this.pc_step_inner(basis, coeffs, dt, protected, cfg)
        })
    }

    fn run_in_pool<R: Send>(
        &self,
        cfg: &PcStepConfig,
        f: impl FnOnce(&Self) -> Result<R, Error> + Send,
    ) -> Result<R, Error> {
        if let Some(n) = cfg.num_threads {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| Error::Internal(format!("rayon pool build: {e}")))?;
            pool.install(|| f(self))
        } else {
            f(self)
        }
    }

    fn pc_step_inner(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word],
        cfg: &PcStepConfig,
    ) -> Result<PcStepTimings, Error> {
        let PcStepConfig {
            max_basis,
            admit_basis,
            drop_tol,
            tau_add,
            ..
        } = *cfg;
        // Admission bound: enrichment may grow the live basis to `admit`
        // >= `max_basis`; the final `cap_basis` then keeps the top-
        // `max_basis` strings by evolved |coeff| over the whole union
        // (retained + admitted) — rank displacement. With `admit_basis =
        // None` admission is bounded by `max_basis` itself, `cap_basis` is
        // a no-op, and membership turnover requires `drop_tol > 0`.
        let admit = admit_basis.unwrap_or(max_basis).max(max_basis);
        let tau_add = tau_add.unwrap_or(0.0);
        let mut t = PcStepTimings::default();

        // 1. First-hop expansion. After this, `coeffs` contains the pre-step
        // coefficients followed by zeros for the newly-added leakage strings.
        // We rely on `coeffs` itself as the pre-step buffer for the corrector
        // — no `.clone()` is needed because `expm_step` only borrows it.
        let t0 = Instant::now();
        let leak = self.leakage_with_prune(basis, coeffs, protected, admit, tau_add)?;
        t.leakage1_us = t0.elapsed().as_micros() as u64;

        let t0 = Instant::now();
        add_leakage_capped(basis, coeffs, leak, admit);
        t.expand1_us = t0.elapsed().as_micros() as u64;

        // 2. Predictor: `expm_step` reads `coeffs` immutably and returns a
        // new owned vector with the predicted state.
        let t0 = Instant::now();
        let coeffs_predict = self.expm_step(basis, dt, coeffs, drop_tol);
        t.expm1_us = t0.elapsed().as_micros() as u64;

        // 3. Second-hop expansion from the predicted state. After leakage2
        // we no longer need `coeffs_predict`. Extend `coeffs` with zeros for
        // any newly-added second-hop strings so it remains a valid input
        // (pre-step state) for the corrector.
        let t0 = Instant::now();
        let leak2 = self.leakage_with_prune(basis, &coeffs_predict, protected, admit, tau_add)?;
        t.leakage2_us = t0.elapsed().as_micros() as u64;
        drop(coeffs_predict);

        let t0 = Instant::now();
        add_leakage_capped(basis, coeffs, leak2, admit);
        t.expand2_us = t0.elapsed().as_micros() as u64;

        // 4. Corrector: redo from pre-step state on the doubly-enlarged basis.
        let t0 = Instant::now();
        *coeffs = self.expm_step(basis, dt, coeffs, drop_tol);
        t.expm2_us = t0.elapsed().as_micros() as u64;

        // 5. Prune basis entries below `drop_tol` (protected words never dropped).
        prune_basis(basis, coeffs, drop_tol, protected);
        cap_basis(basis, coeffs, max_basis, protected);
        Ok(t)
    }

    /// Compute `exp(dt · M) · b` for the in-basis-restricted generator
    /// `M`, matrix-free, via `quspin-expm` (see [`crate::mf_expm`]).
    fn expm_step(&self, basis: &[Word], dt: f64, b: &[f64], drop_tol: f64) -> Vec<f64> {
        mf_expm::expm_apply_mf(self, basis, dt, b, drop_tol)
    }
}
