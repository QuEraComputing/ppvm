// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Predictor-corrector adaptive step `O ← exp(dt·L*) O`.

use crate::sector::Sector;
use crate::spec::LindbladSpec;
use crate::truncate::{add_leakage_capped, cap_basis, prune_basis};
use crate::word::Word;
use crate::{Error, PcStepConfig, mf_expm};
use num::Complex;
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

/// Clock for one `pc_step` phase. Disarmed (`None`) on the untimed path, so
/// [`LindbladSpec::pc_step`] pays no `Instant` syscalls.
struct Phase(Option<Instant>);

impl Phase {
    fn start(timed: bool) -> Self {
        Self(timed.then(Instant::now))
    }

    fn stop(self, slot: &mut u64) {
        if let Some(t0) = self.0 {
            *slot = t0.elapsed().as_micros() as u64;
        }
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
            this.pc_step_inner(basis, coeffs, dt, protected, cfg, false)
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
            this.pc_step_inner(basis, coeffs, dt, protected, cfg, true)
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
        timed: bool,
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
        let p = Phase::start(timed);
        let leak = self.leakage_with_prune(basis, coeffs, protected, admit, tau_add)?;
        p.stop(&mut t.leakage1_us);

        let p = Phase::start(timed);
        add_leakage_capped(basis, coeffs, leak, admit);
        p.stop(&mut t.expand1_us);

        // 2. Predictor: `expm_step` reads `coeffs` immutably and returns a
        // new owned vector with the predicted state.
        let p = Phase::start(timed);
        let coeffs_predict = self.expm_step(basis, dt, coeffs, drop_tol);
        p.stop(&mut t.expm1_us);

        // 3. Second-hop expansion from the predicted state. After leakage2
        // we no longer need `coeffs_predict`. Extend `coeffs` with zeros for
        // any newly-added second-hop strings so it remains a valid input
        // (pre-step state) for the corrector.
        let p = Phase::start(timed);
        let leak2 = self.leakage_with_prune(basis, &coeffs_predict, protected, admit, tau_add)?;
        p.stop(&mut t.leakage2_us);
        drop(coeffs_predict);

        let p = Phase::start(timed);
        add_leakage_capped(basis, coeffs, leak2, admit);
        p.stop(&mut t.expand2_us);

        // 4. Corrector: redo from pre-step state on the doubly-enlarged basis.
        let p = Phase::start(timed);
        *coeffs = self.expm_step(basis, dt, coeffs, drop_tol);
        p.stop(&mut t.expm2_us);

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

    /// One predictor-corrector step in **orbit-rep form** at momentum
    /// `sector`: the same five phases as [`Self::pc_step`], but the basis
    /// holds only canonical translation-orbit representatives, the
    /// coefficients are complex, and the `L*` action is phase-aware (see
    /// [`crate::sector`]). The basis stays ~`|G|`× smaller than the
    /// equivalent full-basis complex evolution, every step.
    ///
    /// `max_basis` is a hard rank cap on the live orbit-rep basis:
    /// enrichment adds at most `admit − basis.len()` of the largest
    /// leakage reps, the leakage map is capped to the same room, and the
    /// post-step basis is trimmed to the top-`max_basis` reps by `|c|`.
    /// Pass a large value (e.g. `usize::MAX`) for the near-exact,
    /// uncapped case. `drop_tol` additionally prunes by magnitude.
    /// `protected` reps are never dropped.
    ///
    /// `basis` is assumed to contain only canonical orbit
    /// representatives. If not, call
    /// [`canonicalize_basis_to_rep`](crate::canonicalize_basis_to_rep)
    /// first.
    ///
    /// Honours `cfg.num_threads` the same way [`Self::pc_step`] does.
    pub fn pc_step_orbit_rep(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<Complex<f64>>,
        dt: f64,
        protected: &[Word],
        sector: Sector<'_>,
        cfg: &PcStepConfig,
    ) -> Result<(), Error> {
        self.run_in_pool(cfg, |this| {
            this.pc_step_orbit_rep_inner(basis, coeffs, dt, protected, sector, cfg)
        })
    }

    fn pc_step_orbit_rep_inner(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<Complex<f64>>,
        dt: f64,
        protected: &[Word],
        sector: Sector<'_>,
        cfg: &PcStepConfig,
    ) -> Result<(), Error> {
        let PcStepConfig {
            max_basis,
            admit_basis,
            drop_tol,
            tau_add,
            ..
        } = *cfg;
        // Admission bound, mirroring `pc_step_inner`: enrichment may grow
        // the live basis to `admit` >= `max_basis`; the final `cap_basis`
        // keeps the top-`max_basis` reps by evolved |coeff| over the whole
        // union (rank displacement). With `admit_basis = None` admission is
        // bounded by `max_basis` itself and membership turnover requires
        // `drop_tol > 0`.
        let admit = admit_basis.unwrap_or(max_basis).max(max_basis);
        let tau_add = tau_add.unwrap_or(0.0);

        // 1. First-hop phase-aware leakage.
        let mut leak = self.leakage_orbit_rep(basis, coeffs, protected, sector, admit)?;
        if tau_add > 0.0 {
            leak.retain(|(_, c)| c.norm() > tau_add);
        }
        add_leakage_capped(basis, coeffs, leak, admit);

        // 2. Predictor: the phase-aware action is built once and reused
        //    across every matvec.
        let coeffs_predict = mf_expm::expm_apply_orbit_rep(self, basis, sector, dt, coeffs);

        // 3. Second-hop leakage from the predicted state.
        let mut leak2 = self.leakage_orbit_rep(basis, &coeffs_predict, protected, sector, admit)?;
        drop(coeffs_predict);
        if tau_add > 0.0 {
            leak2.retain(|(_, c)| c.norm() > tau_add);
        }
        add_leakage_capped(basis, coeffs, leak2, admit);

        // 4. Corrector: redo from the pre-step state (the basis grew).
        *coeffs = mf_expm::expm_apply_orbit_rep(self, basis, sector, dt, coeffs);

        // 5. Prune by magnitude, then rank-cap to max_basis.
        prune_basis(basis, coeffs, drop_tol, protected);
        cap_basis(basis, coeffs, max_basis, protected);
        Ok(())
    }
}
