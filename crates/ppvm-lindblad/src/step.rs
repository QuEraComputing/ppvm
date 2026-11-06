// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Predictor-corrector adaptive step `O ← exp(dt·L*) O`.

use crate::sector::Sector;
use crate::spec::LindbladSpec;
use crate::truncate::{
    add_leakage_capped, cap_basis, cap_basis_by_score, cap_basis_with_scratch, prune_basis,
};
use crate::word::Word;
use crate::{
    CandidateBuffer, CandidateProposer, CandidateStrategy, Error, PcStepConfig, PcStepWorkspace,
    mf_expm,
};
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
pub(crate) struct Phase(Option<Instant>);

impl Phase {
    pub(crate) fn start(timed: bool) -> Self {
        Self(timed.then(Instant::now))
    }

    pub(crate) fn stop(self, slot: &mut u64) {
        if let Some(t0) = self.0 {
            *slot = t0.elapsed().as_micros() as u64;
        }
    }
}

impl<const C: usize> LindbladSpec<C> {
    /// One predictor-corrector step `O ← exp(dt·L*) O` in the adaptive
    /// real-coefficient Pauli basis: first-hop leakage admission, predictor
    /// exponential, second-hop admission from the predicted state, corrector
    /// exponential from the saved pre-step state, then truncation (prune +
    /// rank cap) per [`PcStepConfig`]. Exact in `dt` within the working
    /// basis — the only error is basis truncation.
    ///
    /// When the second hop admits no string (in particular when the first
    /// hop already filled the admission room `admit_basis − |basis|`, the
    /// usual case once the basis has reached `max_basis`), the corrector
    /// would repeat the predictor exactly; the second leakage pass (if
    /// `room = 0`) and the corrector exponential are then skipped, with
    /// bit-identical results. The step's timings report 0 for skipped
    /// phases.
    ///
    /// `protected` words are never dropped. All tuning knobs live in `cfg`.
    ///
    /// The generator is fixed during this call. Piecewise-constant driving
    /// requires a spec for each segment and steps split at every switch.
    /// Forward-looking retention extrapolates only the current step's change;
    /// it does not anticipate subsequent drive or dissipation changes.
    pub fn pc_step(
        &self,
        basis: &mut Vec<Word<C>>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word<C>],
        cfg: &PcStepConfig,
    ) -> Result<(), Error> {
        self.pc_step_pooled(basis, coeffs, dt, protected, cfg, false)
            .map(|_| ())
    }

    /// Same as [`Self::pc_step`] but also returns a per-phase timing
    /// breakdown (microseconds), for profiling parallel scaling and hot
    /// spots.
    pub fn pc_step_timed(
        &self,
        basis: &mut Vec<Word<C>>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word<C>],
        cfg: &PcStepConfig,
    ) -> Result<PcStepTimings, Error> {
        self.pc_step_pooled(basis, coeffs, dt, protected, cfg, true)
    }

    fn pc_step_pooled(
        &self,
        basis: &mut Vec<Word<C>>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word<C>],
        cfg: &PcStepConfig,
        timed: bool,
    ) -> Result<PcStepTimings, Error> {
        // Only checkout and return hold the lock; concurrent trajectories use separate buffers.
        let mut workspace = self
            .workspaces
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .pop()
            .unwrap_or_default();
        let result = self.pc_step_execute(basis, coeffs, dt, protected, cfg, timed, &mut workspace);
        self.workspaces
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .push(workspace);
        result
    }

    /// Adaptive step using buffers owned explicitly by the caller.
    ///
    /// Ordinary step calls reuse an internal workspace automatically. This entry point
    /// lets a caller retain the buffers independently of the model.
    /// Capacity follows the active state; a rank limit does not allocate its full size.
    /// Numerical ordering is identical to [`Self::pc_step_timed`].
    pub fn pc_step_with_workspace(
        &self,
        basis: &mut Vec<Word<C>>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word<C>],
        cfg: &PcStepConfig,
        workspace: &mut PcStepWorkspace<C>,
    ) -> Result<PcStepTimings, Error> {
        self.pc_step_execute(basis, coeffs, dt, protected, cfg, true, workspace)
    }

    #[allow(clippy::too_many_arguments)]
    fn pc_step_execute(
        &self,
        basis: &mut Vec<Word<C>>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word<C>],
        cfg: &PcStepConfig,
        timed: bool,
        workspace: &mut PcStepWorkspace<C>,
    ) -> Result<PcStepTimings, Error> {
        let mut proposer = std::mem::take(&mut workspace.proposer);
        let result = match cfg.proposer {
            CandidateStrategy::Streaming => self.pc_step_execute_with_proposer(
                basis,
                coeffs,
                dt,
                protected,
                cfg,
                timed,
                workspace,
                &mut proposer.streaming,
            ),
            CandidateStrategy::Shortlist => self.pc_step_execute_with_proposer(
                basis,
                coeffs,
                dt,
                protected,
                cfg,
                timed,
                workspace,
                &mut proposer,
            ),
        };
        workspace.proposer = proposer;
        result
    }

    /// Adaptive step with a caller-owned candidate strategy and reusable workspace.
    ///
    /// The strategy receives a borrowed view of the active sum and its spare word
    /// storage. The stepper initializes admitted coefficients, evolves the state,
    /// and applies final truncation. Ordinary step calls use the strategy in [`PcStepConfig::proposer`].
    #[allow(clippy::too_many_arguments)]
    pub fn pc_step_with_proposer<P: CandidateProposer<C>>(
        &self,
        basis: &mut Vec<Word<C>>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word<C>],
        cfg: &PcStepConfig,
        workspace: &mut PcStepWorkspace<C>,
        proposer: &mut P,
    ) -> Result<PcStepTimings, Error> {
        self.pc_step_execute_with_proposer(
            basis, coeffs, dt, protected, cfg, true, workspace, proposer,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn pc_step_execute_with_proposer<P: CandidateProposer<C>>(
        &self,
        basis: &mut Vec<Word<C>>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word<C>],
        cfg: &PcStepConfig,
        timed: bool,
        workspace: &mut PcStepWorkspace<C>,
        proposer: &mut P,
    ) -> Result<PcStepTimings, Error> {
        let ranking_factor = cfg.ranking_factor(dt)?;
        if basis.len() != coeffs.len() {
            return Err(Error::LengthMismatch {
                what: "basis and coeffs",
                a: basis.len(),
                b: coeffs.len(),
            });
        }
        let dimension = 1usize
            .checked_shl((2 * self.n_qubits()) as u32)
            .unwrap_or(usize::MAX);
        let admit = cfg.admit_basis.unwrap_or(cfg.max_basis).max(cfg.max_basis);
        {
            // Reserve from the active state rather than an arbitrarily large rank limit.
            let active = basis.len().min(cfg.max_basis);
            let bound = active.saturating_mul(3).min(admit).min(dimension);
            if let Some(required) = bound.checked_mul(2)
                && workspace.taylor.capacity() < required
            {
                workspace
                    .taylor
                    .reserve_exact(required - workspace.taylor.len());
            }
            if let Some(reserve) = active.checked_mul(4) {
                let reserve = reserve.max(bound).min(dimension);
                if basis.capacity() < reserve {
                    basis.reserve_exact(reserve.saturating_sub(basis.len()));
                }
                if coeffs.capacity() < reserve {
                    coeffs.reserve_exact(reserve.saturating_sub(coeffs.len()));
                }
            }
        }
        self.run_in_pool(cfg, |this| {
            this.pc_step_inner(
                basis,
                coeffs,
                dt,
                protected,
                cfg,
                timed,
                workspace,
                proposer,
                ranking_factor,
            )
        })
    }

    pub(crate) fn run_in_pool<R: Send>(
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

    #[allow(clippy::too_many_arguments)]
    fn pc_step_inner<P: CandidateProposer<C>>(
        &self,
        basis: &mut Vec<Word<C>>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word<C>],
        cfg: &PcStepConfig,
        timed: bool,
        workspace: &mut PcStepWorkspace<C>,
        proposer: &mut P,
        ranking_factor: Option<f64>,
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

        let original_terms = basis.len();

        // No basis permutation occurs before final truncation, so a saved
        // coefficient prefix retains the original word positions for correction.
        // 1. First-hop expansion appends zero coefficients for admitted words.
        let p = Phase::start(timed);
        let mut buffer = CandidateBuffer::new(basis, coeffs, &mut workspace.taylor);
        proposer.propose(self, &mut buffer, protected, admit)?;
        p.stop(&mut t.leakage1_us);

        let p = Phase::start(timed);
        let count = proposer.select(&mut buffer, tau_add);
        assert!(count <= admit.saturating_sub(buffer.basis().len()));
        buffer.finish(count);
        coeffs.resize(basis.len(), 0.0);
        p.stop(&mut t.expand1_us);

        if basis.len() < admit || ranking_factor.is_some() {
            workspace.coefficient_auxiliary.clear();
            workspace
                .coefficient_auxiliary
                .extend_from_slice(&coeffs[..original_terms]);
        }

        // The main coefficient vector is the exponential's output; the prefix backup is only for a possible corrector.
        let p = Phase::start(timed);
        mf_expm::expm_apply_inplace(self, basis, dt, coeffs, drop_tol, workspace);
        p.stop(&mut t.expm1_us);

        let n_predict = basis.len();
        if admit > n_predict {
            let p = Phase::start(timed);
            let mut buffer = CandidateBuffer::new(basis, coeffs, &mut workspace.taylor);
            proposer.propose(self, &mut buffer, protected, admit)?;
            p.stop(&mut t.leakage2_us);
            let p = Phase::start(timed);
            let count = proposer.select(&mut buffer, tau_add);
            assert!(count <= admit.saturating_sub(buffer.basis().len()));
            buffer.finish(count);
            coeffs.resize(basis.len(), 0.0);
            p.stop(&mut t.expand2_us);
        }

        if basis.len() != n_predict {
            let p = Phase::start(timed);
            coeffs.clear();
            coeffs.extend_from_slice(&workspace.coefficient_auxiliary);
            coeffs.resize(basis.len(), 0.0);
            mf_expm::expm_apply_inplace(self, basis, dt, coeffs, drop_tol, workspace);
            p.stop(&mut t.expm2_us);
        }

        if let Some(factor) = ranking_factor {
            // Score before compaction so the saved coefficient prefix stays aligned.
            let history = &workspace.coefficient_auxiliary;
            cap_basis_by_score(
                basis,
                coeffs,
                max_basis,
                protected,
                &mut workspace.taylor,
                |index, value| {
                    if drop_tol > 0.0 && value.abs() < drop_tol {
                        return f64::NEG_INFINITY;
                    }
                    let before = history.get(index).copied().unwrap_or(0.0);
                    (before + factor * (value - before)).abs()
                },
            );
            prune_basis(basis, coeffs, drop_tol, protected);
        } else {
            prune_basis(basis, coeffs, drop_tol, protected);
            cap_basis_with_scratch(basis, coeffs, max_basis, protected, &mut workspace.taylor);
        }
        Ok(t)
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
    /// post-step basis is trimmed to the top-`max_basis` reps by retention
    /// score (`|c|` by default).
    /// Pass a large value (e.g. `usize::MAX`) for the near-exact,
    /// uncapped case. `drop_tol` additionally prunes by magnitude.
    /// `protected` reps are never dropped.
    ///
    /// `basis` is assumed to contain only canonical orbit
    /// representatives. If not, call
    /// [`canonicalize_basis_to_rep`](crate::canonicalize_basis_to_rep)
    /// first.
    ///
    /// Honours `cfg.num_threads` and `cfg.forward_looking` the same way
    /// [`Self::pc_step`] does, using the complex modulus for retention scores.
    pub fn pc_step_orbit_rep(
        &self,
        basis: &mut Vec<Word<C>>,
        coeffs: &mut Vec<Complex<f64>>,
        dt: f64,
        protected: &[Word<C>],
        sector: &Sector<'_>,
        cfg: &PcStepConfig,
    ) -> Result<(), Error> {
        cfg.ranking_factor(dt)?;
        if basis.len() != coeffs.len() {
            return Err(Error::LengthMismatch {
                what: "basis and coeffs",
                a: basis.len(),
                b: coeffs.len(),
            });
        }
        self.run_in_pool(cfg, |this| {
            this.pc_step_orbit_rep_inner(basis, coeffs, dt, protected, sector, cfg)
        })
    }

    fn pc_step_orbit_rep_inner(
        &self,
        basis: &mut Vec<Word<C>>,
        coeffs: &mut Vec<Complex<f64>>,
        dt: f64,
        protected: &[Word<C>],
        sector: &Sector<'_>,
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

        // Appended orbit representatives have zero pre-step coefficients.
        let ranking_factor = cfg.ranking_factor(dt)?;
        let history = ranking_factor.map(|_| coeffs.clone());

        // 1. First-hop phase-aware leakage.
        let mut leak = self.leakage_orbit_rep(basis, coeffs, protected, sector, admit)?;
        if tau_add > 0.0 {
            leak.retain(|(_, c)| c.norm() > tau_add);
        }
        add_leakage_capped(basis, coeffs, leak, admit);

        // 2. Predictor: the phase-aware action is built once and reused
        //    across every matvec.
        let coeffs_predict = mf_expm::expm_apply_orbit_rep(self, basis, sector, dt, coeffs);

        // 3. Second-hop leakage from the predicted state, skipped when the
        //    first hop left no admission room (see `pc_step_inner`).
        let n_predict = basis.len();
        if admit > n_predict {
            let mut leak2 =
                self.leakage_orbit_rep(basis, &coeffs_predict, protected, sector, admit)?;
            if tau_add > 0.0 {
                leak2.retain(|(_, c)| c.norm() > tau_add);
            }
            add_leakage_capped(basis, coeffs, leak2, admit);
        }

        // 4. Corrector: redo from the pre-step state if the basis grew;
        //    otherwise it would reproduce the predictor exactly.
        if basis.len() == n_predict {
            *coeffs = coeffs_predict;
        } else {
            drop(coeffs_predict);
            *coeffs = mf_expm::expm_apply_orbit_rep(self, basis, sector, dt, coeffs);
        }

        // 5. Score before compaction to preserve the pre-step coefficient positions.
        if let (Some(factor), Some(history)) = (ranking_factor, history) {
            cap_basis_by_score(
                basis,
                coeffs,
                max_basis,
                protected,
                &mut Vec::new(),
                |index, value| {
                    if drop_tol > 0.0 && value.norm() < drop_tol {
                        return f64::NEG_INFINITY;
                    }
                    let before = history.get(index).copied().unwrap_or_default();
                    (before + (value - before) * factor).norm()
                },
            );
            prune_basis(basis, coeffs, drop_tol, protected);
            return Ok(());
        }
        prune_basis(basis, coeffs, drop_tol, protected);
        cap_basis(basis, coeffs, max_basis, protected);
        Ok(())
    }
}

#[cfg(test)]
mod default_tests {
    use super::*;
    use crate::{JumpInput, word_from_codes};

    #[test]
    fn default_step_matches_original_with_general_jumps_and_kossakowski() {
        let h = vec![
            ("XXI".into(), 0.3),
            ("IYY".into(), -0.17),
            ("XXI".into(), 0.13),
        ];
        let jumps = vec![JumpInput {
            lincomb: vec![
                ("XII".into(), Complex::new(0.5, 0.0)),
                ("YII".into(), Complex::new(0.0, 0.5)),
            ],
            rate: 0.4,
        }];
        let mut reference = LindbladSpec::<2>::new(3, &h, &jumps).unwrap();
        let mut optimized = LindbladSpec::<2>::new(3, &h, &jumps).unwrap();
        let operators = vec![
            vec![("IXI".into(), Complex::new(1.0, 0.0))],
            vec![("IYI".into(), Complex::new(1.0, 0.0))],
        ];
        let k = vec![
            vec![Complex::new(0.1, 0.0), Complex::new(0.0, 0.02)],
            vec![Complex::new(0.0, -0.02), Complex::new(0.2, 0.0)],
        ];
        reference.add_kossakowski(&operators, &k).unwrap();
        optimized.add_kossakowski(&operators, &k).unwrap();
        let seed = word_from_codes::<2>(&[2, 0, 0]).unwrap();
        let (mut a, mut b) = (vec![seed], vec![seed]);
        let (mut ca, mut cb) = (vec![1.0], vec![1.0]);
        for step in 0..12 {
            let cfg = PcStepConfig {
                max_basis: 8 + step % 3,
                admit_basis: Some(24 + step % 4),
                drop_tol: 1e-7,
                num_threads: Some(2),
                ..Default::default()
            };
            reference
                .pc_step_reference(&mut a, &mut ca, 0.1, &[seed], &cfg)
                .unwrap();
            optimized
                .pc_step(&mut b, &mut cb, 0.1, &[seed], &cfg)
                .unwrap();
            assert_eq!(a, b);
            assert_eq!(ca, cb);
        }
        let pool = optimized.workspaces.lock().unwrap();
        assert_eq!(pool.len(), 1);
        assert!(pool[0].taylor.capacity() > 0);
        assert!(pool[0].generator.is_none());
    }

    #[test]
    fn concurrent_default_trajectories_share_model_without_sharing_scratch() {
        let spec =
            LindbladSpec::<2>::new(3, &[("XXI".into(), 0.3), ("IYY".into(), -0.17)], &[]).unwrap();
        let barrier = std::sync::Barrier::new(4);
        std::thread::scope(|scope| {
            for trajectory in 0..4 {
                let spec = &spec;
                let barrier = &barrier;
                scope.spawn(move || {
                    let seed = word_from_codes::<2>(&[trajectory % 3 + 1, 0, 0]).unwrap();
                    let cfg = PcStepConfig {
                        max_basis: 8 + trajectory as usize,
                        admit_basis: Some(24),
                        ..Default::default()
                    };
                    let (mut a, mut b) = (vec![seed], vec![seed]);
                    let (mut ca, mut cb) = (vec![1.0], vec![1.0]);
                    barrier.wait();
                    for _ in 0..6 {
                        spec.pc_step_reference(&mut a, &mut ca, 0.1, &[seed], &cfg)
                            .unwrap();
                        spec.pc_step(&mut b, &mut cb, 0.1, &[seed], &cfg).unwrap();
                        assert_eq!(a, b);
                        assert_eq!(ca, cb);
                    }
                });
            }
        });
        let pool = spec.workspaces.lock().unwrap();
        assert!(!pool.is_empty());
        assert!(pool.len() <= 4);
    }
    #[test]
    fn large_rank_limit_does_not_eagerly_allocate_for_a_small_state() {
        let n = 129;
        let mut z = vec!['I'; n];
        z[0] = 'Z';
        let spec = LindbladSpec::<4>::new(n, &[(z.into_iter().collect(), 0.37)], &[]).unwrap();
        let mut codes = vec![0; n];
        codes[0] = 2;
        let seed = word_from_codes::<4>(&codes).unwrap();
        let mut basis = vec![seed];
        let mut coeffs = vec![1.0];
        let cfg = PcStepConfig {
            max_basis: 1 << 30,
            admit_basis: Some(3 << 30),
            ..Default::default()
        };
        spec.pc_step(&mut basis, &mut coeffs, 0.1, &[seed], &cfg)
            .unwrap();
        assert_eq!(basis, vec![seed]);
        assert_eq!(coeffs, vec![1.0]);
        assert!(basis.capacity() <= 8);
        assert!(coeffs.capacity() <= 8);
        let pool = spec.workspaces.lock().unwrap();
        assert!(pool[0].taylor.capacity() <= 16);
        assert!(pool[0].coefficient_auxiliary.capacity() <= 8);
    }
}
