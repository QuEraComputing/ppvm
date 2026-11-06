// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Configuration objects for the predictor-corrector stepper.

/// Built-in candidate generation for real-coefficient predictor-corrector steps.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum CandidateStrategy {
    /// Chunk-capped signed leakage accumulation.
    #[default]
    Streaming,
    /// Shortlist individual contributions, then score their complete signed leakage.
    Shortlist,
}

/// Truncation and execution policy for a single predictor-corrector step
/// ([`crate::LindbladSpec::pc_step`], [`crate::LindbladSpec::pc_step_timed`]).
///
/// These are the per-run *tuning knobs*, kept separate from the per-call data
/// (`basis`, `coeffs`, `dt`, `protected`).
///
/// `max_basis` is the primary accuracy/cost dial; `admit_basis` selects the
/// displacement scheme; `drop_tol` is the churn valve of the admission-bound
/// scheme; `tau_add` is a wall optimization at most.
#[derive(Debug, Clone, Copy)]
pub struct PcStepConfig {
    /// Candidate generation for real-coefficient steps. Explicit proposer calls
    /// use their supplied strategy instead of this field.
    pub proposer: CandidateStrategy,
    /// Rank coefficients over a fixed interval when `Some(h)`, using
    /// `|c_before + h/dt * (c_after-c_before)|`. Actual evolved values are kept.
    /// `None` uses ordinary magnitude ranking; `Some(dt)` is identical to it.
    /// The interval starts at the beginning of this step, not its endpoint.
    /// Prediction uses only this step's coefficient change, with no history
    /// across calls and no evaluation of future generator changes.
    pub forward_looking: Option<f64>,
    /// Hard rank cap on the retained basis: after the corrector, only the
    /// top-`max_basis` strings by retention score are kept (`|coeff|` by
    /// default; protected words always survive). The primary convergence
    /// dial — verify by re-running at 2×.
    pub max_basis: usize,
    /// Working-set (admission) bound. When `Some(a)` with `a > max_basis`,
    /// enrichment may grow the live basis to `a` and the final cap performs a
    /// genuine top-`max_basis`-of-union rank displacement (the analog of
    /// two-site TDVP truncation at `χ_max`); `drop_tol` is then not needed
    /// for membership turnover. `None` bounds admission by `max_basis`
    /// itself — the valve scheme, which requires `drop_tol > 0` to keep the
    /// basis adapting once it fills.
    pub admit_basis: Option<usize>,
    /// Magnitude prune applied after the corrector: basis entries whose
    /// `|coeff|` is below `drop_tol` are discarded (protected words are
    /// always kept). `<= 0.0` disables pruning — valid only with
    /// `admit_basis` set, otherwise the basis freezes once it fills the cap.
    pub drop_tol: f64,
    /// Optional absolute rate threshold on leakage admission: a candidate is
    /// admitted only if its inflow rate exceeds `tau_add`. This is the
    /// natural (dt- and drop_tol-independent) parameterization — the
    /// admission accuracy cliff sits at a fixed `tau_add`. `None` = no
    /// filter, the recommended default with cap-based truncation.
    pub tau_add: Option<f64>,
    /// When `Some(n)`, run the entire step inside a freshly built rayon
    /// thread pool of `n` threads (useful for benchmarking parallel
    /// scaling). When `None`, the global rayon pool is used.
    pub num_threads: Option<usize>,
}

impl PcStepConfig {
    pub(crate) fn ranking_factor(&self, dt: f64) -> Result<Option<f64>, crate::Error> {
        match self.forward_looking {
            Some(horizon) => {
                if !dt.is_finite()
                    || dt <= 0.0
                    || !horizon.is_finite()
                    || horizon <= 0.0
                    || !(horizon / dt).is_finite()
                {
                    return Err(crate::Error::InvalidRankingHorizon);
                }
                Ok((horizon != dt).then_some(horizon / dt))
            }
            None => Ok(None),
        }
    }

    /// Set the working-set cap to `ceil(ratio * max_basis)`.
    ///
    /// The ratio scales the total admission bound, including the retained
    /// basis. For example, `max_basis = 1025` and `ratio = 1.5` admit at most
    /// 1538 strings before truncation. Existing absolute caps are overwritten.
    ///
    /// # Errors
    ///
    /// Returns an error for nonfinite ratios, ratios below one, or caps that
    /// cannot be represented as a `usize`.
    pub fn with_admit_ratio(mut self, ratio: f64) -> Result<Self, crate::Error> {
        if !ratio.is_finite() || ratio < 1.0 {
            return Err(crate::Error::InvalidAdmissionRatio { ratio });
        }
        // Integer multipliers preserve exact counts even beyond f64's integer precision.
        if ratio.fract() == 0.0 && ratio < usize::MAX as f64 {
            self.admit_basis = Some(
                self.max_basis
                    .checked_mul(ratio as usize)
                    .ok_or(crate::Error::AdmissionCapOverflow)?,
            );
            return Ok(self);
        }
        let cap = (self.max_basis as f64 * ratio).ceil();
        if !cap.is_finite() || cap as u128 > usize::MAX as u128 {
            return Err(crate::Error::AdmissionCapOverflow);
        }
        self.admit_basis = Some(cap as usize);
        Ok(self)
    }
}

impl Default for PcStepConfig {
    /// Uncapped, unfiltered, no pruning: the near-exact reference
    /// configuration. Production runs should set `max_basis` (and usually
    /// `admit_basis ≈ 2-3×` it).
    fn default() -> Self {
        Self {
            proposer: CandidateStrategy::Streaming,
            forward_looking: None,
            max_basis: usize::MAX,
            admit_basis: None,
            drop_tol: 0.0,
            tau_add: None,
            num_threads: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::PcStepConfig;

    #[test]
    fn fractional_admission_rounds_up() {
        for (basis, ratio, expected) in [
            (1025, 1.5, 1538),
            (1025, 1.234, 1265),
            (1024, 3.0, 3072),
            (1024, 6.0, 6144),
            (1024, 1.0, 1024),
        ] {
            let cfg = PcStepConfig {
                max_basis: basis,
                ..Default::default()
            }
            .with_admit_ratio(ratio)
            .unwrap();
            assert_eq!(cfg.admit_basis, Some(expected));
            assert_eq!(cfg.max_basis, basis);
        }
    }

    #[test]
    fn invalid_admission_ratios_and_overflow_are_rejected() {
        for ratio in [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            -1.0,
            0.0,
            0.5,
            f64::MAX,
        ] {
            assert!(
                PcStepConfig {
                    max_basis: 1024,
                    ..Default::default()
                }
                .with_admit_ratio(ratio)
                .is_err()
            );
        }
        assert!(PcStepConfig::default().with_admit_ratio(2.0).is_err());
        assert_eq!(
            PcStepConfig::default()
                .with_admit_ratio(1.0)
                .unwrap()
                .admit_basis,
            Some(usize::MAX)
        );
    }

    #[test]
    fn admission_ratio_preserves_other_policy_settings() {
        let cfg = PcStepConfig {
            max_basis: 10,
            admit_basis: Some(40),
            drop_tol: 0.1,
            tau_add: Some(0.2),
            num_threads: Some(2),
            ..Default::default()
        }
        .with_admit_ratio(1.5)
        .unwrap();
        assert_eq!(cfg.admit_basis, Some(15));
        assert_eq!(cfg.drop_tol, 0.1);
        assert_eq!(cfg.tau_add, Some(0.2));
        assert_eq!(cfg.num_threads, Some(2));
    }
}
