//! Configuration objects for the predictor-corrector stepper.

/// Truncation and execution policy for a single predictor-corrector step
/// ([`crate::LindbladSpec::pc_step`], [`crate::LindbladSpec::pc_step_timed`],
/// [`crate::orbit_rep::pc_step_orbit_rep`]).
///
/// These are the per-run *tuning knobs*, kept separate from the per-call data
/// (`basis`, `coeffs`, `dt`, `protected`, and — on the orbit path — the
/// translation group and momentum).
///
/// `max_basis` is the primary accuracy/cost dial; `admit_basis` selects the
/// displacement scheme; `drop_tol` is the churn valve of the admission-bound
/// scheme; `tau_add` is a wall optimization at most.
#[derive(Debug, Clone, Copy)]
pub struct PcStepConfig {
    /// Hard rank cap on the retained basis: after the corrector, only the
    /// top-`max_basis` strings by `|coeff|` are kept (protected words always
    /// survive). The primary convergence dial — verify by re-running at 2×.
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

impl Default for PcStepConfig {
    /// Uncapped, unfiltered, no pruning: the near-exact reference
    /// configuration. Production runs should set `max_basis` (and usually
    /// `admit_basis ≈ 2-3×` it).
    fn default() -> Self {
        Self {
            max_basis: usize::MAX,
            admit_basis: None,
            drop_tol: 0.0,
            tau_add: None,
            num_threads: None,
        }
    }
}
