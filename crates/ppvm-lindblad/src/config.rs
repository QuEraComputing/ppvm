// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Configuration objects for the predictor-corrector stepper.

/// Policy options for a single predictor-corrector step
/// ([`crate::LindbladSpec::pc_step`] / [`crate::LindbladSpec::pc_step_timed`]).
///
/// These are the per-run *tuning knobs*, kept separate from the per-call data
/// (`basis`, `coeffs`, `dt`, `protected`) so the step functions don't grow an
/// unwieldy positional-argument list.
#[derive(Debug, Clone, Copy)]
pub struct PcStepConfig {
    /// Leakage-rate threshold for basis growth: a leakage Pauli is appended to
    /// the basis only if its `|coeff|` exceeds `tau_add`. Larger values keep the
    /// basis smaller (faster, less accurate); it also drives the streaming
    /// Cauchy-Schwarz prune in [`crate::LindbladSpec::leakage_with_prune`].
    pub tau_add: f64,
    /// Magnitude prune applied after the corrector: basis entries whose
    /// `|coeff|` is below `drop_tol` are discarded (protected words are always
    /// kept). `drop_tol <= 0.0` disables pruning.
    pub drop_tol: f64,
    /// When `Some(n)`, run the entire step inside a freshly built rayon thread
    /// pool of `n` threads (useful for benchmarking parallel scaling). When
    /// `None`, the global rayon pool is used.
    pub num_threads: Option<usize>,
}
