// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! End-to-end guard for the Trotter benchmark circuit
//! (`crates/ppvm-runtime/benches/trotter.rs`).
//!
//! The benchmark itself only times `trotter_func` — it never reads the
//! propagated observable back, so nothing there guards the *result*. These
//! tests run the same gate sequence (`rx` + `rzz` + `pauli_error`, truncating
//! after every gate as the bench does) on a small chain and check:
//!
//! 1. **Math regression** — the exact (untruncated, `NoStrategy`) expectation
//!    value matches a frozen golden constant.
//! 2. **Truncation fidelity** — the `CoefficientThreshold(1e-6)` run stays
//!    within a tight tolerance of the exact run, confirming magnitude
//!    truncation drops only genuinely negligible mass: an over-aggressive
//!    approximation (e.g. dropping sub-threshold contributions to surviving
//!    terms) would drift past the bound and fail here.

use ppvm_runtime::prelude::*;
use ppvm_runtime::strategy::CoefficientThreshold;

const N: usize = 4;
const STEPS: usize = 10;
const THETA_X: f64 = 0.1; // dt * external_field, with dt = 0.1, h = 1.0
const THETA_ZZ: f64 = 0.0125; // dt * interaction_strength, with j = 1/8
const NOISE: [f64; 3] = [2.5e-5; 3]; // 1e-4 / 4 per channel, as in the bench

/// Apply `STEPS` first-order Trotter steps of the noisy XZZ Ising chain,
/// truncating after every gate (identical structure to the benchmark's
/// `trotter_func`). A macro rather than a generic fn so it works verbatim on
/// both the `NoStrategy` and `CoefficientThreshold` configs without replaying
/// the full method-bound list.
macro_rules! trotter_evolve {
    ($state:expr) => {{
        for _ in 0..STEPS {
            for i in 0..N {
                $state.rx(i, THETA_X);
                $state.truncate();
                $state.pauli_error(i, NOISE);
                $state.truncate();
            }
            for i in 0..N - 1 {
                $state.rzz(i, i + 1, THETA_ZZ);
                $state.truncate();
                $state.pauli_error(i, NOISE);
                $state.truncate();
                $state.pauli_error(i + 1, NOISE);
                $state.truncate();
            }
        }
    }};
}

/// `⟨0…0| O |0…0⟩` for the propagated observable `O`: on the all-zero state
/// only Paulis built from `I`/`Z` survive (each contributes its coefficient).
macro_rules! expect_on_zero {
    ($state:expr) => {{
        let mut acc = 0.0_f64;
        for (k, v) in $state.data().iter() {
            if k.to_string().chars().all(|c| c == 'I' || c == 'Z') {
                acc += *v;
            }
        }
        acc
    }};
}

/// Seed the observable `Σ_i Z_i` (matching the benchmark's initial state).
macro_rules! seed_sum_z {
    ($state:expr) => {{
        for j in 0..N {
            let term: String = (0..N).map(|i| if i == j { 'Z' } else { 'I' }).collect();
            $state += (term.as_str(), 1.0);
        }
    }};
}

#[test]
fn trotter_result_is_stable_and_truncation_is_faithful() {
    // Exact reference: no truncation strategy, so the full Pauli support is
    // kept (bounded by 4^N terms for a fixed N).
    let mut exact: PauliSum<config::indexmap::ByteFxHashF64<4>> =
        PauliSum::builder().n_qubits(N).build();
    seed_sum_z!(exact);
    trotter_evolve!(exact);
    let exact_val = expect_on_zero!(exact);

    // Truncated run: identical circuit, magnitude truncation at 1e-6.
    let mut approx: PauliSum<config::indexmap::ByteFxHashF64<4, CoefficientThreshold>> =
        PauliSum::builder()
            .n_qubits(N)
            .strategy(CoefficientThreshold(1e-6))
            .build();
    seed_sum_z!(approx);
    trotter_evolve!(approx);
    let approx_val = expect_on_zero!(approx);

    // (1) Math-regression guard: frozen golden for the exact circuit. Any
    // change to the gate math (a wrong sign, bit-flip, or addressing bug in
    // rx/rzz/pauli_error) moves this by O(1e-3) or more; the bound is far
    // tighter than that yet far looser than cross-platform trig last-bit noise.
    const GOLDEN: f64 = 2.161056303575631;
    assert!(
        (exact_val - GOLDEN).abs() < 1e-9,
        "exact Trotter expectation {exact_val} drifted from golden {GOLDEN}"
    );

    // (2) Truncation-fidelity guard: the truncated run must stay within one
    // truncation unit (1e-6) of exact. Standard magnitude truncation drifts
    // only ~2.4e-8 here. An over-aggressive approximation that drops
    // sub-threshold contributions to *surviving* terms drifts ~2.2e-6 and
    // trips this bound, which is the regression we guard.
    const TOL: f64 = 1e-6;
    let drift = (exact_val - approx_val).abs();
    assert!(
        drift < TOL,
        "truncated result {approx_val} drifted from exact {exact_val} by {drift:e} (> {TOL:e})"
    );
}
