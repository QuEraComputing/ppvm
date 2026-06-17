// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Regression guard for `CoefficientThreshold` truncation semantics under
//! *deferred* truncation.
//!
//! `PauliSum::map_insert` must never drop a gate-produced term on its own;
//! terms are only ever removed by an explicit `truncate()`. A since-reverted
//! optimization (the `Strategy::discard` insertion-pruning) skipped inserting
//! a produced term whose coefficient was individually below the threshold and
//! whose key was not yet present, on the assumption that the *next* call is
//! always `truncate()`. That assumption breaks when a caller applies several
//! gates between truncations: two sub-threshold contributions to the same key
//! can sum to an above-threshold coefficient (`|a + b| ≥ τ` while
//! `|a|, |b| < τ`), so dropping either at production loses a genuinely
//! surviving term — and biases the terms it would have fed.
//!
//! This test pins the correct behavior: with truncation deferred across
//! gates, the `CoefficientThreshold` result must equal a non-truncating
//! reference filtered once at the end (standard "insert everything, truncate
//! last" semantics).

use ppvm_runtime::prelude::*;
use ppvm_runtime::strategy::CoefficientThreshold;

macro_rules! sorted_terms {
    ($state:expr) => {{
        let mut v: Vec<(String, f64)> = $state.iter().map(|(k, c)| (k.to_string(), *c)).collect();
        v.sort_by(|a, b| a.0.cmp(&b.0));
        v
    }};
}

#[test]
fn deferred_truncation_keeps_terms_whose_merged_coefficient_survives() {
    // Two rx(θ) on a single Z with no truncate in between. The first rx
    // produces Y at ~sin(θ) ≈ 0.030 (< τ); the second adds another ~0.030, so
    // the accumulated Y ≈ sin(2θ) ≈ 0.060 is *above* τ and must survive the
    // final truncate. Insertion-pruning would drop the first (absent,
    // sub-threshold) Y and never accumulate it — losing the term entirely.
    const TAU: f64 = 0.05;
    const THETA: f64 = 0.03;

    // CoefficientThreshold, truncation deferred until the very end.
    let mut approx: PauliSum<config::indexmap::ByteFxHashF64<2, CoefficientThreshold>> =
        PauliSum::builder()
            .n_qubits(1)
            .strategy(CoefficientThreshold(TAU))
            .build();
    approx += ("Z", 1.0);
    approx.rx(0, THETA);
    approx.rx(0, THETA);
    approx.truncate();

    // Exact reference: NoStrategy never truncates; filter once at τ at the end.
    let mut exact: PauliSum<config::indexmap::ByteFxHashF64<2>> =
        PauliSum::builder().n_qubits(1).build();
    exact += ("Z", 1.0);
    exact.rx(0, THETA);
    exact.rx(0, THETA);
    let reference: Vec<(String, f64)> = sorted_terms!(exact)
        .into_iter()
        .filter(|(_, c)| c.abs() >= TAU)
        .collect();

    // Sanity: the merged Y term really is above threshold (so this guards a
    // dropped surviving term, not just a coefficient nudge).
    assert!(
        reference.iter().any(|(k, c)| k == "Y" && c.abs() >= TAU),
        "test setup broken: reference should contain an above-threshold Y term, got {reference:?}"
    );

    assert_eq!(
        sorted_terms!(approx),
        reference,
        "CoefficientThreshold with deferred truncation dropped or biased a term \
         that standard insert-then-truncate keeps"
    );
}
