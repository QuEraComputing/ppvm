// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Tests for `GeneralizedTableau::project` (post-selected Z measurement).

use std::f64::consts::FRAC_PI_2;

use ppvm_pauli_sum::config::dashmap::ByteFxHashF64;
use ppvm_tableau::measure::ProjectError;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<ByteFxHashF64<1>, usize>;

const TOL: f64 = 1e-12;

/// Full observable snapshot of a tableau: stabilizer frame (via `Display`),
/// coefficients sorted by index, and the measurement record.
type Snapshot = (String, Vec<(usize, f64, f64)>, Vec<Option<bool>>);

fn snapshot(tab: &Tab) -> Snapshot {
    let mut coeffs: Vec<_> = tab
        .coefficients
        .iter()
        .map(|(c, i)| (*i, c.re, c.im))
        .collect();
    coeffs.sort_by_key(|e| e.0);
    (
        format!("{}", tab.tableau),
        coeffs,
        tab.measurement_record.clone(),
    )
}

/// `ry(θ)|0⟩`, with `Z_0` either in the stabilizer group or anticommuting
/// with a stabilizer (via `H`); measurement collapses these differently.
fn ry_state(theta: f64, z_anticommutes: bool) -> Tab {
    let mut tab = Tab::new_with_seed(1, 1e-12, 0);
    if z_anticommutes {
        // H|0⟩ = ry(π/2)|0⟩, so H then ry(θ − π/2) is ry(θ)|0⟩.
        tab.h(0);
        tab.ry(0, theta - FRAC_PI_2);
    } else {
        tab.ry(0, theta);
    }
    let (_, stab_anticomm_bits, _) = tab.compute_decomposition(0, Pauli::Z);
    assert_eq!(stab_anticomm_bits != 0, z_anticommutes);
    tab
}

#[test]
fn project_returns_born_probability_and_collapses() {
    for z_anticommutes in [false, true] {
        for theta in [0.3f64, 1.0, 2.0, 2.9] {
            let p1 = (theta / 2.0).sin().powi(2);

            let mut tab = ry_state(theta, z_anticommutes);
            let p = tab.project(0, true, TOL).unwrap();
            assert!((p - p1).abs() < 1e-10, "θ={theta}: {p} vs {p1}");
            assert!((tab.z_expectation(0) + 1.0).abs() < 1e-10);

            let mut tab = ry_state(theta, z_anticommutes);
            let p = tab.project(0, false, TOL).unwrap();
            assert!((p - (1.0 - p1)).abs() < 1e-10);
            assert!((tab.z_expectation(0) - 1.0).abs() < 1e-10);
        }
    }
}

#[test]
fn project_matches_measure_that_landed_on_same_outcome() {
    // A post-selected projection must leave exactly the state that a sampled
    // measurement with the same outcome leaves.
    let mut base = Tab::new_with_seed(3, 1e-12, 7);
    base.h(0);
    base.t(0);
    base.cnot(0, 1);
    base.ry(1, 0.7);
    base.h(2);
    base.t(2);
    base.cz(1, 2);

    for q in 0..3 {
        for seed in 0..20 {
            let mut measured = base.fork(Some(seed));
            let outcome = measured.measure(q).unwrap();
            let mut projected = base.fork(Some(seed));
            projected.project(q, outcome, TOL).unwrap();
            assert_eq!(
                snapshot(&measured),
                snapshot(&projected),
                "q={q} seed={seed}"
            );
        }
    }
}

#[test]
fn ghz_correlations_and_impossible_outcome() {
    let mut tab = Tab::new_with_seed(3, 1e-12, 0);
    tab.h(0);
    tab.t(0);
    tab.cnot(0, 1);
    tab.cnot(1, 2);

    assert!((tab.project(0, true, TOL).unwrap() - 0.5).abs() < 1e-12);
    assert!((tab.project(1, true, TOL).unwrap() - 1.0).abs() < 1e-12);

    assert!(matches!(
        tab.project(2, false, TOL),
        Err(ProjectError::ImpossibleOutcome { .. })
    ));
    assert!((tab.project(2, true, TOL).unwrap() - 1.0).abs() < 1e-12);
    assert_eq!(
        tab.measurement_record,
        vec![Some(true), Some(true), Some(true)]
    );
}

#[test]
fn impossible_outcome_is_rejected_without_mutation() {
    // |0⟩ in an `H`-rotated frame, where collapse rebuilds the coefficient
    // vector; a rejection must hand back the original one.
    let mut tab = ry_state(0.0, true);
    tab.t(0);
    let before = snapshot(&tab);
    assert!(matches!(
        tab.project(0, true, TOL),
        Err(ProjectError::ImpossibleOutcome { .. })
    ));
    assert_eq!(snapshot(&tab), before);
    assert!((tab.project(0, false, TOL).unwrap() - 1.0).abs() < 1e-10);
}

#[test]
fn tol_threshold_is_respected() {
    // P(1) = sin²(θ/2) = 1e-8: accepted with tol = 1e-12, rejected with 1e-6.
    let theta = 2.0 * 1e-4f64.asin();
    let mut tab = ry_state(theta, false);
    assert!(matches!(
        tab.project(0, true, 1e-6),
        Err(ProjectError::ImpossibleOutcome { .. })
    ));
    let p = tab.project(0, true, TOL).unwrap();
    assert!((p - 1e-8).abs() < 1e-14);
}

#[test]
fn lost_qubit_is_rejected() {
    let mut tab = Tab::new_with_seed(1, 1e-12, 0);
    tab.loss_channel(0, 1.0);
    assert_eq!(
        tab.project(0, false, TOL),
        Err(ProjectError::LostQubit { addr0: 0 })
    );
}

#[test]
fn project_after_truncation_divides_by_norm() {
    // Same setup as `test_measure_after_truncation_divides_by_norm`: rx drops
    // a branch without renormalizing, but qubit 0 is still exactly |0⟩.
    let theta = 2.0 * 0.25f64.asin();
    let mut tab = Tab::new_with_seed(1, 0.3, 0);
    tab.rx(0, theta);
    assert!((tab.project(0, false, TOL).unwrap() - 1.0).abs() < 1e-12);
}

#[test]
fn chain_rule_matches_z_string_reference() {
    let n = 3;
    let mut base = Tab::new_with_seed(n, 1e-12, 0);
    base.h(0);
    base.t(0);
    base.cnot(0, 1);
    base.ry(1, 0.9);
    base.rx(2, 1.3);
    base.cz(1, 2);
    base.h(2);
    base.t(2);
    base.cnot(2, 0);

    let mut total = 0.0;
    for z in 0..1usize << n {
        // P(z) = Π_q P(z_q | z_<q); an impossible factor makes P(z) = 0.
        let mut tab = base.fork(None);
        let p: f64 = (0..n)
            .map(|q| tab.project(q, (z >> q) & 1 == 1, TOL).unwrap_or(0.0))
            .product();

        // Reference: ⟨ψ|z⟩⟨z|ψ⟩ = ⟨Π_q (I + (-1)^{z_q} Z_q)/2⟩, expanded into
        // Z-string expectations.
        let reference: f64 = (0..1usize << n)
            .map(|subset| {
                let word: String = (0..n)
                    .map(|q| if (subset >> q) & 1 == 1 { 'Z' } else { 'I' })
                    .collect();
                let word: PauliWord<<ByteFxHashF64<1> as Config>::Storage> = word.into();
                let sign = if (z & subset).count_ones() % 2 == 0 {
                    1.0
                } else {
                    -1.0
                };
                sign * base.expectation(&word)
            })
            .sum::<f64>()
            / (1usize << n) as f64;

        assert!(
            (p - reference).abs() < 1e-10,
            "z={z:03b}: {p} vs {reference}"
        );
        total += p;
    }
    assert!((total - 1.0).abs() < 1e-10, "Σ P(z) = {total}");
}
