// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::f64::consts::FRAC_PI_2;

use ppvm_pauli_sum::config::fxhash::ByteF64;
use ppvm_tableau::prelude::*;

type TestTableau = GeneralizedTableau<ByteF64<1>>;

fn word(s: &str) -> PauliWord<u64> {
    s.into()
}

fn assert_close(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() < tol,
        "expected {expected}, got {actual} (|Δ| = {})",
        (actual - expected).abs()
    );
}

/// All `4^n` Pauli strings on `n` qubits.
fn all_paulis(n: usize) -> Vec<String> {
    let mut out = vec![String::new()];
    for _ in 0..n {
        out = out
            .into_iter()
            .flat_map(|s| ["I", "X", "Y", "Z"].map(|p| format!("{s}{p}")))
            .collect();
    }
    out
}

/// `P(b) = 2^-n Σ_T (-1)^{b·T} ⟨Z_T⟩`, computed without collapsing the state.
fn bitstring_prob_from_expectations(tab: &TestTableau, bits: &[bool]) -> f64 {
    let n = bits.len();
    let mut sum = 0.0;
    for mask in 0..(1usize << n) {
        let w: String = (0..n)
            .map(|q| if mask >> q & 1 == 1 { 'Z' } else { 'I' })
            .collect();
        let sign = (0..n).filter(|&q| mask >> q & 1 == 1 && bits[q]).count() % 2;
        let e = tab.expectation(&word(&w));
        sum += if sign == 0 { e } else { -e };
    }
    sum / (1usize << n) as f64
}

/// Chain-rule bitstring probability via successive projections.
fn bitstring_prob_by_projection(tab: &TestTableau, bits: &[bool]) -> f64 {
    let mut t = tab.fork(Some(0));
    let mut prob = 1.0;
    for (q, &b) in bits.iter().enumerate() {
        match t.project(q, b) {
            Ok(p) => prob *= p,
            Err(ProjectError::ZeroProbability { .. }) => return 0.0,
            Err(e) => panic!("unexpected error: {e}"),
        }
    }
    prob
}

fn all_bitstrings(n: usize) -> Vec<Vec<bool>> {
    (0..(1usize << n))
        .map(|m| (0..n).map(|q| m >> q & 1 == 1).collect())
        .collect()
}

/// A small entangled, non-Clifford circuit exercising both measurement cases.
fn magic_circuit() -> TestTableau {
    let mut tab = TestTableau::new(3, 1e-12);
    tab.h(0);
    tab.t(0);
    tab.cnot(0, 1);
    tab.h(2);
    tab.t(2);
    tab.cnot(1, 2);
    tab.ry(1, 0.7);
    tab.t(1);
    tab.h(1);
    tab.rx(2, 0.3);
    tab
}

#[test]
fn project_zero_state() {
    let mut tab = TestTableau::new(1, 1e-12);
    assert_close(tab.project(0, false).unwrap(), 1.0, 1e-12);
    assert_eq!(tab.current_measurement_record(), &[Some(false)]);
}

/// True if projecting `addr0` takes the case-a path (Z anticommutes with some
/// stabilizer), false for case b (Z is a stabilizer up to sign).
fn is_case_a(tab: &TestTableau, addr0: usize) -> bool {
    let (_, stab_anticomm_bits, _) = tab.compute_decomposition(addr0, Pauli::Z);
    stab_anticomm_bits != 0
}

/// Projects `addr0` onto an outcome that must have zero probability, and checks
/// that the error is returned and the tableau, coefficients, Pauli expectations,
/// and measurement record are all unchanged.
fn assert_zero_probability_preserves_state(tab: &mut TestTableau, addr0: usize, outcome: bool) {
    let n = tab.n_qubits();
    let paulis = all_paulis(n);
    let before_display = format!("{tab}");
    let before_expectations: Vec<f64> = paulis.iter().map(|p| tab.expectation(&word(p))).collect();
    let before_record = tab.current_measurement_record().to_vec();

    let err = tab.project(addr0, outcome).unwrap_err();
    assert!(
        matches!(err, ProjectError::ZeroProbability { addr0: a, outcome: o, .. } if a == addr0 && o == outcome),
        "unexpected error: {err:?}"
    );

    assert_eq!(format!("{tab}"), before_display);
    for (p, before) in paulis.iter().zip(before_expectations) {
        assert_close(tab.expectation(&word(p)), before, 1e-12);
    }
    assert_eq!(tab.current_measurement_record(), before_record.as_slice());
}

#[test]
fn project_zero_probability_case_b() {
    // |0⟩: Z is a stabilizer, P(1) = 0.
    let mut tab = TestTableau::new(1, 1e-12);
    assert!(!is_case_a(&tab, 0));
    assert_zero_probability_preserves_state(&mut tab, 0, true);

    // |+⟩ with a T phase, projected onto 0: the projection makes Z a
    // stabilizer, so the impossible second projection also goes through case b.
    let mut tab = TestTableau::new(1, 1e-12);
    tab.h(0);
    tab.t(0);
    assert!(is_case_a(&tab, 0));
    tab.project(0, false).unwrap();
    assert!(!is_case_a(&tab, 0));
    assert_zero_probability_preserves_state(&mut tab, 0, true);
}

#[test]
fn project_zero_probability_case_a() {
    // H then RY(-π/2) returns to |0⟩, but the rotation only reweights the
    // coefficients: the stabilizer frame still holds X, so Z is not a stabilizer
    // and the projection takes case a even though P(1) = 0.
    let mut tab = TestTableau::new(1, 1e-12);
    tab.h(0);
    tab.ry(0, -FRAC_PI_2);
    assert!(is_case_a(&tab, 0));
    assert_close(tab.expectation(&word("Z")), 1.0, 1e-12);
    assert_zero_probability_preserves_state(&mut tab, 0, true);
    // The possible outcome still projects, with probability 1.
    assert_close(tab.project(0, false).unwrap(), 1.0, 1e-12);

    // Same on 3 qubits, with qubit 0 in |0⟩ (case a) next to an entangled,
    // non-Clifford pair on qubits 1 and 2.
    let mut tab = TestTableau::new(3, 1e-12);
    tab.h(0);
    tab.h(1);
    tab.h(2);
    tab.t(1);
    tab.cz(1, 2);
    tab.ry(0, -FRAC_PI_2);
    assert!(is_case_a(&tab, 0));
    assert!(tab.coefficients.len() > 1);
    assert_zero_probability_preserves_state(&mut tab, 0, true);
    assert_close(tab.project(0, false).unwrap(), 1.0, 1e-12);
}

#[test]
fn project_plus_state() {
    for outcome in [false, true] {
        let mut tab = TestTableau::new(1, 1e-12);
        tab.h(0);
        assert_close(tab.project(0, outcome).unwrap(), 0.5, 1e-12);
        let z = if outcome { -1.0 } else { 1.0 };
        assert_close(tab.expectation(&word("Z")), z, 1e-12);
        // The outcome is now deterministic.
        assert_close(tab.project(0, outcome).unwrap(), 1.0, 1e-12);
        assert_eq!(
            tab.current_measurement_record(),
            &[Some(outcome), Some(outcome)]
        );
    }
}

#[test]
fn project_ry_matches_cos_squared() {
    for theta in [0.1, 0.7, 1.3, 2.0, 2.9] {
        for outcome in [false, true] {
            let mut tab = TestTableau::new(1, 1e-12);
            tab.ry(0, theta);
            let p0 = (theta / 2.0_f64).cos().powi(2);
            let expected = if outcome { 1.0 - p0 } else { p0 };
            assert_close(tab.project(0, outcome).unwrap(), expected, 1e-10);
            let z = if outcome { -1.0 } else { 1.0 };
            assert_close(tab.expectation(&word("Z")), z, 1e-10);
        }
    }
}

#[test]
fn project_matches_measure() {
    let base = magic_circuit();
    let paulis = all_paulis(3);
    for q in 0..3 {
        for outcome in [false, true] {
            let mut projected = base.fork(Some(0));
            if projected.project(q, outcome).is_err() {
                continue;
            }
            let mut measured = (0..1000u64)
                .map(|s| {
                    let mut t = base.fork(Some(s));
                    let m = t.measure(q);
                    (t, m)
                })
                .find(|(_, m)| *m == Some(outcome))
                .expect("outcome with nonzero probability should be sampled")
                .0;
            for p in &paulis {
                assert_close(
                    projected.expectation(&word(p)),
                    measured.expectation(&word(p)),
                    1e-9,
                );
            }
            // Subsequent measurements on the two states agree in distribution;
            // check that the next qubit's probabilities agree exactly.
            let next = (q + 1) % 3;
            let p_proj = projected.fork(Some(1)).project(next, false);
            let p_meas = measured.project(next, false);
            match (p_proj, p_meas) {
                (Ok(a), Ok(b)) => assert_close(a, b, 1e-9),
                (Err(_), Err(_)) => {}
                other => panic!("mismatch: {other:?}"),
            }
        }
    }
}

#[test]
fn project_chain_rule_matches_expectations() {
    let tab = magic_circuit();
    let mut total = 0.0;
    for bits in all_bitstrings(3) {
        let by_projection = bitstring_prob_by_projection(&tab, &bits);
        let by_expectation = bitstring_prob_from_expectations(&tab, &bits);
        assert_close(by_projection, by_expectation, 1e-9);
        total += by_projection;
    }
    assert_close(total, 1.0, 1e-9);
}

#[test]
fn project_chain_rule_order_independent() {
    let tab = magic_circuit();
    for bits in all_bitstrings(3) {
        let mut t = tab.fork(Some(0));
        let mut prob = 1.0;
        for q in [2, 0, 1] {
            match t.project(q, bits[q]) {
                Ok(p) => prob *= p,
                Err(_) => {
                    prob = 0.0;
                    break;
                }
            }
        }
        assert_close(prob, bitstring_prob_from_expectations(&tab, &bits), 1e-9);
    }
}

#[test]
fn project_ghz() {
    let mut tab = TestTableau::new(3, 1e-12);
    tab.h(0);
    tab.cnot(0, 1);
    tab.cnot(1, 2);
    for bits in all_bitstrings(3) {
        let expected = if bits.iter().all(|&b| b) || bits.iter().all(|&b| !b) {
            0.5
        } else {
            0.0
        };
        assert_close(bitstring_prob_by_projection(&tab, &bits), expected, 1e-12);
    }
}

#[test]
fn project_lost_qubit_errors() {
    let mut tab = TestTableau::new(2, 1e-12);
    tab.h(0);
    tab.loss_channel(0, 1.0);
    assert_eq!(tab.project(0, false), Err(ProjectError::QubitLost(0)));
    assert!(tab.current_measurement_record().is_empty());
    // Other qubits are unaffected.
    assert_close(tab.project(1, false).unwrap(), 1.0, 1e-12);
}
