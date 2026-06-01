// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bnum::types::U256;
use itertools::Itertools;
use num::complex::Complex;
use ppvm_runtime::config::dashmap::ByteFxHashF64;
use ppvm_tableau::prelude::*;

#[test]
fn test_tableau() {
    let mut tableau: Tableau<ByteFxHashF64<1>> = Tableau::new(2);

    tableau.h(0);
    tableau.cnot(0, 1);

    assert_eq!(tableau.stabilizers()[0].to_string(), "+XX");
    assert_eq!(tableau.stabilizers()[1].to_string(), "+ZZ");

    let mut tableau: Tableau<ByteFxHashF64<1>> = Tableau::new(1);
    tableau.h(0);
    // test nonhermitian forward prop
    tableau.s(0);

    assert_eq!(tableau.stabilizers()[0].to_string(), "+Y");
    assert_eq!(tableau.destabilizers()[0].to_string(), "+Z");
}

#[test]
fn generalized_tableau() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(2, 1e-12);

    tableau.h(0);
    tableau.cnot(0, 1);
    tableau.t(0);

    assert_eq!(tableau.coefficients.len(), 2);
    let idx: Vec<_> = tableau
        .coefficients
        .clone()
        .into_iter()
        .map(|(_, i)| i)
        .sorted()
        .collect();
    assert_eq!(idx, vec![0, 1]);

    tableau.t_adj(0);

    assert_eq!(tableau.coefficients.len(), 1);

    tableau.t(0);
    tableau.t(1);

    println!("{}", tableau);

    // NOTE: since IZ|psi> = (IZ) * ZZ |psi> = ZI|psi>, we don't branch again
    assert_eq!(tableau.coefficients.len(), 2);

    let mut sorted_coefficients = tableau.coefficients.clone();
    sorted_coefficients.sort_by(|entry1, entry2| entry1.1.cmp(&entry2.1));

    const PI: f64 = std::f64::consts::PI;
    let cos_pi_8: f64 = (PI / 8.0).cos();
    let sin_pi_8: f64 = (PI / 8.0).sin();
    let expected_coefficients = [
        Complex {
            re: (PI / 4.0).cos() * (cos_pi_8 * cos_pi_8 - sin_pi_8 * sin_pi_8),
            im: (PI / 4.0).sin() * (cos_pi_8 * cos_pi_8 - sin_pi_8 * sin_pi_8),
        },
        Complex {
            re: (PI / 4.0).cos() * 2.0 * sin_pi_8 * cos_pi_8,
            im: (PI / 4.0).sin() * -2.0 * sin_pi_8 * cos_pi_8,
        },
    ];

    for ((val1, idx1), (idx2, val2)) in sorted_coefficients
        .iter()
        .zip(expected_coefficients.iter().enumerate())
    {
        assert_eq!(idx1, &(idx2 as u128));
        assert!((val1.re - val2.re).abs() < 1e-11);
        assert!((val1.im - val2.im).abs() < 1e-11);
    }

    println!("{}", tableau);
}

#[test]
fn test_generalized_tableau_phase() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(1, 1e-12);

    tableau.h(0);
    tableau.t(0);
    tableau.t(0);

    let mut sorted_coefficients = tableau.coefficients.clone();
    sorted_coefficients.sort_by(|entry1, entry2| entry1.1.cmp(&entry2.1));

    let expected_coefficients = [Complex { re: 0.5, im: 0.5 }, Complex { re: 0.5, im: -0.5 }];

    println!("{}", tableau);

    for ((val1, idx1), (idx2, val2)) in sorted_coefficients
        .iter()
        .zip(expected_coefficients.iter().enumerate())
    {
        assert_eq!(idx1, &(idx2 as u128));
        assert!((val1.re - val2.re).abs() < 1e-9);
        assert!((val1.im - val2.im).abs() < 1e-9);
    }

    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(1, 1e-12);

    tableau.h(0);
    tableau.x(0);
    tableau.t(0);

    let mut sorted_coefficients = tableau.coefficients.clone();
    sorted_coefficients.sort_by(|entry1, entry2| entry1.1.cmp(&entry2.1));

    let expected_coefficients = [
        Complex {
            re: 0.8535533905932737,
            im: 0.3535533905932738,
        },
        Complex {
            re: -0.14644660940672624,
            im: 0.3535533905932738,
        },
    ];

    for ((val1, idx1), (idx2, val2)) in sorted_coefficients
        .iter()
        .zip(expected_coefficients.iter().enumerate())
    {
        assert_eq!(idx1, &(idx2 as u128));
        assert!((val1.re - val2.re).abs() < 1e-9);
        assert!((val1.im - val2.im).abs() < 1e-9);
    }

    println!("{}", tableau);
}

#[test]
fn test_generalized_tableau_multiple_ts() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(1, 1e-12);

    tableau.h(0);

    tableau.t(0);
    tableau.t(0);
    tableau.t(0);
    tableau.t(0);

    // four T gates should be equivalent to a Z
    assert_eq!(tableau.coefficients.len(), 1);
}

#[test]
fn test_generalized_tableau_multiple_ts2() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(2, 1e-12);

    tableau.h(0);
    tableau.h(1);

    tableau.t(0);
    tableau.t(0);
    tableau.t(0);
    tableau.t(0);

    tableau.t(1);
    tableau.t(1);
    tableau.t(1);
    tableau.t(1);

    // four T gates should be equivalent to a Z
    assert_eq!(tableau.coefficients.len(), 1);
}

#[test]
fn test_generalized_tableau_multiqubit_branching() {
    let n = 18;
    let mut tableau: GeneralizedTableau<ByteFxHashF64<3>, u128> = GeneralizedTableau::new(n, 1e-12);

    for i in 0..n {
        tableau.h(i);
    }

    // make sure to branch, but watch out since we have 2 ^ t scaling
    let mut tgate_counter: u32 = 0;
    for i in (0..10).step_by(2) {
        tableau.t(i);
        tgate_counter += 1;
    }

    assert_eq!(tableau.coefficients.len(), 2_usize.pow(tgate_counter));

    // test random measurement
    let outcome = tableau.measure(0);

    // should remove a branch
    assert_eq!(tableau.coefficients.len(), 2_usize.pow(tgate_counter - 1));

    // let's move it back
    if outcome.unwrap() {
        tableau.x(0);
    }

    tableau.h(0);
    tableau.t(0);
    assert_eq!(tableau.coefficients.len(), 2_usize.pow(tgate_counter));
}

#[test]
fn test_multiqubit_ghz_state() {
    let n = 18;
    let mut tableau: GeneralizedTableau<ByteFxHashF64<3>, u128> = GeneralizedTableau::new(n, 1e-12);

    tableau.h(0);
    tableau.t(0);
    // Let's generate a GHZ state
    for i in 0..n - 1 {
        tableau.cnot(i, i + 1);
    }

    assert_eq!(tableau.coefficients.len(), 2);

    let outcome = tableau.measure(0);
    println!("{}", tableau);
    println!("{}", tableau.coefficients.len());

    for i in 0..n {
        let outcome_i = tableau.measure(i);
        assert_eq!(outcome, outcome_i)
    }
}

/// T on a Z-eigenstate (|0⟩) should not branch — Z commutes with the stabilizer +Z,
/// so the T gate simply applies a global phase to the single coefficient.
#[test]
fn test_t_on_computational_basis_no_branching() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(1, 1e-12);

    // |0⟩ is stabilized by +Z; T|0⟩ = |0⟩ (up to global phase)
    tableau.t(0);
    assert_eq!(tableau.coefficients.len(), 1, "T on |0⟩ should not branch");

    // |1⟩ is stabilized by -Z; T|1⟩ = e^{iπ/4}|1⟩
    tableau.x(0);
    tableau.t(0);
    assert_eq!(tableau.coefficients.len(), 1, "T on |1⟩ should not branch");
}

/// Verify that T†T = I: applying T then T† should leave the state unchanged.
#[test]
fn test_t_adj_cancels_t() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(2, 1e-12);

    tableau.h(0);
    tableau.h(1);
    tableau.cnot(0, 1);

    let coefficients_before = tableau.coefficients.clone();

    tableau.t(0);
    assert_eq!(
        tableau.coefficients.len(),
        2,
        "T should branch on |+⟩-like state"
    );
    tableau.t_adj(0);
    assert_eq!(
        tableau.coefficients.len(),
        1,
        "T†T should cancel back to 1 branch"
    );

    // Coefficient should match the original (up to floating point)
    assert!((tableau.coefficients[0].0.re - coefficients_before[0].0.re).abs() < 1e-10);
    assert!((tableau.coefficients[0].0.im - coefficients_before[0].0.im).abs() < 1e-10);
}

/// Clifford gates (H, X, Y, Z, S, CNOT, CZ) must never change the number of branches.
#[test]
fn test_clifford_gates_do_not_branch() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(2, 1e-12);

    // Start with a non-trivial state that has 2 branches
    tableau.h(0);
    tableau.t(0);
    assert_eq!(tableau.coefficients.len(), 2);

    // Apply every Clifford gate and verify branch count stays at 2
    tableau.x(0);
    assert_eq!(
        tableau.coefficients.len(),
        2,
        "X should not change branch count"
    );
    tableau.y(0);
    assert_eq!(
        tableau.coefficients.len(),
        2,
        "Y should not change branch count"
    );
    tableau.z(0);
    assert_eq!(
        tableau.coefficients.len(),
        2,
        "Z should not change branch count"
    );
    tableau.h(1);
    assert_eq!(
        tableau.coefficients.len(),
        2,
        "H should not change branch count"
    );
    tableau.s(0);
    assert_eq!(
        tableau.coefficients.len(),
        2,
        "S should not change branch count"
    );
    tableau.cnot(0, 1);
    assert_eq!(
        tableau.coefficients.len(),
        2,
        "CNOT should not change branch count"
    );
    tableau.cz(0, 1);
    assert_eq!(
        tableau.coefficients.len(),
        2,
        "CZ should not change branch count"
    );
}

/// After any sequence of gates, the coefficient norm should be 1.
#[test]
fn test_normalization_preserved() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(3, 1e-12);

    tableau.h(0);
    tableau.h(1);
    tableau.h(2);
    tableau.t(0);
    tableau.t(1);
    tableau.cnot(0, 1);
    tableau.t(2);
    tableau.cz(1, 2);

    let norm_sq: f64 = tableau
        .coefficients
        .clone()
        .into_iter()
        .map(|(c, _)| c.re * c.re + c.im * c.im)
        .sum();

    assert!(
        (norm_sq - 1.0).abs() < 1e-8,
        "Coefficient norm² should be 1, got {}",
        norm_sq
    );
}

/// TT produces the same quantum state as S, but with a different representation:
/// S is Clifford (modifies tableau, 1 branch), TT is non-Clifford (2 branches).
/// Verify that TT gives the known analytical coefficients for S|+⟩ = (|0⟩+i|1⟩)/√2.
#[test]
fn test_two_t_gates_coefficients() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(1, 1e-12);
    tableau.h(0);
    tableau.t(0);
    tableau.t(0);

    // TT on |+⟩ produces 2 branches
    assert_eq!(tableau.coefficients.len(), 2);

    let mut sorted = tableau.coefficients.clone();
    sorted.sort_by(|a, b| a.1.cmp(&b.1));

    // Expected: TT|+⟩ represented as two branches with these coefficients
    let expected = [Complex { re: 0.5, im: 0.5 }, Complex { re: 0.5, im: -0.5 }];

    println!("{}", tableau);

    for ((val, idx), (exp_idx, exp_val)) in sorted.iter().zip(expected.iter().enumerate()) {
        assert_eq!(*idx, exp_idx as u128);
        assert!(
            (val.re - exp_val.re).abs() < 1e-11,
            "re mismatch at index {}: {} vs {}",
            idx,
            val.re,
            exp_val.re
        );
        assert!(
            (val.im - exp_val.im).abs() < 1e-11,
            "im mismatch at index {}: {} vs {}",
            idx,
            val.im,
            exp_val.im
        );
    }

    // Norm should be 1
    let norm_sq: f64 = sorted.iter().map(|(c, _)| c.re * c.re + c.im * c.im).sum();
    assert!((norm_sq - 1.0).abs() < 1e-10);
}

/// 8 T gates on a superposition state should be equivalent to identity (T⁸ = I).
#[test]
fn test_eight_t_gates_is_identity() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(1, 1e-12);

    tableau.h(0);

    let coeff_before = tableau.coefficients.clone();

    for _ in 0..8 {
        tableau.t(0);
    }

    // T⁸ = I, so we should be back to 1 branch with the same coefficient
    assert_eq!(
        tableau.coefficients.len(),
        1,
        "T⁸ should collapse back to 1 branch"
    );
    assert!((tableau.coefficients[0].0.re - coeff_before[0].0.re).abs() < 1e-10);
    assert!((tableau.coefficients[0].0.im - coeff_before[0].0.im).abs() < 1e-10);
}

/// CZ gate: verify that CZ on |++⟩ produces the correct entangled state.
/// CZ|++⟩ = (|00⟩ + |01⟩ + |10⟩ - |11⟩)/2, which is not a stabilizer state
/// in the X basis. But from the stabilizer perspective, CZ maps
/// +XI → +XZ and +IX → +ZX, so applying CZ and then T should branch correctly.
#[test]
fn test_cz_gate_with_t() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(2, 1e-12);

    tableau.h(0);
    tableau.h(1);
    tableau.cz(0, 1);

    // The stabilizers after H H CZ should be +XZ and +ZX
    // Applying T on qubit 0: the stabilizer +XZ has X on qubit 0, so it should branch
    tableau.t(0);
    assert_eq!(
        tableau.coefficients.len(),
        2,
        "T after CZ|++⟩ should branch"
    );

    // Normalization check
    let norm_sq: f64 = tableau
        .coefficients
        .clone()
        .into_iter()
        .map(|(c, _)| c.re * c.re + c.im * c.im)
        .sum();
    assert!(
        (norm_sq - 1.0).abs() < 1e-8,
        "Norm should be preserved after CZ + T"
    );
}

/// Measurement on a generalized tableau with multiple branches should
/// preserve normalization and produce consistent outcomes on repeated measurement.
#[test]
fn test_measurement_idempotent() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(2, 1e-12);

    tableau.h(0);
    tableau.h(1);
    tableau.t(0);
    tableau.t(1);

    let outcome = tableau.measure(0);
    let branches_after_first = tableau.coefficients.len();

    // Measuring the same qubit again should be deterministic and give the same result
    let outcome2 = tableau.measure(0);
    assert_eq!(
        outcome, outcome2,
        "Repeated measurement should be deterministic"
    );
    assert_eq!(
        tableau.coefficients.len(),
        branches_after_first,
        "Second measurement should not change branch count"
    );
}

/// Verify that T gates on independent qubits produce the expected exponential branching.
#[test]
fn test_independent_t_gates_exponential_branching() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(4, 1e-12);

    for i in 0..4 {
        tableau.h(i);
    }

    // Each T on an independent qubit doubles the branch count
    for i in 0..4 {
        tableau.t(i);
        assert_eq!(
            tableau.coefficients.len(),
            2_usize.pow((i + 1) as u32),
            "After {} T gates, expected {} branches",
            i + 1,
            2_usize.pow((i + 1) as u32)
        );
    }
}

/// T|+⟩ = (|0⟩ + e^{iπ/4}|1⟩)/√2 has equal amplitudes, so Z-measurement is 50/50.
/// This verifies the generalized tableau measurement correctly computes <ψ|Z|ψ> = 0.
#[test]
fn test_t_gate_measurement_statistics() {
    let trials = 2000;
    let mut count_one = 0;

    for _ in 0..trials {
        let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> =
            GeneralizedTableau::new(1, 1e-12);
        tableau.h(0);
        tableau.t(0);

        if tableau.measure(0).unwrap() {
            count_one += 1;
        }
    }

    let prob_one = count_one as f64 / trials as f64;

    assert!(
        (prob_one - 0.5).abs() < 0.05,
        "T|+⟩ measurement should be ~50/50, got P(1)={:.3}",
        prob_one
    );
}

/// sqrt_y should implement Ry(+π/2): sqrt_y|0⟩ = |+⟩ (stabilized by +X).
/// With the bug (s, sqrt_x, s_adj order), it gives Ry(-π/2)|0⟩ = |−⟩ (stabilized by −X).
/// After H: |+⟩ → |0⟩ (measure 0), |−⟩ → |1⟩ (measure 1).
/// This circuit is purely Clifford so the measurement is deterministic.
#[test]
fn test_sqrt_y_direction() {
    use ppvm_runtime::traits::CliffordExtensions;

    // sqrt_y|0⟩ should be |+⟩
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(1, 1e-12);
    tableau.sqrt_y(0);
    tableau.h(0);
    assert!(
        !tableau.measure(0).unwrap(),
        "sqrt_y|0⟩ should be |+⟩; after H measurement must be 0"
    );

    // sqrt_y_adj|0⟩ should be |−⟩
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, u128> = GeneralizedTableau::new(1, 1e-12);
    tableau.sqrt_y_adj(0);
    tableau.h(0);
    assert!(
        tableau.measure(0).unwrap(),
        "sqrt_y_adj|0⟩ should be |−⟩; after H measurement must be 1"
    );
}

#[test]
fn test_buint_index() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<32>, U256> =
        GeneralizedTableau::new(130, 1e-12);

    tableau.h(0);

    // NOTE: would overflow at this point for u128 and more than 128 qubits
    tableau.t(1);
}
