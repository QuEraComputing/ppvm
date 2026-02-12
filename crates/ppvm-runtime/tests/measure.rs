use num::complex::Complex64;
use ppvm_runtime::{config::dashmap::ByteFxHashF64, prelude::*};

#[test]
fn test_measure_deterministic() {
    // Test deterministic measurement: |0⟩ state
    let mut tableau: Tableau<ByteFxHashF64<1>> = Tableau::new(1);

    // Initial state is |0⟩, stabilizer is +Z
    println!("Initial tableau:\n{}", tableau);

    let outcome = tableau.measure(0);
    assert_eq!(
        outcome, false,
        "Measuring |0⟩ should give outcome 0 (false)"
    );

    println!("\nAfter measurement:\n{}", tableau);
}

#[test]
fn test_measure_after_x() {
    // Test deterministic measurement: |1⟩ state
    let mut tableau: Tableau<ByteFxHashF64<1>> = Tableau::new(1);

    // Apply X to get |1⟩, stabilizer becomes -Z
    tableau.x(0);
    println!("After X gate:\n{}", tableau);

    let outcome = tableau.measure(0);
    assert_eq!(outcome, true, "Measuring |1⟩ should give outcome 1 (true)");

    println!("\nAfter measurement:\n{}", tableau);
}

#[test]
fn test_measure_after_hadamard() {
    // Test random measurement: |+⟩ state
    let mut tableau: Tableau<ByteFxHashF64<1>> = Tableau::new(1);

    // Apply H to get |+⟩ = (|0⟩ + |1⟩)/√2
    // Stabilizer becomes +X
    tableau.h(0);
    println!("After H gate (|+⟩ state):\n{}", tableau);
    println!(
        "Stabilizer has X: {}",
        tableau.stabilizers()[0].word.xbits[0]
    );

    // Measurement is random (50/50 chance of |0⟩ or |1⟩)
    let outcome = tableau.measure(0);
    println!("\nMeasurement outcome: {}", outcome);
    println!("After measurement:\n{}", tableau);

    // After measurement, the stabilizer should be ±Z
    assert!(
        !tableau.stabilizers()[0].word.xbits[0],
        "After measurement, stabilizer should not have X"
    );
    assert!(
        tableau.stabilizers()[0].word.zbits[0],
        "After measurement, stabilizer should have Z"
    );

    // Phase should match the outcome: phase 0 for outcome=false, phase 2 for outcome=true
    let expected_phase = if outcome { 2 } else { 0 };
    assert_eq!(
        tableau.stabilizers()[0].phase,
        expected_phase,
        "Stabilizer phase should match measurement outcome"
    );
}

#[test]
fn test_measure_two_qubits_bell_state() {
    // Test measurement on Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2
    let mut tableau: Tableau<ByteFxHashF64<1>> = Tableau::new(2);

    // Create Bell state: H on qubit 0, then CNOT(0,1)
    tableau.h(0);
    tableau.cnot(0, 1);
    println!("Bell state tableau:\n{}", tableau);

    // Measuring qubit 0 is random (50/50 for |0⟩ or |1⟩)
    let outcome0 = tableau.measure(0);
    println!("\nFirst measurement outcome: {}", outcome0);
    println!("After first measurement:\n{}", tableau);

    // After measuring qubit 0, qubit 1 should be in the same state (deterministic)
    // This is because of entanglement in the Bell state
    let outcome1 = tableau.measure(1);
    println!("\nSecond measurement outcome: {}", outcome1);
    println!("After second measurement:\n{}", tableau);

    // Both measurements should give the same result (perfect correlation in Bell state)
    assert_eq!(
        outcome0, outcome1,
        "Bell state measurements should be perfectly correlated"
    );
}

#[test]
fn test_measure_statistics() {
    // Test that random measurements give approximately 50/50 distribution
    let trials = 1000;
    let mut count_zero = 0;
    let mut count_one = 0;

    for _ in 0..trials {
        let mut tableau: Tableau<ByteFxHashF64<1>> = Tableau::new(1);
        tableau.h(0); // Create |+⟩ state

        let outcome = tableau.measure(0);
        if outcome {
            count_one += 1;
        } else {
            count_zero += 1;
        }
    }

    println!("Statistics over {} trials:", trials);
    println!(
        "  Outcome 0: {} ({:.1}%)",
        count_zero,
        100.0 * count_zero as f64 / trials as f64
    );
    println!(
        "  Outcome 1: {} ({:.1}%)",
        count_one,
        100.0 * count_one as f64 / trials as f64
    );

    // Check that distribution is roughly 50/50 (within 3 sigma for binomial distribution)
    // For 1000 trials, standard deviation is sqrt(1000 * 0.5 * 0.5) ≈ 15.8
    // 3 sigma ≈ 47.4, so we expect outcomes in range [450, 550]
    assert!(
        count_zero >= 400 && count_zero <= 600,
        "Measurement statistics should be approximately 50/50 (got {} zeros out of {})",
        count_zero,
        trials
    );
    assert!(
        count_one >= 400 && count_one <= 600,
        "Measurement statistics should be approximately 50/50 (got {} ones out of {})",
        count_one,
        trials
    );
}

#[test]
fn test_measure_generalized_tableau_bell() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(2, 1e-12);

    // Create Bell state: H on qubit 0, then CNOT(0,1)
    tableau.h(0);
    tableau.cnot(0, 1);

    tableau.t(0);
    // tableau.t(1);

    let outcome = tableau.measure(0);
    println!("{}", tableau);

    println!("Outcome: {}", outcome);
    assert!(tableau.coefficients.len() == 1);

    let tableau_outcome = tableau.tableau.measure(0);
    assert_eq!(
        tableau_outcome, outcome,
        "Tableau measurement outcome should match sampled outcome"
    );

    let outcome2 = tableau.measure(1);

    assert_eq!(
        outcome2, outcome,
        "Bell state measurement must be consistents"
    )
}

#[test]
fn test_measure_generalized_tableau_deterministic() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1, 1e-12);

    // Create |+⟩ state
    let outcome = tableau.measure(0);
    assert_eq!(tableau.coefficients.len(), 1);
    assert_eq!(outcome, false);
    assert!((tableau.coefficients[0].0 - 1.0).re.abs() < 1e-10);
    assert!(tableau.coefficients[0].0.im.abs() < 1e-10);

    tableau.x(0);
    let outcome = tableau.measure(0);
    assert_eq!(tableau.coefficients.len(), 1);
    assert_eq!(outcome, true);
    assert!((tableau.coefficients[0].0 - 1.0).re.abs() < 1e-10);
    assert!(tableau.coefficients[0].0.im.abs() < 1e-10);
}

#[test]
fn test_measure_generalized_random() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1, 1e-12);

    // Create |+⟩ state
    tableau.h(0);
    tableau.t(0);
    println!("Original tableau: {}", tableau);
    let outcome = tableau.measure(0);
    assert_eq!(tableau.coefficients.len(), 1);
    let r = tableau.coefficients[0].0.re.abs();
    let i = tableau.coefficients[0].0.im.abs();
    assert!(((r * r + i * i).sqrt() - 1.0) < 1e-10);
    println!("{}", tableau);
    println!("{}", outcome);
}

#[test]
fn test_measure_generalized_tableau_statistics() {
    let trials = 1000;
    let mut count_zero = 0;
    let mut count_one = 0;

    for _ in 0..trials {
        let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
            GeneralizedTableau::new(2, 1e-12);

        // Create Bell state: H on qubit 0, then CNOT(0,1)
        tableau.h(0);
        tableau.cnot(0, 1);

        tableau.t(0);
        tableau.t(1);

        let outcome = tableau.measure(0);
        if outcome {
            count_one += 1;
        } else {
            count_zero += 1;
        }

        assert_eq!(tableau.coefficients.len(), 1);
        let outcome2 = tableau.measure(1);
        assert_eq!(outcome, outcome2, "Bell measurements must be consistent");
    }

    println!("Statistics over {} trials:", trials);
    println!(
        "  Outcome 0: {} ({:.1}%)",
        count_zero,
        100.0 * count_zero as f64 / trials as f64
    );
    println!(
        "  Outcome 1: {} ({:.1}%)",
        count_one,
        100.0 * count_one as f64 / trials as f64
    );

    assert!(
        count_zero >= 400 && count_zero <= 600,
        "Measurement statistics should be approximately 50/50 (got {} zeros out of {})",
        count_zero,
        trials
    );
    assert!(
        count_one >= 400 && count_one <= 600,
        "Measurement statistics should be approximately 50/50 (got {} ones out of {})",
        count_one,
        trials
    );
}

/// Coefficients must be normalized (Σ|c|² = 1) after measurement on a multi-branch state.
#[test]
fn test_measure_generalized_normalization() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(3, 1e-12);

    tableau.h(0);
    tableau.h(1);
    tableau.h(2);
    tableau.t(0);
    tableau.t(1);
    tableau.t(2);

    // 8 branches before measurement
    assert_eq!(tableau.coefficients.len(), 8);

    tableau.measure(0);

    let norm_sq: f64 = tableau
        .coefficients
        .clone()
        .into_iter()
        .map(|(c, _)| c.re * c.re + c.im * c.im)
        .sum();
    assert!(
        (norm_sq - 1.0).abs() < 1e-8,
        "Norm² should be 1 after measurement, got {}",
        norm_sq
    );

    tableau.measure(1);

    let norm_sq: f64 = tableau
        .coefficients
        .clone()
        .into_iter()
        .map(|(c, _)| c.re * c.re + c.im * c.im)
        .sum();
    assert!(
        (norm_sq - 1.0).abs() < 1e-8,
        "Norm² should be 1 after second measurement, got {}",
        norm_sq
    );
}

/// T on |0⟩ doesn't change the state (only global phase), so measurement is still deterministic.
#[test]
fn test_measure_generalized_deterministic_with_t() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1, 1e-12);

    // T|0⟩ = |0⟩ (no branching, Z eigenstate)
    tableau.t(0);
    assert_eq!(tableau.coefficients.len(), 1);

    let outcome = tableau.measure(0);
    assert_eq!(outcome, false, "T|0⟩ should measure as 0");
    assert_eq!(tableau.coefficients.len(), 1);

    // T|1⟩ = e^{iπ/4}|1⟩ (no branching, Z eigenstate)
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1, 1e-12);
    tableau.x(0);
    tableau.t(0);
    assert_eq!(tableau.coefficients.len(), 1);

    let outcome = tableau.measure(0);
    assert_eq!(outcome, true, "T|1⟩ should measure as 1");
    assert_eq!(tableau.coefficients.len(), 1);
}

/// Measurement halves the branch count when the measured qubit caused branching.
#[test]
fn test_measure_reduces_branches() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(3, 1e-12);

    tableau.h(0);
    tableau.h(1);
    tableau.h(2);
    tableau.t(0);
    tableau.t(1);
    tableau.t(2);
    assert_eq!(tableau.coefficients.len(), 8);

    tableau.measure(0);
    assert_eq!(
        tableau.coefficients.len(),
        4,
        "Measuring 1 of 3 T-branched qubits: 8 → 4"
    );

    tableau.measure(1);
    assert_eq!(tableau.coefficients.len(), 2, "Measuring 2nd: 4 → 2");

    tableau.measure(2);
    assert_eq!(tableau.coefficients.len(), 1, "Measuring 3rd: 2 → 1");
}

/// On a product state with independent T gates, measuring one qubit
/// should not affect the other qubit's measurement statistics.
#[test]
fn test_measure_product_state_independence() {
    let trials = 1000;
    let mut count_q1_one = 0;

    for _ in 0..trials {
        let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
            GeneralizedTableau::new(2, 1e-12);

        tableau.h(0);
        tableau.h(1);
        tableau.t(0);
        tableau.t(1);

        // Measure qubit 0 first (discard result)
        tableau.measure(0);

        // Qubit 1 should still be 50/50 (T|+⟩ has equal amplitudes)
        if tableau.measure(1) {
            count_q1_one += 1;
        }
    }

    let prob = count_q1_one as f64 / trials as f64;
    assert!(
        (prob - 0.5).abs() < 0.06,
        "Qubit 1 should be ~50/50 regardless of qubit 0 outcome, got P(1)={:.3}",
        prob
    );
}

/// Measuring all qubits of a GHZ-like state (with T) should give perfectly correlated outcomes.
#[test]
fn test_measure_generalized_ghz_correlation() {
    for _ in 0..50 {
        let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
            GeneralizedTableau::new(4, 1e-12);

        tableau.h(0);
        tableau.t(0);
        for i in 0..3 {
            tableau.cnot(i, i + 1);
        }

        let first = tableau.measure(0);
        for i in 1..4 {
            let outcome = tableau.measure(i);
            assert_eq!(
                outcome, first,
                "GHZ qubit {} should match qubit 0 (trial outcome={})",
                i, first
            );
        }
    }
}

/// After measurement, re-measuring the same qubit must always return the same outcome
/// and leave coefficients unchanged.
#[test]
fn test_measure_generalized_idempotent() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(2, 1e-12);

    tableau.h(0);
    tableau.h(1);
    tableau.t(0);
    tableau.t(1);

    let outcome = tableau.measure(0);
    let coeffs_after: Vec<_> = tableau.coefficients.clone().into_iter().collect();

    // Re-measure same qubit multiple times
    for _ in 0..5 {
        let repeated = tableau.measure(0);
        assert_eq!(
            repeated, outcome,
            "Repeated measurement must be deterministic"
        );

        let coeffs_now: Vec<_> = tableau.coefficients.clone().into_iter().collect();
        assert_eq!(coeffs_now.len(), coeffs_after.len());
        for ((c1, i1), (c2, i2)) in coeffs_after.iter().zip(coeffs_now.iter()) {
            assert_eq!(i1, i2);
            assert!((c1.re - c2.re).abs() < 1e-12);
            assert!((c1.im - c2.im).abs() < 1e-12);
        }
    }
}

/// Measure on a 4-qubit entangled state with T gates interspersed.
/// Verifies that measurement collapses branches and maintains valid state.
#[test]
fn test_measure_generalized_entangled_chain() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(4, 1e-12);

    tableau.h(0);
    tableau.t(0);
    tableau.cnot(0, 1);
    tableau.h(2);
    tableau.t(2);
    tableau.cnot(2, 3);
    tableau.cnot(1, 2);

    let branches_before = tableau.coefficients.len();
    assert!(branches_before > 1, "State should have multiple branches");

    // Measure qubits one by one; each should reduce or maintain branch count
    let mut prev_branches = branches_before;
    for i in 0..4 {
        tableau.measure(i);
        assert!(
            tableau.coefficients.len() <= prev_branches,
            "Branch count should not increase after measurement"
        );

        let norm_sq: f64 = tableau
            .coefficients
            .clone()
            .into_iter()
            .map(|(c, _)| c.re * c.re + c.im * c.im)
            .sum();
        assert!(
            (norm_sq - 1.0).abs() < 1e-8,
            "Norm² should be 1 after measuring qubit {}, got {}",
            i,
            norm_sq
        );

        prev_branches = tableau.coefficients.len();
    }

    // After measuring all qubits, should have exactly 1 branch
    assert_eq!(tableau.coefficients.len(), 1);
}

/// Verify that the generalized tableau and its inner tableau agree on
/// deterministic measurement outcomes.
#[test]
fn test_measure_generalized_agrees_with_inner_tableau() {
    for _ in 0..20 {
        let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
            GeneralizedTableau::new(2, 1e-12);

        tableau.h(0);
        tableau.cnot(0, 1);
        tableau.t(0);

        let outcome0 = tableau.measure(0);

        // The inner tableau should now deterministically agree
        let inner_outcome0 = tableau.tableau.measure(0);
        assert_eq!(
            outcome0, inner_outcome0,
            "Inner tableau must agree with generalized measurement"
        );

        // Qubit 1 should be correlated (Bell state)
        let outcome1 = tableau.measure(1);
        assert_eq!(outcome0, outcome1, "Bell state qubits must be correlated");
    }
}

#[test]
fn test_measure_generalized_tableau_t_gate_deterministic() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1, 1e-12);

    tableau.h(0);

    // Two T gates + 3S gates should rotate around Z
    tableau.t(0);
    tableau.t(0);
    tableau.s(0);
    tableau.s(0);
    tableau.s(0);

    // Another H and we should be back to |0⟩ (deterministic)
    tableau.h(0);

    println!("Tableau before measurement:\n{}", tableau);
    let outcome = tableau.measure(0);
    assert_eq!(
        outcome, false,
        "State should be |0⟩ after T and S rotations"
    );
}

#[test]
fn test_measure_generalized_tableau_t_gate_random() {
    let mut tableau: GeneralizedTableau<ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1, 1e-12);

    tableau.h(0);

    // Two T gates + 3S gates should rotate around Z
    tableau.t(0);
    tableau.t(0);
    tableau.s(0);
    tableau.s(0);
    tableau.s(0);

    let trials = 1000;
    let mut count_q1_one = 0;

    for _ in 0..trials {
        let mut copy = tableau.clone();
        if copy.measure(0) {
            count_q1_one += 1;
        }
    }

    let probability = count_q1_one as f64 / trials as f64;
    println!("Probability of measuring |1⟩ on qubit 0: {}", probability);
    assert!(
        (probability - 0.5).abs() < 0.06,
        "Measurement should be approximately 50/50 after T and S rotations, got P(1)={:.3}",
        probability
    );
}
