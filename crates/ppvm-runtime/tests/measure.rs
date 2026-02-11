use num::complex::Complex64;
use ppvm_runtime::{config::dashmap::ByteFxHashF64, prelude::*};

#[test]
fn test_measure_deterministic() {
    // Test deterministic measurement: |0⟩ state
    let mut tableau: Tableau<1, ByteFxHashF64<1>> = Tableau::new();

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
    let mut tableau: Tableau<1, ByteFxHashF64<1>> = Tableau::new();

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
    let mut tableau: Tableau<1, ByteFxHashF64<1>> = Tableau::new();

    // Apply H to get |+⟩ = (|0⟩ + |1⟩)/√2
    // Stabilizer becomes +X
    tableau.h(0);
    println!("After H gate (|+⟩ state):\n{}", tableau);
    println!("Stabilizer has X: {}", tableau.stabilizers[0].word.xbits[0]);

    // Measurement is random (50/50 chance of |0⟩ or |1⟩)
    let outcome = tableau.measure(0);
    println!("\nMeasurement outcome: {}", outcome);
    println!("After measurement:\n{}", tableau);

    // After measurement, the stabilizer should be ±Z
    assert!(
        !tableau.stabilizers[0].word.xbits[0],
        "After measurement, stabilizer should not have X"
    );
    assert!(
        tableau.stabilizers[0].word.zbits[0],
        "After measurement, stabilizer should have Z"
    );

    // Phase should match the outcome: phase 0 for outcome=false, phase 2 for outcome=true
    let expected_phase = if outcome { 2 } else { 0 };
    assert_eq!(
        tableau.stabilizers[0].phase, expected_phase,
        "Stabilizer phase should match measurement outcome"
    );
}

#[test]
fn test_measure_two_qubits_bell_state() {
    // Test measurement on Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2
    let mut tableau: Tableau<2, ByteFxHashF64<1>> = Tableau::new();

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
        let mut tableau: Tableau<1, ByteFxHashF64<1>> = Tableau::new();
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
fn test_measure_generalized_tableau() {
    let mut tableau: GeneralizedTableau<2, ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1e-12);

    // Create Bell state: H on qubit 0, then CNOT(0,1)
    tableau.h(0);
    tableau.cnot(0, 1);

    tableau.t(0);
    // tableau.t(1);

    let outcome = tableau.measure(0);

    println!("Outcome: {}", outcome);
    println!("{}", tableau);

    assert!(tableau.coefficients.len() == 1);

    let tableau_outcome = tableau.tableau.measure(0);
    assert_eq!(
        tableau_outcome, outcome,
        "Tableau measurement outcome should match sampled outcome"
    );
}
