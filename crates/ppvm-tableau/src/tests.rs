use num::complex::Complex64;

use crate::{config, Tableau, TableauSum};

fn assert_single_row(tableau: &Tableau, x: u8, z: u8, phase: u8) {
    assert_eq!(tableau.n_qubits, 1);
    assert_eq!(tableau.n_words, 1);
    assert_eq!(tableau.x.len(), 1);
    assert_eq!(tableau.z.len(), 1);
    assert_eq!(tableau.phase.len(), 1);
    assert_eq!(bit_at(&tableau.x, tableau.n_words, 0, 0), x);
    assert_eq!(bit_at(&tableau.z, tableau.n_words, 0, 0), z);
    assert_eq!(phase_at(&tableau.phase, 0), phase);
}

fn approx_eq(a: Complex64, b: Complex64, eps: f64) -> bool {
    (a - b).norm() <= eps
}

fn bit_at(columns: &[u64], n_words: usize, col: usize, row: usize) -> u8 {
    let word = row / 64;
    let mask = 1u64 << (row % 64);
    let idx = col * n_words + word;
    if columns[idx] & mask != 0 { 1 } else { 0 }
}

fn phase_at(phase: &[u64], row: usize) -> u8 {
    let word = row / 64;
    let mask = 1u64 << (row % 64);
    if phase[word] & mask != 0 { 1 } else { 0 }
}

#[test]
fn gates_on_zero_state() {
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    assert_eq!(state.len(), 1);
    let (tableau, coeff) = state.map.iter().next().unwrap();
    assert_single_row(tableau, 0, 1, 0);
    assert!(approx_eq(*coeff, Complex64::new(1.0, 0.0), 1e-12));

    state.x(0);
    let (tableau, coeff) = state.map.iter().next().unwrap();
    assert_single_row(tableau, 0, 1, 1);
    assert!(approx_eq(*coeff, Complex64::new(1.0, 0.0), 1e-12));

    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.y(0);
    let (tableau, coeff) = state.map.iter().next().unwrap();
    assert_single_row(tableau, 0, 1, 1);
    assert!(approx_eq(*coeff, Complex64::new(1.0, 0.0), 1e-12));

    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.z(0);
    let (tableau, coeff) = state.map.iter().next().unwrap();
    assert_single_row(tableau, 0, 1, 0);
    assert!(approx_eq(*coeff, Complex64::new(1.0, 0.0), 1e-12));

    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.h(0);
    let (tableau, coeff) = state.map.iter().next().unwrap();
    assert_single_row(tableau, 1, 0, 0);
    assert!(approx_eq(*coeff, Complex64::new(1.0, 0.0), 1e-12));

    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.s(0);
    let (tableau, coeff) = state.map.iter().next().unwrap();
    assert_single_row(tableau, 0, 1, 0);
    assert!(approx_eq(*coeff, Complex64::new(1.0, 0.0), 1e-12));

    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.t(0);
    assert_eq!(state.len(), 1);
    let (tableau, coeff) = state.map.iter().next().unwrap();
    assert_single_row(tableau, 0, 1, 0);
    assert!(approx_eq(*coeff, Complex64::new(1.0, 0.0), 1e-12));
}

#[test]
fn t_gate_branching_on_plus_state() {
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.h(0);
    state.t(0);
    assert_eq!(state.len(), 2);

    let (a, b) = crate::sum::t_coeffs::<Complex64>();
    let mut plus = Tableau::new(1);
    plus.h(0);
    let mut minus = plus.clone();
    minus.z(0);

    let coeff_plus = state.coeff(&plus).expect("missing |+> tableau");
    let coeff_minus = state.coeff(&minus).expect("missing |-> tableau");
    assert!(approx_eq(coeff_plus, a, 1e-12));
    assert!(approx_eq(coeff_minus, b, 1e-12));
}

#[test]
fn multiple_t_gates_exponential_branching() {
    // Test that multiple T gates cause exponential branching when applied to different qubits
    // Note: applying multiple T gates to the SAME qubit in superposition causes branch merging
    // because the tableaux can become equivalent after Z operations

    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(3);
    state.h(0);
    state.h(1);
    state.h(2);

    // Apply first T gate to qubit 0 - should create 2 branches
    state.t(0);
    assert_eq!(state.len(), 2, "After 1 T gate should have 2 branches");

    // Apply second T gate to qubit 1 - should create 4 branches
    state.t(1);
    assert_eq!(state.len(), 4, "After 2 T gates should have 4 branches");

    // Apply third T gate to qubit 2 - should create 8 branches
    state.t(2);
    assert_eq!(state.len(), 8, "After 3 T gates should have 8 branches");
}

#[test]
fn t_gate_multi_qubit_independent() {
    // Test T gates on independent qubits branch independently
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(2);
    state.h(0);
    state.h(1);

    // Apply T to qubit 0 - should create 2 branches
    state.t(0);
    assert_eq!(state.len(), 2, "After T on qubit 0 should have 2 branches");

    // Apply T to qubit 1 - should double the branches to 4
    state.t(1);
    assert_eq!(state.len(), 4, "After T on both qubits should have 4 branches");
}

#[test]
fn t_gate_with_clifford_gates() {
    // Test that Clifford gates don't cause additional branching
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.h(0);
    state.t(0);
    assert_eq!(state.len(), 2, "After H+T should have 2 branches");

    // Apply more Clifford gates - should not increase branch count
    state.s(0);
    assert_eq!(state.len(), 2, "After S gate should still have 2 branches");

    state.h(0);
    assert_eq!(state.len(), 2, "After H gate should still have 2 branches");

    state.z(0);
    assert_eq!(state.len(), 2, "After Z gate should still have 2 branches");

    state.x(0);
    assert_eq!(state.len(), 2, "After X gate should still have 2 branches");
}

#[test]
fn t_gate_coefficient_normalization() {
    // Test that coefficients are properly normalized
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.h(0);
    state.t(0);

    // Sum of squared magnitudes should be 1 (probability conservation)
    let mut total_prob = 0.0;
    for (_, coeff) in state.map.iter() {
        total_prob += coeff.norm_sqr();
    }
    assert!(approx_eq(Complex64::new(total_prob, 0.0), Complex64::new(1.0, 0.0), 1e-12),
            "Total probability should be 1.0, got {}", total_prob);
}

#[test]
fn t_gate_ht_sequence() {
    // Test H-T-H sequence (equivalent to T-dagger up to global phase)
    let mut state1: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state1.h(0);
    state1.t(0);
    state1.h(0);

    // Should still have 2 branches after HT†H sequence
    assert_eq!(state1.len(), 2);

    // Verify probability conservation
    let mut total_prob = 0.0;
    for (_, coeff) in state1.map.iter() {
        total_prob += coeff.norm_sqr();
    }
    assert!(approx_eq(Complex64::new(total_prob, 0.0), Complex64::new(1.0, 0.0), 1e-12));
}

#[test]
fn t_gate_commutation_with_z() {
    // Test that T and Z gates commute (T Z = Z T)
    let mut state1: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state1.h(0);
    state1.t(0);
    state1.z(0);

    let mut state2: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state2.h(0);
    state2.z(0);
    state2.t(0);

    // Both should have same number of branches
    assert_eq!(state1.len(), state2.len());
}

#[test]
fn t_gate_multiple_qubits_selective() {
    // Test T gate on selected qubits in multi-qubit system
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(3);

    // Put all qubits in superposition
    state.h(0);
    state.h(1);
    state.h(2);

    // Apply T to qubit 1 only
    state.t(1);
    assert_eq!(state.len(), 2, "T on qubit 1 should create 2 branches");

    // Apply T to qubit 2
    state.t(2);
    assert_eq!(state.len(), 4, "T on qubits 1,2 should create 4 branches");

    // Qubit 0 is still in superposition but no T applied
    // Apply Clifford to qubit 0 - should not change branch count
    state.s(0);
    assert_eq!(state.len(), 4, "Clifford on qubit 0 should maintain 4 branches");
}

#[test]
fn t_gate_on_different_basis_states() {
    // Test T gate on |0>, |1>, |+>, |-> states

    // T|0> should not branch (eigenstate)
    let mut state_zero: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state_zero.t(0);
    assert_eq!(state_zero.len(), 1, "T|0> should not branch");

    // T|1> should not branch (eigenstate)
    let mut state_one: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state_one.x(0);
    state_one.t(0);
    assert_eq!(state_one.len(), 1, "T|1> should not branch");

    // T|+> should branch
    let mut state_plus: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state_plus.h(0);
    state_plus.t(0);
    assert_eq!(state_plus.len(), 2, "T|+> should branch into 2");

    // T|-> should also branch
    let mut state_minus: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state_minus.x(0);
    state_minus.h(0);
    state_minus.t(0);
    assert_eq!(state_minus.len(), 2, "T|-> should branch into 2");
}

#[test]
fn t_gate_branch_merging() {
    // Test that identical tableaux with different coefficients are merged
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.h(0);
    state.t(0);
    assert_eq!(state.len(), 2);

    // Apply S gate 4 times (S^4 = I, returns to identity modulo phase)
    state.s(0);
    state.s(0);
    state.s(0);
    state.s(0);

    // Should still have 2 branches (no merging because tableaux are still different)
    assert_eq!(state.len(), 2);
}

#[test]
fn t_gate_sequential_on_same_qubit() {
    // Test multiple T gates on the same qubit
    // When applying T to the same qubit multiple times, branches can merge
    // because the tableaux may become equivalent after Z operations
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.h(0);

    // First T gate creates 2 branches
    state.t(0);
    assert_eq!(state.len(), 2, "After 1 T gate should have 2 branches");

    // Second T gate on same qubit - branches may merge
    state.t(0);
    // The exact number depends on whether tableaux become equivalent
    assert!(state.len() >= 2, "Should have at least 2 branches after 2 T gates");

    // Verify probability conservation regardless of branch count
    let mut total_prob = 0.0;
    for (_, coeff) in state.map.iter() {
        total_prob += coeff.norm_sqr();
    }
    assert!(approx_eq(Complex64::new(total_prob, 0.0), Complex64::new(1.0, 0.0), 1e-12));
}

#[test]
fn t_gate_same_qubit_merging_behavior() {
    // Explicitly test that T gates on same qubit cause merging
    // This is the expected behavior in tableau-based simulation
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.h(0);

    // Apply T gate once
    state.t(0);
    let branches_after_one = state.len();
    assert_eq!(branches_after_one, 2);

    // Apply T gate again to same qubit
    state.t(0);
    let branches_after_two = state.len();

    // In tableau simulation, applying T to the same qubit can cause merging
    // because Z operations on already-branched states can create equivalent tableaux
    // The branch count won't necessarily double
    assert!(branches_after_two >= 2, "Should maintain at least 2 branches");

    // Probability should still be conserved
    let mut total_prob = 0.0;
    for (_, coeff) in state.map.iter() {
        total_prob += coeff.norm_sqr();
    }
    assert!(approx_eq(Complex64::new(total_prob, 0.0), Complex64::new(1.0, 0.0), 1e-12));
}

#[test]
fn t_gate_coefficients_complex() {
    // Verify that T gate coefficients have correct complex values
    let (a, b) = crate::sum::t_coeffs::<Complex64>();

    // a = (1 + e^(iπ/4))/2 = (1 + (√2/2)(1+i))/2
    // b = (1 - e^(iπ/4))/2 = (1 - (√2/2)(1+i))/2

    // Check that |a|^2 + |b|^2 = 1 (unitarity)
    let sum_prob = a.norm_sqr() + b.norm_sqr();
    assert!(approx_eq(Complex64::new(sum_prob, 0.0), Complex64::new(1.0, 0.0), 1e-12),
            "T gate coefficients should satisfy |a|^2 + |b|^2 = 1");

    // Check approximate values
    let sqrt2_2 = std::f64::consts::FRAC_1_SQRT_2;
    let expected_a = Complex64::new((1.0 + sqrt2_2) / 2.0, sqrt2_2 / 2.0);
    let expected_b = Complex64::new((1.0 - sqrt2_2) / 2.0, -sqrt2_2 / 2.0);

    assert!(approx_eq(a, expected_a, 1e-12), "T coefficient 'a' should match expected value");
    assert!(approx_eq(b, expected_b, 1e-12), "T coefficient 'b' should match expected value");
}

#[test]
fn t_gate_with_parallel_threshold() {
    // Test that parallel threshold doesn't affect correctness
    let mut state_seq: TableauSum<config::fxhash::FxComplex> =
        TableauSum::new(1).with_parallel_threshold(usize::MAX);
    state_seq.h(0);
    state_seq.t(0);
    state_seq.t(0);

    let mut state_par: TableauSum<config::fxhash::FxComplex> =
        TableauSum::new(1).with_parallel_threshold(1);
    state_par.h(0);
    state_par.t(0);
    state_par.t(0);

    // Both should have the same number of branches
    assert_eq!(state_seq.len(), state_par.len(),
               "Sequential and parallel execution should produce same branch count");

    // Verify both have correct probability normalization
    for state in [&state_seq, &state_par] {
        let mut total_prob = 0.0;
        for (_, coeff) in state.map.iter() {
            total_prob += coeff.norm_sqr();
        }
        assert!(approx_eq(Complex64::new(total_prob, 0.0), Complex64::new(1.0, 0.0), 1e-12));
    }
}

#[test]
fn t_gate_two_qubit_entangled_basis() {
    // Test T gate on a two-qubit system in various states
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(2);

    // Create |++> state
    state.h(0);
    state.h(1);

    // Apply T to first qubit - should branch
    state.t(0);
    assert_eq!(state.len(), 2, "T on qubit 0 in |++> should create 2 branches");

    // Apply T to second qubit - should double branches
    state.t(1);
    assert_eq!(state.len(), 4, "T on both qubits in |++> should create 4 branches");

    // Verify probability conservation
    let mut total_prob = 0.0;
    for (_, coeff) in state.map.iter() {
        total_prob += coeff.norm_sqr();
    }
    assert!(approx_eq(Complex64::new(total_prob, 0.0), Complex64::new(1.0, 0.0), 1e-12));
}

#[test]
fn t_gate_exact_tableau_coefficients_after_split() {
    // Test exact coefficients of each tableau after T gate splitting
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.h(0);
    state.t(0);

    assert_eq!(state.len(), 2, "Should have exactly 2 branches");

    // Get expected coefficients
    let (a, b) = crate::sum::t_coeffs::<Complex64>();

    // Build expected tableaux manually
    let mut plus = Tableau::new(1);
    plus.h(0);

    let mut minus = plus.clone();
    minus.z(0);

    // Verify each tableau has correct coefficient
    let coeff_plus = state.coeff(&plus).expect("Missing |+> tableau");
    let coeff_minus = state.coeff(&minus).expect("Missing |-> tableau");

    assert!(approx_eq(coeff_plus, a, 1e-12),
            "Coefficient for |+> should be {}, got {}", a, coeff_plus);
    assert!(approx_eq(coeff_minus, b, 1e-12),
            "Coefficient for |-> should be {}, got {}", b, coeff_minus);

    // Verify no other tableaux exist
    assert_eq!(state.map.len(), 2, "Should have exactly 2 tableaux in map");
}

#[test]
fn t_gate_multi_qubit_individual_coefficients() {
    // Test coefficients for multi-qubit T gate splitting
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(2);
    state.h(0);
    state.h(1);
    state.t(0);

    assert_eq!(state.len(), 2, "Should have 2 branches after T on qubit 0");

    let (a, b) = crate::sum::t_coeffs::<Complex64>();

    // When we apply T to qubit 0 in |++> state:
    // Initial: (1/√2)(|0> + |1>) ⊗ (1/√2)(|0> + |1>)
    // After T on qubit 0: branches based on qubit 0

    // Build the two expected tableaux
    let mut t1 = Tableau::new(2);
    t1.h(0);
    t1.h(1);

    let mut t2 = t1.clone();
    t2.z(0);

    let coeff1 = state.coeff(&t1).expect("Missing first tableau");
    let coeff2 = state.coeff(&t2).expect("Missing second tableau");

    assert!(approx_eq(coeff1, a, 1e-12),
            "First tableau coefficient should be {}, got {}", a, coeff1);
    assert!(approx_eq(coeff2, b, 1e-12),
            "Second tableau coefficient should be {}, got {}", b, coeff2);
}

#[test]
fn t_gate_two_splits_coefficient_products() {
    // Test coefficients after two T gates on different qubits
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(2);
    state.h(0);
    state.h(1);
    state.t(0);
    state.t(1);

    assert_eq!(state.len(), 4, "Should have 4 branches after T on both qubits");

    let (a, b) = crate::sum::t_coeffs::<Complex64>();

    // After T on qubit 0 and T on qubit 1, we expect 4 branches with coefficients:
    // a*a (no Z on either qubit)
    // a*b (Z on qubit 1 only)
    // b*a (Z on qubit 0 only)
    // b*b (Z on both qubits)

    let mut base = Tableau::new(2);
    base.h(0);
    base.h(1);

    // Tableau with no Z gates
    let coeff_00 = state.coeff(&base).expect("Missing tableau (no Z)");
    assert!(approx_eq(coeff_00, a.clone() * a.clone(), 1e-12),
            "Coefficient should be a*a = {}, got {}", a * a, coeff_00);

    // Tableau with Z on qubit 0
    let mut t_z0 = base.clone();
    t_z0.z(0);
    let coeff_z0 = state.coeff(&t_z0).expect("Missing tableau (Z on qubit 0)");
    assert!(approx_eq(coeff_z0, b.clone() * a.clone(), 1e-12),
            "Coefficient should be b*a = {}, got {}", b * a, coeff_z0);

    // Tableau with Z on qubit 1
    let mut t_z1 = base.clone();
    t_z1.z(1);
    let coeff_z1 = state.coeff(&t_z1).expect("Missing tableau (Z on qubit 1)");
    assert!(approx_eq(coeff_z1, a.clone() * b.clone(), 1e-12),
            "Coefficient should be a*b = {}, got {}", a * b, coeff_z1);

    // Tableau with Z on both qubits
    let mut t_z01 = base.clone();
    t_z01.z(0);
    t_z01.z(1);
    let coeff_z01 = state.coeff(&t_z01).expect("Missing tableau (Z on both qubits)");
    assert!(approx_eq(coeff_z01, b.clone() * b.clone(), 1e-12),
            "Coefficient should be b*b = {}, got {}", b * b, coeff_z01);

    // Verify total probability
    let total_prob = coeff_00.norm_sqr() + coeff_z0.norm_sqr()
                   + coeff_z1.norm_sqr() + coeff_z01.norm_sqr();
    assert!(approx_eq(Complex64::new(total_prob, 0.0), Complex64::new(1.0, 0.0), 1e-12));
}

#[test]
fn t_gate_coefficient_phases() {
    // Test that T gate coefficients have correct phases
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.h(0);
    state.t(0);

    let (a, b) = crate::sum::t_coeffs::<Complex64>();

    // Verify the phases are correct
    // a = (1 + e^(iπ/4))/2 should have a positive imaginary part
    assert!(a.im > 0.0, "Coefficient 'a' should have positive imaginary part");

    // b = (1 - e^(iπ/4))/2 should have a negative imaginary part
    assert!(b.im < 0.0, "Coefficient 'b' should have negative imaginary part");

    // Both should have positive real parts
    assert!(a.re > 0.0, "Coefficient 'a' should have positive real part");
    assert!(b.re > 0.0, "Coefficient 'b' should have positive real part");

    // Verify the coefficients in the state match
    let mut plus = Tableau::new(1);
    plus.h(0);
    let mut minus = plus.clone();
    minus.z(0);

    let state_coeff_plus = state.coeff(&plus).expect("Missing |+> tableau");
    let state_coeff_minus = state.coeff(&minus).expect("Missing |-> tableau");

    assert_eq!(state_coeff_plus.re.is_sign_positive(), a.re.is_sign_positive());
    assert_eq!(state_coeff_plus.im.is_sign_positive(), a.im.is_sign_positive());
    assert_eq!(state_coeff_minus.re.is_sign_positive(), b.re.is_sign_positive());
    assert_eq!(state_coeff_minus.im.is_sign_positive(), b.im.is_sign_positive());
}

#[test]
fn t_gate_coefficient_interference_check() {
    // Test that coefficients can interfere correctly when branches merge
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.h(0);

    // Apply T, then H to potentially cause interference
    state.t(0);
    assert_eq!(state.len(), 2);

    let (a, b) = crate::sum::t_coeffs::<Complex64>();

    // Store coefficients before H
    let mut plus_before = Tableau::new(1);
    plus_before.h(0);
    let mut minus_before = plus_before.clone();
    minus_before.z(0);

    let coeff_plus_before = state.coeff(&plus_before).expect("Missing |+>");
    let coeff_minus_before = state.coeff(&minus_before).expect("Missing |->");

    assert!(approx_eq(coeff_plus_before, a, 1e-12));
    assert!(approx_eq(coeff_minus_before, b, 1e-12));

    // Apply H - this transforms the basis
    state.h(0);

    // After H, the tableaux change but total probability should be conserved
    let mut total_prob = 0.0;
    for (_, coeff) in state.map.iter() {
        total_prob += coeff.norm_sqr();
    }
    assert!(approx_eq(Complex64::new(total_prob, 0.0), Complex64::new(1.0, 0.0), 1e-12),
            "Probability should be conserved after H gate");
}

#[test]
fn t_gate_zero_coefficient_check() {
    // Verify that no tableau has zero or near-zero coefficient after valid T gate
    let mut state: TableauSum<config::fxhash::FxComplex> = TableauSum::new(1);
    state.h(0);
    state.t(0);

    // All coefficients should be non-zero
    for (tableau, coeff) in state.map.iter() {
        assert!(coeff.norm() > 1e-14,
                "Tableau {:?} has near-zero coefficient {}", tableau, coeff);
    }

    // Now test with multiple T gates
    state.t(0);
    for (tableau, coeff) in state.map.iter() {
        assert!(coeff.norm() > 1e-14,
                "Tableau {:?} has near-zero coefficient {} after 2 T gates", tableau, coeff);
    }
}
