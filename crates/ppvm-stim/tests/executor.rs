use ppvm_runtime::config::indexmap::ByteFxHashF64;
use ppvm_stim::{execute, parse_extended};
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<ByteFxHashF64<1>, usize>;

fn run(src: &str, n_qubits: usize) -> (Vec<Option<bool>>, Tab) {
    let prog = parse_extended(src).expect("parse_extended");
    let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
    let results = execute(&prog, &mut tab).expect("execute");
    (results, tab)
}

#[test]
fn x_then_measure_returns_one() {
    let (results, _) = run("X 0\nM 0", 1);
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn fresh_state_measures_zero() {
    let (results, _) = run("M 0", 1);
    assert_eq!(results, vec![Some(false)]);
}

#[test]
fn h_h_returns_zero() {
    let (results, _) = run("H 0\nH 0\nM 0", 1);
    assert_eq!(results, vec![Some(false)]);
}

#[test]
fn rx_pi_flips_qubit() {
    let (results, _) = run("I[R_X(theta=1.0*pi)] 0\nM 0", 1);
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn u3_pi_flip_via_y_axis() {
    let (results, _) = run("I[U3(theta=1.0*pi, phi=0.0, lambda=0.0)] 0\nM 0", 1);
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn t_gate_via_s_t_tag_no_op_on_zero() {
    let (results, _) = run("S[T] 0\nM 0", 1);
    assert_eq!(results, vec![Some(false)]);
}

#[test]
fn cnot_creates_bell_correlation() {
    for _ in 0..32 {
        let (results, _) = run("H 0\nCX 0 1\nM 0 1", 2);
        assert_eq!(results.len(), 2);
        assert_eq!(results[0], results[1]);
    }
}

#[test]
fn cnot_alias_equivalents() {
    let (cx, _) = run("X 0\nCX 0 1\nM 1", 2);
    let (cnot, _) = run("X 0\nCNOT 0 1\nM 1", 2);
    let (zcx, _) = run("X 0\nZCX 0 1\nM 1", 2);
    assert_eq!(cx, vec![Some(true)]);
    assert_eq!(cnot, vec![Some(true)]);
    assert_eq!(zcx, vec![Some(true)]);
}

#[test]
fn mr_resets_qubit_after_measure() {
    let (results, _) = run("X 0\nMR 0\nM 0", 1);
    assert_eq!(results, vec![Some(true), Some(false)]);
}

#[test]
fn loss_channel_with_p1_marks_qubit_lost() {
    let prog = parse_extended("I_ERROR[loss](1.0) 0").unwrap();
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    execute(&prog, &mut tab).unwrap();
    assert!(tab.is_lost[0]);
}

#[test]
fn repeat_executes_body_n_times() {
    let (results, _) = run("REPEAT 2 { X 0 }\nM 0", 1);
    assert_eq!(results, vec![Some(false)]);
    let (results, _) = run("REPEAT 3 { X 0 }\nM 0", 1);
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn comments_and_annotations_are_no_ops() {
    let (results, _) = run("# c\nQUBIT_COORDS(0,0) 0\nX 0\nTICK\nM 0", 1);
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn measurement_buffer_is_pre_sized() {
    let prog = parse_extended("X 0\nM 0 1 2 3 4").unwrap();
    assert_eq!(prog.measurement_count(), 5);
    let mut tab: Tab = GeneralizedTableau::new(5, 1e-10);
    let results = execute(&prog, &mut tab).unwrap();
    assert_eq!(results.len(), 5);
}

#[test]
fn sample_runs_n_shots_each_with_fresh_tableau() {
    use ppvm_stim::sample;
    let prog = parse_extended("X 0\nM 0").unwrap();
    let shots = sample::<_, _, _, _>(&prog, 5, || {
        GeneralizedTableau::<ByteFxHashF64<1>, usize>::new(1, 1e-10)
    })
    .unwrap();
    assert_eq!(shots.len(), 5);
    for shot in &shots {
        assert_eq!(shot, &vec![Some(true)]);
    }
}

#[test]
fn sample_zero_shots_returns_empty() {
    use ppvm_stim::sample;
    let prog = parse_extended("X 0\nM 0").unwrap();
    let shots = sample::<_, _, _, _>(&prog, 0, || {
        GeneralizedTableau::<ByteFxHashF64<1>, usize>::new(1, 1e-10)
    })
    .unwrap();
    assert!(shots.is_empty());
}

#[test]
fn sample_random_h_distribution_within_3_sigma() {
    // H 0; M 0 — over 4096 shots, expect ≈50% ones, allow 3σ slack.
    use ppvm_stim::sample;
    let prog = parse_extended("H 0\nM 0").unwrap();
    let n = 4096;
    let mut seed_counter: u64 = 0;
    let shots = sample::<_, _, _, _>(&prog, n, || {
        seed_counter += 1;
        GeneralizedTableau::<ByteFxHashF64<1>, usize>::new_with_seed(1, 1e-10, seed_counter)
    })
    .unwrap();
    let ones = shots.iter().filter(|s| s[0] == Some(true)).count();
    let mean = (n / 2) as f64;
    let std = (n as f64 * 0.25).sqrt();
    assert!(
        (ones as f64 - mean).abs() < 3.0 * std,
        "got {ones} ones, mean={mean}, std={std}"
    );
}

// ============================================================
// Migrated from ppvm-tableau/tests/gates.rs section 8
// ============================================================

fn run_str(src: &str, tab: &mut Tab) -> Vec<Option<bool>> {
    ppvm_stim::run_string(src, tab).unwrap()
}

#[test]
fn test_stim_mr_measure_and_reset() {
    // MR should measure then reset to |0⟩
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    let results = run_str("X 0\nMR 0", &mut tab);
    // Measurement should give 1 (since qubit was |1⟩)
    assert_eq!(results, vec![Some(true)]);
    // After MR, qubit should be reset to |0⟩
    let results2 = run_str("M 0", &mut tab);
    assert_eq!(results2, vec![Some(false)]);
}

#[test]
fn test_stim_mr_zero_state() {
    // MR on |0⟩ should give 0 and leave in |0⟩
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    let results = run_str("MR 0", &mut tab);
    assert_eq!(results, vec![Some(false)]);
    let results2 = run_str("M 0", &mut tab);
    assert_eq!(results2, vec![Some(false)]);
}

#[test]
fn test_stim_cy_gate() {
    // CY should entangle qubits like CX but with Y-basis on target
    // CY|10⟩ = |1⟩ ⊗ Y|0⟩ = |1⟩ ⊗ i|1⟩ (up to phase)
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-10);
    let results = run_str("X 0\nCY 0 1\nM 0 1", &mut tab);
    // Control was |1⟩, so CY flips target
    assert_eq!(results[0], Some(true));
    assert_eq!(results[1], Some(true));
}

#[test]
fn test_stim_cy_control_zero() {
    // CY|00⟩ = |00⟩ (control is 0, no action)
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-10);
    let results = run_str("CY 0 1\nM 0 1", &mut tab);
    assert_eq!(results, vec![Some(false), Some(false)]);
}

#[test]
fn test_stim_cz_gate() {
    // CZ|11⟩ = -|11⟩ (phase flip, but measurement outcome same)
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-10);
    let results = run_str("X 0\nX 1\nCZ 0 1\nM 0 1", &mut tab);
    assert_eq!(results, vec![Some(true), Some(true)]);
}

#[test]
fn test_stim_cz_on_zero() {
    // CZ|00⟩ = |00⟩
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-10);
    let results = run_str("CZ 0 1\nM 0 1", &mut tab);
    assert_eq!(results, vec![Some(false), Some(false)]);
}

#[test]
fn test_stim_s_dag() {
    // S_DAG on |0⟩: Z stabilizer unchanged (Z phase invariant)
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    let results = run_str("S_DAG 0\nM 0", &mut tab);
    assert_eq!(results, vec![Some(false)]);
}

#[test]
fn test_stim_s_dag_t_is_t_adj() {
    // S_DAG[T] should be T†. T†T = I on |+⟩ should leave 1 branch.
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    run_str("H 0\nS[T] 0\nS_DAG[T] 0", &mut tab);
    assert_eq!(tab.coefficients.len(), 1, "T†T should cancel to 1 branch");
}

#[test]
fn test_stim_sqrt_x_dag() {
    // SQRT_X_DAG then SQRT_X should compose to identity
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    let results = run_str("SQRT_X_DAG 0\nSQRT_X 0\nM 0", &mut tab);
    assert_eq!(results, vec![Some(false)], "SQRT_X_DAG · SQRT_X = I");
}

#[test]
fn test_stim_sqrt_y_dag() {
    // SQRT_Y_DAG then SQRT_Y should compose to identity
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    let results = run_str("SQRT_Y_DAG 0\nSQRT_Y 0\nM 0", &mut tab);
    assert_eq!(results, vec![Some(false)], "SQRT_Y_DAG · SQRT_Y = I");
}

#[test]
fn test_stim_sqrt_z_is_s() {
    // SQRT_Z should be the same as S
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    let results = run_str("SQRT_Z 0\nM 0", &mut tab);
    assert_eq!(results, vec![Some(false)]);
}

#[test]
fn test_stim_sqrt_z_dag_is_s_adj() {
    // SQRT_Z_DAG then SQRT_Z = I
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    let results = run_str("SQRT_Z_DAG 0\nSQRT_Z 0\nM 0", &mut tab);
    assert_eq!(results, vec![Some(false)]);
}

#[test]
fn test_stim_correlated_loss_simple() {
    // I_ERROR[correlated_loss](1.0) should lose both qubits
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-10);
    run_str("I_ERROR[correlated_loss](1.0) 0 1", &mut tab);
    assert!(
        tab.is_lost[0] && tab.is_lost[1],
        "Both qubits should be lost"
    );
}

#[test]
fn test_stim_correlated_loss_zero_prob() {
    // I_ERROR[correlated_loss](0.0) should not lose any qubits
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-10);
    run_str("I_ERROR[correlated_loss](0.0) 0 1", &mut tab);
    assert!(
        !tab.is_lost[0] && !tab.is_lost[1],
        "No qubits should be lost"
    );
}

#[test]
fn test_stim_comments_and_empty_lines() {
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    let results = run_str(
        "# This is a comment\n\nX 0\n# Another comment\nM 0",
        &mut tab,
    );
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn test_stim_noop_instructions() {
    // TICK, DETECTOR, QUBIT_COORDS, SHIFT_COORDS, OBSERVABLE_INCLUDE should not crash
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    let results = run_str(
        "TICK\nDETECTOR\nQUBIT_COORDS\nSHIFT_COORDS\nX 0\nM 0",
        &mut tab,
    );
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn test_stim_zcx_alias() {
    // ZCX should be equivalent to CX/CNOT
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-10);
    let results = run_str("X 0\nZCX 0 1\nM 0 1", &mut tab);
    assert_eq!(results, vec![Some(true), Some(true)]);
}

#[test]
fn test_stim_zcy_alias() {
    // ZCY should be equivalent to CY
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-10);
    let results = run_str("X 0\nZCY 0 1\nM 0 1", &mut tab);
    assert_eq!(results[0], Some(true));
    assert_eq!(results[1], Some(true));
}

#[test]
fn test_stim_zcz_alias() {
    // ZCZ should be equivalent to CZ
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-10);
    let results = run_str("X 0\nX 1\nZCZ 0 1\nM 0 1", &mut tab);
    assert_eq!(results, vec![Some(true), Some(true)]);
}

// ============================================================
// Measurement readout noise tests
// ============================================================

#[test]
fn measure_noise_zero_equals_noiseless() {
    // X 0; MZ(0.0) 0 should always give Some(true).
    let (results, _) = run("X 0\nMZ(0.0) 0", 1);
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn measure_noise_one_always_flips() {
    // X 0; MZ(1.0) 0 — the true outcome is 1, but readout always flips it to 0.
    for _ in 0..16 {
        let (results, _) = run("X 0\nMZ(1.0) 0", 1);
        assert_eq!(results, vec![Some(false)]);
    }
    // Also check on |0>: MZ(1.0) of |0> always records 1.
    for _ in 0..16 {
        let (results, _) = run("MZ(1.0) 0", 1);
        assert_eq!(results, vec![Some(true)]);
    }
}

#[test]
fn measure_noise_does_not_affect_state() {
    // After X 0; MZ(1.0) 0; MZ 0:
    //   First measurement records 0 (true outcome 1, noise flipped).
    //   Second measurement reads the *true* state (still |1>) and records 1.
    let (results, _) = run("X 0\nMZ(1.0) 0\nMZ 0", 1);
    assert_eq!(results, vec![Some(false), Some(true)]);
}

#[test]
fn mr_noise_one_flips_recorded_but_resets_correctly() {
    // X 0; MR(1.0) 0; MZ 0:
    //   MR(1.0): measures (true outcome 1), records flipped 0, resets to |0>.
    //   MZ 0: measures |0>, records 0.
    let (results, _) = run("X 0\nMR(1.0) 0\nMZ 0", 1);
    assert_eq!(results, vec![Some(false), Some(false)]);
}

#[test]
fn measure_noise_distribution_within_3_sigma() {
    use ppvm_stim::sample;
    // X 0; MZ(0.3) 0 — true outcome is 1, recorded bit flips with prob 0.3.
    // So recorded == 0 with probability 0.3 over many shots.
    let prog = parse_extended("X 0\nMZ(0.3) 0").unwrap();
    let n = 4096usize;
    let mut seed_counter: u64 = 0;
    let shots = sample::<_, _, _, _>(&prog, n, || {
        seed_counter += 1;
        GeneralizedTableau::<ByteFxHashF64<1>, usize>::new_with_seed(1, 1e-10, seed_counter)
    })
    .unwrap();
    let zeros = shots.iter().filter(|s| s[0] == Some(false)).count();
    let mean = (n as f64) * 0.3;
    let std = ((n as f64) * 0.3 * 0.7).sqrt();
    assert!(
        ((zeros as f64) - mean).abs() < 3.0 * std,
        "got {zeros} zeros, expected mean {mean} +/- 3*{std}"
    );
}

#[test]
fn mpad_single_zero_appends_some_false() {
    let (results, _) = run("MPAD 0", 1);
    assert_eq!(results, vec![Some(false)]);
}

#[test]
fn mpad_single_one_appends_some_true() {
    let (results, _) = run("MPAD 1", 1);
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn mpad_multi_bit_in_order() {
    let (results, _) = run("MPAD 0 1 0 1", 1);
    assert_eq!(
        results,
        vec![Some(false), Some(true), Some(false), Some(true)]
    );
}

#[test]
fn mpad_interleaved_with_measurement() {
    let (results, _) = run("X 0\nMPAD 1\nM 0\nMPAD 0", 1);
    assert_eq!(results, vec![Some(true), Some(true), Some(false)]);
}

#[test]
fn mpad_inside_repeat_block_executes_each_iteration() {
    let (results, _) = run("REPEAT 3 {\n    MPAD 1\n}", 1);
    assert_eq!(results, vec![Some(true), Some(true), Some(true)]);
}

#[test]
fn mpad_noise_distribution_within_3_sigma() {
    use ppvm_stim::sample;
    // MPAD(0.3) 0 — pad value is 0; recorded bit flips to 1 with prob 0.3.
    let prog = parse_extended("MPAD(0.3) 0").unwrap();
    let n = 4096usize;
    let mut seed_counter: u64 = 0;
    let shots = sample::<_, _, _, _>(&prog, n, || {
        seed_counter += 1;
        GeneralizedTableau::<ByteFxHashF64<1>, usize>::new_with_seed(1, 1e-10, seed_counter)
    })
    .unwrap();
    let ones = shots.iter().filter(|s| s[0] == Some(true)).count();
    let mean = (n as f64) * 0.3;
    let std = ((n as f64) * 0.3 * 0.7).sqrt();
    assert!(
        ((ones as f64) - mean).abs() < 3.0 * std,
        "got {ones} ones, expected mean {mean} +/- 3*{std}"
    );
}
