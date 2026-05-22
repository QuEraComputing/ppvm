//! End-to-end fixture coverage: parse + resolve + run each `.sst` file in
//! this directory via the public `PPVM` API.

use ppvm_vihaco::composite::PPVM;

#[test]
fn bell_sst_runs_and_records_two_measurements() {
    let mut machine = PPVM::default();
    machine
        .run_file("tests/bell.sst")
        .unwrap_or_else(|e| panic!("run bell.sst: {e:?}"));
    assert_eq!(machine.measurement_record().len(), 2);
}

#[test]
fn hello_circuit_sst_parses_and_runs() {
    let mut machine = PPVM::default();
    machine
        .run_file("tests/hello_circuit.sst")
        .unwrap_or_else(|e| panic!("run hello_circuit.sst: {e:?}"));
    // hello_circuit.sst applies H + CNOT + RX(0.1); no measurements.
    assert_eq!(machine.measurement_record().len(), 0);
}

#[test]
fn run_file_via_library_helper() {
    let machine =
        ppvm_vihaco::run_file("tests/bell.sst").unwrap_or_else(|e| panic!("run bell.sst: {e:?}"));
    assert_eq!(machine.measurement_record().len(), 2);
}

#[test]
fn function_call_jumps_into_callee_body() {
    // `function_call.sst` has main `call` into `@run_circuit`, which puts q1
    // in |+>, measures it, and `halt`s. Verifies CallPatch resolves the
    // symbolic target and op_call actually transfers control there.
    let machine = ppvm_vihaco::run_file("tests/function_call.sst")
        .unwrap_or_else(|e| panic!("run function_call.sst: {e:?}"));
    let record = machine.measurement_record();
    assert_eq!(record.len(), 1, "expected exactly one measurement");
    assert_eq!(record[0].len(), 1);
    assert!(record[0][0].is_some(), "measurement should not be lost");
}

#[test]
#[ignore]
fn function_call_returns() {
    let machine = ppvm_vihaco::run_file("tests/function_call_ret.sst")
        .unwrap_or_else(|e| panic!("run function_call.sst: {e:?}"));
    let record = machine.measurement_record();
    assert_eq!(record.len(), 1, "expected exactly one measurement");
    assert_eq!(record[0].len(), 1);
    assert!(record[0][0].is_some(), "measurement should not be lost");
}

#[test]
fn branch_on_outcome_deterministic_x_path() {
    // `branch_on_outcome_x.sst` applies X to q0 instead of H, so the outcome
    // is deterministically 1. The cond_br must therefore take the @one path,
    // which flips q1 before measuring it, yielding m1 = 1 as well.
    let machine = ppvm_vihaco::run_file("tests/branch_on_outcome_x.sst")
        .unwrap_or_else(|e| panic!("run branch_on_outcome_x.sst: {e:?}"));
    let record = machine.measurement_record();
    assert_eq!(record.len(), 2, "expected exactly two measurements");
    assert_eq!(record[0], vec![Some(true)], "X-prepared q0 must measure 1");
    assert_eq!(record[1], vec![Some(true)], "branch must have flipped q1");
}

#[test]
fn branch_on_outcome_statistics_balanced_and_invariant_holds() {
    // `branch_on_outcome.sst` puts q0 in |+>, so its measurement is a fair
    // coin. The branch then flips q1 iff the outcome was 1, making m1 == m0
    // an invariant on every shot. Run many shots and check both properties.
    const SHOTS: usize = 400;
    let mut ones = 0usize;
    for _ in 0..SHOTS {
        let machine = ppvm_vihaco::run_file("tests/branch_on_outcome.sst")
            .unwrap_or_else(|e| panic!("run branch_on_outcome.sst: {e:?}"));
        let record = machine.measurement_record();
        assert_eq!(record.len(), 2);
        let m0 = record[0][0].expect("q0 should not be lost");
        let m1 = record[1][0].expect("q1 should not be lost");
        assert_eq!(m0, m1, "branch must steer q1 to match q0 on every shot");
        if m0 {
            ones += 1;
        }
    }
    // Fair coin with SHOTS=400: mean=200, stddev=10. ±6σ window catches a
    // truly broken RNG without flaking on a healthy one.
    let lo = SHOTS / 2 - 60;
    let hi = SHOTS / 2 + 60;
    assert!(
        (lo..=hi).contains(&ones),
        "expected {lo}..={hi} ones out of {SHOTS}, got {ones}"
    );
}
