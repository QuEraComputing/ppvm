//! End-to-end fixture coverage: parse + resolve + run each `.sst` file in
//! this directory via the public `PPVM` API.

use ppvm_vihaco::composite::PPVM;
use ppvm_vihaco::measurements::MeasurementOutcome;

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
    assert!(
        record[0][0] != MeasurementOutcome::Lost,
        "measurement should not be lost"
    );
}

#[test]
fn function_call_returns() {
    let machine = ppvm_vihaco::run_file("tests/function_call_ret.sst")
        .unwrap_or_else(|e| panic!("run function_call.sst: {e:?}"));
    let record = machine.measurement_record();
    assert_eq!(record.len(), 1, "expected exactly one measurement");
    assert_eq!(record[0].len(), 1);
    assert!(
        record[0][0] != MeasurementOutcome::Lost,
        "measurement should not be lost"
    );
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
    assert_eq!(
        record[0].as_slice(),
        &[MeasurementOutcome::One],
        "X-prepared q0 must measure 1"
    );
    assert_eq!(
        record[1].as_slice(),
        &[MeasurementOutcome::One],
        "branch must have flipped q1"
    );
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
        let m0 = record[0][0];
        let m1 = record[1][0];
        assert_eq!(m0, m1, "branch must steer q1 to match q0 on every shot");
        assert!(
            m0 != MeasurementOutcome::Lost,
            "measurement should not be lost"
        );
        if m0 == MeasurementOutcome::One {
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

#[test]
fn function_call_branch_on_both_returned_values() {
    // `function_call_branch_both.sst`: helper returns the tri-state outcome
    // (0/1/Lost) via `ret 1`. Main first branches on is_lost, then on the
    // 0/1 outcome, steering q1 to |1> on the lost path and on the
    // kept-outcome=1 path, leaving q1 in |0> only on the kept-outcome=0
    // path. With loss prob 0.5 and a |+> prep:
    //   P(m1 = 1) = P(lost) + P(kept) · P(outcome = 1 | kept)
    //             = 0.5     + 0.5 · 0.5  = 0.75
    //   P(m0 = lost) = 0.5
    const SHOTS: usize = 400;
    let mut q0_lost = 0usize;
    let mut q1_ones = 0usize;
    for _ in 0..SHOTS {
        let machine = ppvm_vihaco::run_file("tests/function_call_branch_both.sst")
            .unwrap_or_else(|e| panic!("run function_call_branch_both.sst: {e:?}"));
        let record = machine.measurement_record();
        assert_eq!(record.len(), 2, "expected exactly two measurements");
        assert_eq!(record[0].len(), 1);
        assert_eq!(record[1].len(), 1);
        if record[0][0] == MeasurementOutcome::Lost {
            q0_lost += 1;
        }
        if record[1][0] == MeasurementOutcome::One {
            q1_ones += 1;
        }
    }
    // P(lost) = 0.5, SHOTS=400: mean=200, stddev=10. ±6σ window.
    assert!(
        (140..=260).contains(&q0_lost),
        "expected ~200 lost shots, got {q0_lost}"
    );
    // P(m1=1) = 0.75, SHOTS=400: mean=300, stddev≈8.66. ±6σ → ~248..352.
    assert!(
        (240..=360).contains(&q1_ones),
        "expected ~300 q1=true shots, got {q1_ones}"
    );
}
