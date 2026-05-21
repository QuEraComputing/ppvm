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
    let machine = ppvm_vihaco::run_file("tests/bell.sst")
        .unwrap_or_else(|e| panic!("run bell.sst: {e:?}"));
    assert_eq!(machine.measurement_record().len(), 2);
}
