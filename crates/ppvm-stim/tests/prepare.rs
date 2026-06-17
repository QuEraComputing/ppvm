// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_stim::{ExecError, parse_extended, prepare};

fn err_from_src(src: &str) -> ExecError {
    let prog = parse_extended(src).expect("parse_extended");
    prepare(&prog).expect_err("must reject")
}

#[test]
fn unsupported_swap_rejected() {
    let ExecError::Unsupported { name, line } = err_from_src("SWAP 0 1") else {
        panic!("expected ExecError::Unsupported");
    };
    assert_eq!(name, "SWAP");
    assert_eq!(line, 1);
}

#[test]
fn unsupported_mzz_rejected() {
    // Two-qubit Pauli-product measurements remain unsupported.
    let e = err_from_src("MZZ 0 1");
    assert!(matches!(e, ExecError::Unsupported { .. }));
}

#[test]
fn xy_basis_measure_and_reset_accepted() {
    // MX/MY/MRX/MRY and RX/RY are now executed via basis-change decomposition.
    let prog = parse_extended("RX 0\nRY 1\nMX 0\nMY 1\nMRX 0\nMRY 1\n").expect("parse_extended");
    assert_eq!(prepare(&prog), Ok(()));
}

#[test]
fn record_control_accepted_on_controlled_paulis() {
    let prog =
        parse_extended("M 0\nCX rec[-1] 1\nCY rec[-1] 2\nCZ rec[-1] 3\n").expect("parse_extended");
    assert_eq!(prepare(&prog), Ok(()));
}

#[test]
fn record_control_rejected_on_single_qubit_gate() {
    let ExecError::InvalidRecordControl { name, .. } = err_from_src("M 0\nH rec[-1]") else {
        panic!("expected ExecError::InvalidRecordControl");
    };
    assert_eq!(name, "H");
}

#[test]
fn record_in_target_slot_rejected() {
    // The Pauli target may never be a measurement record (Stim: "measurement
    // record editing is not supported").
    let ExecError::InvalidRecordControl { name, .. } = err_from_src("M 0\nCX 1 rec[-1]") else {
        panic!("expected ExecError::InvalidRecordControl");
    };
    assert_eq!(name, "CX");
}

#[test]
fn unsupported_heralded_erase_rejected() {
    let e = err_from_src("HERALDED_ERASE(0.1) 0");
    assert!(matches!(e, ExecError::Unsupported { .. }));
}

#[test]
fn unsupported_swap_inside_repeat_rejected() {
    let ExecError::Unsupported { name, line } = err_from_src("REPEAT 3 {\n    SWAP 0 1\n}\n")
    else {
        panic!("expected ExecError::Unsupported");
    };
    assert_eq!(name, "SWAP");
    assert_eq!(line, 2);
}

#[test]
fn supported_structural_instructions_are_not_rejected_by_prepare() {
    let prog =
        parse_extended("MPAD 0 1\nI_ERROR[correlated_loss](0.5) 0 1\n").expect("parse_extended");
    assert_eq!(prepare(&prog), Ok(()));
}

#[test]
fn invalid_measure_probability_rejected() {
    let ExecError::InvalidProbability { name, line, value } = err_from_src("M(2.0) 0") else {
        panic!("expected ExecError::InvalidProbability");
    };
    assert_eq!(name, "M");
    assert_eq!(line, 1);
    assert_eq!(value, 2.0);
}

#[test]
fn invalid_mpad_probability_rejected() {
    let ExecError::InvalidProbability { name, line, value } = err_from_src("MPAD(-0.1) 0") else {
        panic!("expected ExecError::InvalidProbability");
    };
    assert_eq!(name, "MPAD");
    assert_eq!(line, 1);
    assert_eq!(value, -0.1);
}
