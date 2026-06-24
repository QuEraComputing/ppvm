// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_stim::{ExecError, parse_extended, validate};

fn err_from_src(src: &str) -> ExecError {
    let prog = parse_extended(src).expect("parse_extended");
    validate(&prog).expect_err("must reject")
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
    assert_eq!(validate(&prog), Ok(()));
}

#[test]
fn record_control_accepted_on_controlled_paulis() {
    let prog =
        parse_extended("M 0\nCX rec[-1] 1\nCY rec[-1] 2\nCZ rec[-1] 3\n").expect("parse_extended");
    assert_eq!(validate(&prog), Ok(()));
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
fn record_control_with_no_prior_measurement_rejected() {
    // Stim raises IndexError for `CX rec[-1] 0` with an empty measurement record;
    // we reject it rather than silently applying no correction.
    let ExecError::InvalidRecordControl { name, line, .. } = err_from_src("CX rec[-1] 0") else {
        panic!("expected ExecError::InvalidRecordControl");
    };
    assert_eq!(name, "CX");
    assert_eq!(line, 1);
}

#[test]
fn record_control_looking_back_past_record_start_rejected() {
    // One measurement recorded; `rec[-2]` looks back before the record start.
    let ExecError::InvalidRecordControl { name, .. } = err_from_src("M 0\nCX rec[-2] 1") else {
        panic!("expected ExecError::InvalidRecordControl");
    };
    assert_eq!(name, "CX");
}

#[test]
fn in_range_record_control_accepted() {
    // The boundary case `rec[-1]` with exactly one prior measurement is valid.
    let prog = parse_extended("M 0\nCX rec[-1] 1\n").expect("parse_extended");
    assert_eq!(validate(&prog), Ok(()));
}

#[test]
fn record_control_out_of_range_on_first_repeat_iteration_rejected() {
    // Inside a REPEAT the first iteration sees the shortest record. `rec[-2]`
    // needs two prior measurements but the first iteration has only produced
    // one when the CX runs, so Stim would error on iteration 1 — reject it.
    let ExecError::InvalidRecordControl { .. } =
        err_from_src("REPEAT 3 {\n    M 0\n    CX rec[-2] 1\n}\n")
    else {
        panic!("expected ExecError::InvalidRecordControl");
    };
}

#[test]
fn record_control_referencing_prior_iteration_accepted() {
    // `rec[-2]` after two measurements per iteration is in range from the first
    // iteration on, so the loop validates.
    let prog = parse_extended("REPEAT 3 {\n    M 0\n    M 1\n    CX rec[-2] 2\n}\n")
        .expect("parse_extended");
    assert_eq!(validate(&prog), Ok(()));
}

#[test]
fn mpp_repeated_qubit_rejected() {
    // `MPP X0*X0` repeats qubit 0; Stim folds it to identity, but our CX-ladder
    // gadget assumes distinct qubits, so we reject it with a clear error.
    let ExecError::InvalidPauliProduct { line, message } = err_from_src("MPP X0*X0") else {
        panic!("expected ExecError::InvalidPauliProduct");
    };
    assert_eq!(line, 1);
    assert!(message.contains("distinct qubits"), "message: {message}");
    assert!(
        message.contains("issue"),
        "message should point to the issue tracker: {message}"
    );
}

#[test]
fn mpp_anti_hermitian_product_rejected() {
    // `MPP Z0*X0` = iY0 is anti-Hermitian (Stim raises ValueError); it also
    // repeats qubit 0, so the distinct-qubit check rejects it.
    assert!(matches!(
        err_from_src("MPP Z0*X0"),
        ExecError::InvalidPauliProduct { .. }
    ));
}

#[test]
fn mpp_repeated_qubit_in_later_product_rejected() {
    // The check runs per product; a clean first product does not mask a repeat
    // in a later one.
    assert!(matches!(
        err_from_src("MPP Z0*Z1 X2*X2"),
        ExecError::InvalidPauliProduct { .. }
    ));
}

#[test]
fn mpp_distinct_qubit_products_accepted() {
    let prog = parse_extended("MPP X0*Y1*Z2 Z3*Z4\n").expect("parse_extended");
    assert_eq!(validate(&prog), Ok(()));
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
fn supported_structural_instructions_are_not_rejected_by_validate() {
    let prog =
        parse_extended("MPAD 0 1\nI_ERROR[correlated_loss](0.5) 0 1\n").expect("parse_extended");
    assert_eq!(validate(&prog), Ok(()));
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
