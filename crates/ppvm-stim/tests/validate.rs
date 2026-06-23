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
fn unsupported_mx_rejected() {
    let e = err_from_src("MX 0");
    assert!(matches!(e, ExecError::Unsupported { .. }));
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
