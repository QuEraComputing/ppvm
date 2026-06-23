// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use stim_parser_2::prelude::*;

#[test]
fn unknown_instruction_returns_error() {
    let err = parse("FROBNICATE 0").unwrap_err();
    let d = err.iter().next().unwrap();
    assert_eq!(d.code, Some("unknown-instruction"));
    assert!(d.message.contains("FROBNICATE"), "message: {}", d.message);
    assert_eq!(
        err.to_string(),
        "error at line 1, col 1: unknown instruction 'FROBNICATE'"
    );
}

#[test]
fn arg_count_mismatch() {
    let err = parse("DEPOLARIZE1(0.1, 0.2) 0").unwrap_err();
    let d = err.iter().next().unwrap();
    assert_eq!(d.code, Some("arg-count"));
    assert!(d.message.contains("DEPOLARIZE1"), "message: {}", d.message);
    // expected 1, found 2
    assert!(
        d.message.contains('1') && d.message.contains('2'),
        "message: {}",
        d.message
    );
    assert!(
        err.to_string().starts_with("error at line 1,"),
        "display: {err}"
    );
}

#[test]
fn target_pair_violation() {
    let err = parse("CX 0 1 2").unwrap_err();
    let d = err.iter().next().unwrap();
    assert_eq!(d.code, Some("target-count"));
    assert!(d.message.contains("CX"), "message: {}", d.message);
    // divisor 2, found 3
    assert!(
        d.message.contains('2') && d.message.contains('3'),
        "message: {}",
        d.message
    );
}

#[test]
fn at_least_one_target_required_for_h() {
    let err = parse("H").unwrap_err();
    let code = err.iter().next().unwrap().code;
    assert!(
        code == Some("target-count") || code == Some("syntax"),
        "code was {code:?}"
    );
}

#[test]
fn invalid_target_yields_syntax_error() {
    // Reference grouped this under `ParseError::Syntax`; the new pipeline
    // tokenizes `abc` as a target then rejects it in the validator, so the
    // code is the more specific `invalid-target`. Same accept/reject decision.
    let err = parse("H abc").unwrap_err();
    let code = err.iter().next().unwrap().code;
    assert!(
        code == Some("syntax") || code == Some("invalid-target"),
        "code was {code:?}"
    );
}

#[test]
fn unclosed_bracket_yields_syntax_error() {
    let err = parse("S[T 0").unwrap_err();
    assert_eq!(err.iter().next().unwrap().code, Some("syntax"));
}

#[test]
fn line_numbers_in_errors_are_correct() {
    let err = parse("X 0\nY 0\nFROBNICATE 0").unwrap_err();
    let s = err.to_string();
    assert!(s.starts_with("error at line 3,"), "message: {s}");
    assert_eq!(err.iter().next().unwrap().code, Some("unknown-instruction"));
}

#[test]
fn syntax_error_includes_line_and_col() {
    // Pin the new Display behavior: line and col both appear in the message.
    let err = parse("H 0\nH abc").unwrap_err();
    let s = err.to_string();
    assert!(s.contains("line 2"), "message was: {s}");
    assert!(s.contains("col"), "message was: {s}");
}

#[test]
fn mpad_zero_targets_is_target_count_error() {
    let err = parse("MPAD").expect_err("must reject");
    let d = err.iter().next().unwrap();
    assert_eq!(d.code, Some("target-count"));
    assert!(d.message.contains("MPAD"), "message: {}", d.message);
}

#[test]
fn mpad_two_args_is_arg_count_error() {
    let err = parse("MPAD(0.1, 0.2) 0").expect_err("must reject");
    let d = err.iter().next().unwrap();
    assert_eq!(d.code, Some("arg-count"));
    assert!(d.message.contains("MPAD"), "message: {}", d.message);
}

#[test]
fn mpad_with_rec_target_is_syntax_error() {
    let err = parse("MPAD rec[-1]").expect_err("must reject");
    let code = err.iter().next().unwrap().code;
    assert!(
        code == Some("syntax") || code == Some("invalid-target"),
        "code was {code:?}"
    );
}
