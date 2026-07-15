// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Printer-fixpoint property checked over generated inputs.
//!
//! The hand-written corpus in `tests/roundtrip.rs` covers ~22 cases.
//! This file lifts the same property to a property-based test: anything
//! the parser accepts must round-trip through `print → parse → print` to
//! a byte-identical fixpoint.
//!
//! See `tests/roundtrip.rs` for the prose explanation of the property.

use proptest::prelude::*;
use stim_parser::prelude::parse;
use stim_parser::prelude::parse_extended;

/// Strategy: build a string out of plausible Stim fragments. Most of
/// these are syntactically valid; some are slightly off so we exercise
/// both the accept and reject branches. The fixpoint property only
/// asserts on the accept side.
fn instruction_fragment() -> impl Strategy<Value = String> {
    prop_oneof![
        // Bare Clifford gates
        Just("H 0\n".to_string()),
        Just("X 0\n".to_string()),
        Just("Y 0\n".to_string()),
        Just("Z 0\n".to_string()),
        Just("S 0\n".to_string()),
        Just("S_DAG 1\n".to_string()),
        Just("I 0\n".to_string()),
        Just("CX 0 1\n".to_string()),
        Just("CZ 1 2\n".to_string()),
        Just("CNOT 0 3\n".to_string()),
        // Reset / measure
        Just("R 0 1\n".to_string()),
        Just("M 0\n".to_string()),
        Just("MZ 0 1 2\n".to_string()),
        Just("MR 0\n".to_string()),
        Just("M(0.001) 0\n".to_string()),
        // Tagged sugar
        Just("S[T] 0\n".to_string()),
        Just("S_DAG[T] 1\n".to_string()),
        Just("I[R_X(theta=0.5*pi)] 0\n".to_string()),
        Just("I[R_Y(theta=1.25*pi)] 1\n".to_string()),
        Just("I[R_Z(theta=-0.5*pi)] 2\n".to_string()),
        Just("I[U3(theta=0.5*pi, phi=1.0*pi, lambda=1.5*pi)] 0\n".to_string()),
        // Noise
        Just("DEPOLARIZE1(0.05) 0\n".to_string()),
        Just("DEPOLARIZE2(0.05) 0 1\n".to_string()),
        Just("PAULI_CHANNEL_1(0.01, 0.02, 0.03) 0\n".to_string()),
        Just("X_ERROR(0.1) 0\n".to_string()),
        Just("I_ERROR[loss](0.01) 0\n".to_string()),
        Just("I_ERROR[correlated_loss](0.1, 0.05, 0.05) 0 1\n".to_string()),
        // MPAD
        Just("MPAD 0 1 0\n".to_string()),
        Just("MPAD(0.1) 1 1 0\n".to_string()),
        // Annotations (with and without rec[-1])
        Just("TICK\n".to_string()),
        Just("DETECTOR\n".to_string()),
        Just("DETECTOR rec[-1]\n".to_string()),
        Just("OBSERVABLE_INCLUDE(0) rec[-1]\n".to_string()),
        Just("QUBIT_COORDS(0, 0) 0\n".to_string()),
        // REPEAT
        Just("REPEAT 3 {\n    H 0\n    M 0\n}\n".to_string()),
        Just("REPEAT 2 {\n    REPEAT 3 {\n        X 0\n    }\n}\n".to_string()),
        // Stylistic noise the printer normalizes away
        Just("# leading\n".to_string()),
        Just("\n".to_string()),
        Just("H 0  # trail\n".to_string()),
    ]
}

fn program_source() -> impl Strategy<Value = String> {
    prop::collection::vec(instruction_fragment(), 0..16).prop_map(|frags| frags.concat())
}

/// Assert the printer-fixpoint property at the raw-AST level.
fn check_raw_fixpoint(src: &str) {
    let Ok(ast1) = parse(src) else { return };
    let s1 = format!("{ast1}");
    let ast2 = parse(&s1)
        .unwrap_or_else(|e| panic!("reparse of printed output failed: {e}\n--printed--\n{s1}"));
    let s2 = format!("{ast2}");
    assert_eq!(s1, s2, "raw printer is not a fixpoint");
}

/// Assert the printer-fixpoint property at the extended-AST level.
fn check_extended_fixpoint(src: &str) {
    let Ok(ast1) = parse_extended(src) else {
        return;
    };
    let s1 = format!("{ast1}");
    let ast2 = parse_extended(&s1)
        .unwrap_or_else(|e| panic!("reparse of printed output failed: {e}\n--printed--\n{s1}"));
    let s2 = format!("{ast2}");
    assert_eq!(s1, s2, "extended printer is not a fixpoint");
}

/// A decimal with 0–4 fractional digits in `[-4, 4]`, rendered as a string.
/// Mirrors the angles a user actually writes; many are not binary-friendly,
/// so `theta/PI` would print a rounding tail without the printer's recovery.
fn decimal_coeff() -> impl Strategy<Value = String> {
    (-40_000i32..=40_000, 0u32..=4).prop_map(|(n, scale)| {
        let v = f64::from(n) / 10f64.powi(scale as i32);
        format!("{v}")
    })
}

proptest! {
    #[test]
    fn raw_printer_is_fixpoint_on_fragments(src in program_source()) {
        check_raw_fixpoint(&src);
    }

    #[test]
    fn extended_printer_is_fixpoint_on_fragments(src in program_source()) {
        check_extended_fixpoint(&src);
    }

    /// Rotation/U3 angles are stored in radians as `c*PI` and re-emitted in
    /// the `<c>*pi` form. For any decimal coefficient a user might write, the
    /// printer must (a) round-trip losslessly — `parse → print → parse` must
    /// recover the exact stored radians — and (b) stay a byte-for-byte
    /// fixpoint, never degrading into a `0.7599999999999999*pi` tail.
    #[test]
    fn rotation_pi_coeff_round_trips(c in decimal_coeff()) {
        let src = format!("I[R_Z(theta={c}*pi)] 0\n");
        let ast1 = parse_extended(&src).expect("parse");
        let theta1 = rotation_theta(&ast1);

        let s1 = format!("{ast1}");
        prop_assert!(
            !s1.contains("999999") && !s1.contains("000000"),
            "printed a rounding tail for theta={c}*pi: {s1}"
        );

        let ast2 = parse_extended(&s1).expect("reparse");
        prop_assert_eq!(theta1.to_bits(), rotation_theta(&ast2).to_bits(),
            "theta not recovered exactly for theta={}*pi", c);
        prop_assert_eq!(format!("{ast2}"), s1, "printer not a fixpoint");
    }
}

/// Extract the single rotation's `theta` from a one-instruction program.
fn rotation_theta(ast: &stim_parser::prelude::ExtendedProgram) -> f64 {
    use stim_parser::prelude::ExtendedInstruction::Rotation;
    match &ast.instructions[0] {
        Rotation { theta, .. } => *theta,
        other => panic!("expected a rotation, got {other:?}"),
    }
}
