// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Parsing of measurement-record gate controls (`rec[-k]`) and the X/Y-basis
//! reset gates (`RX`/`RY`).

use stim_parser::prelude::*;

fn first_gate(src: &str) -> (GateName, Vec<Target>) {
    let prog = parse(src).expect("must parse");
    match &prog.instructions[0] {
        RawInstruction::Gate { name, targets, .. } => (*name, targets.clone()),
        other => panic!("expected Gate, got {other:?}"),
    }
}

#[test]
fn cx_with_record_control_parses() {
    let (name, targets) = first_gate("CX rec[-1] 1");
    assert_eq!(name, GateName::CX);
    assert_eq!(targets, vec![Target::Rec(1), Target::Qubit(1)]);
}

#[test]
fn deeper_record_lookback_parses() {
    let (_, targets) = first_gate("CZ rec[-5] 2");
    assert_eq!(targets, vec![Target::Rec(5), Target::Qubit(2)]);
}

#[test]
fn rx_ry_parse_as_basis_resets() {
    assert_eq!(first_gate("RX 0").0, GateName::ResetX);
    assert_eq!(first_gate("RY 3").0, GateName::ResetY);
}

#[test]
fn zero_and_positive_record_lookbacks_are_rejected() {
    // Stim record references must be strictly negative.
    assert!(parse("CX rec[0] 1").is_err());
    assert!(parse("CX rec[1] 1").is_err());
    assert!(parse("CX rec[] 1").is_err());
}

#[test]
fn record_target_on_measurement_is_rejected() {
    // `rec[-k]` is only a gate control; measurements take qubit targets.
    assert!(parse("M rec[-1]").is_err());
}

#[test]
fn record_control_round_trips_through_display() {
    let src = "CX rec[-1] 1";
    let printed = parse(src).unwrap().to_string();
    assert_eq!(printed.trim(), "CX rec[-1] 1");
    // And the printed form re-parses identically.
    assert_eq!(parse(printed.trim()).unwrap(), parse(src).unwrap());
}
