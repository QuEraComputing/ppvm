// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

use stim_parser::prelude::*;

fn approx_eq(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} vs {b}");
}

#[test]
fn parse_depolarize1_one_arg() {
    let p = parse("DEPOLARIZE1(0.5) 0").unwrap();
    match &p.instructions[0] {
        RawInstruction::Noise {
            name,
            args,
            targets,
            ..
        } => {
            assert_eq!(*name, NoiseName::Depolarize1);
            assert_eq!(args.len(), 1);
            approx_eq(args[0], 0.5);
            assert_eq!(targets, &[0]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_pauli_channel_1_three_args() {
    let p = parse("PAULI_CHANNEL_1(0.1, 0.2, 0.3) 5").unwrap();
    match &p.instructions[0] {
        RawInstruction::Noise {
            name,
            args,
            targets,
            ..
        } => {
            assert_eq!(*name, NoiseName::PauliChannel1);
            assert_eq!(args.len(), 3);
            approx_eq(args[0], 0.1);
            approx_eq(args[1], 0.2);
            approx_eq(args[2], 0.3);
            assert_eq!(targets, &[5]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_depolarize2_pairs() {
    let p = parse("DEPOLARIZE2(0.5) 0 1 2 3").unwrap();
    match &p.instructions[0] {
        RawInstruction::Noise { args, targets, .. } => {
            approx_eq(args[0], 0.5);
            assert_eq!(targets, &[0, 1, 2, 3]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_pauli_channel_2_fifteen_args() {
    let src = "PAULI_CHANNEL_2(0.1, 0.12, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001, 0.0005, 0.0003, 0.0002) 4 3";
    let p = parse(src).unwrap();
    match &p.instructions[0] {
        RawInstruction::Noise { args, targets, .. } => {
            assert_eq!(args.len(), 15);
            approx_eq(args[0], 0.1);
            approx_eq(args[14], 0.0002);
            assert_eq!(targets, &[4, 3]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_args_accept_pi_expr() {
    // Pi expressions are usually inside tags, but the spec says the args parser
    // accepts them too, matching today's stim.rs (which uses parse_pi_expr).
    let p = parse("X_ERROR(0.5*pi) 0").unwrap();
    match &p.instructions[0] {
        RawInstruction::Noise { args, .. } => {
            approx_eq(args[0], 0.5 * std::f64::consts::PI);
        }
        other => panic!("{other:?}"),
    }
}
