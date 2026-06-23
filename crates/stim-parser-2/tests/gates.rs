// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use stim_parser_2::prelude::*;

fn instructions(p: &Program) -> &[Instruction] {
    &p.instructions
}

#[test]
fn parse_single_h_one_target() {
    let p = parse("H 0").expect("must parse");
    assert_eq!(p.instructions.len(), 1);
    match &instructions(&p)[0] {
        Instruction::Gate(GateOp {
            name,
            tags,
            args,
            targets,
            span,
        }) => {
            assert_eq!(*name, GateName::H);
            assert!(tags.is_empty());
            assert!(args.is_empty());
            assert_eq!(targets, &[Target::Qubit(0)]);
            assert_eq!(span.line(&p.line_map), 1);
        }
        other => panic!("expected Gate, got {other:?}"),
    }
}

#[test]
fn parse_h_with_multiple_targets() {
    let p = parse("H 0 1 2 3").expect("must parse");
    match &instructions(&p)[0] {
        Instruction::Gate(GateOp { name, targets, .. }) => {
            assert_eq!(*name, GateName::H);
            assert_eq!(targets, &[0, 1, 2, 3]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_cx_pairs() {
    let p = parse("CX 0 1 2 3").expect("must parse");
    match &instructions(&p)[0] {
        Instruction::Gate(GateOp { name, targets, .. }) => {
            assert_eq!(*name, GateName::CX);
            assert_eq!(targets, &[0, 1, 2, 3]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_each_basic_clifford() {
    for (src, expected) in [
        ("X 0", GateName::X),
        ("Y 0", GateName::Y),
        ("Z 0", GateName::Z),
        ("H 0", GateName::H),
        ("H_XZ 0", GateName::HXZ),
        ("S 0", GateName::S),
        ("S_DAG 0", GateName::SDag),
        ("SQRT_X 0", GateName::SqrtX),
        ("SQRT_Y 0", GateName::SqrtY),
        ("SQRT_Z 0", GateName::SqrtZ),
        ("R 0", GateName::Reset),
        ("RZ 0", GateName::ResetZ),
        ("CY 0 1", GateName::CY),
        ("CZ 0 1", GateName::CZ),
        ("CNOT 0 1", GateName::CNot),
        ("ZCX 0 1", GateName::ZCX),
    ] {
        let p = parse(src).unwrap_or_else(|e| panic!("failed on {src:?}: {e}"));
        match &p.instructions[0] {
            Instruction::Gate(GateOp { name, .. }) => {
                assert_eq!(*name, expected, "src = {src:?}");
            }
            other => panic!("src {src:?} produced {other:?}"),
        }
    }
}
