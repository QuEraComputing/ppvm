// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use stim_parser::prelude::*;

fn approx_eq(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} vs {b}");
}

#[test]
fn parse_tag_bare_ident() {
    let p = parse("S[T] 0").unwrap();
    match &p.instructions[0] {
        Instruction::Gate(GateOp { name, tag, .. }) => {
            assert_eq!(*name, GateName::S);
            assert_eq!(tag, "T");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn empty_and_absent_tags_are_empty_strings() {
    let p = parse("H[] 0\nH 1").unwrap();
    for instruction in &p.instructions {
        let Instruction::Gate(GateOp { tag, .. }) = instruction else {
            panic!("{instruction:?}");
        };
        assert!(tag.is_empty());
    }
    assert_eq!(p.to_string(), "H 0\nH 1\n");
}

#[test]
fn parse_tag_named_param_pi_expr() {
    for (src, expected) in [
        ("I[R_X(theta=0.5*pi)] 0", "R_X(theta=0.5*pi)"),
        ("I[R_X(theta=0.5 * pi)] 0", "R_X(theta=0.5 * pi)"),
        ("I[R_X(theta=0.5pi)] 0", "R_X(theta=0.5pi)"),
    ] {
        let p = parse(src).unwrap_or_else(|e| panic!("{src}: {e:?}"));
        match &p.instructions[0] {
            Instruction::Gate(GateOp { name, tag, .. }) => {
                assert_eq!(*name, GateName::Identity);
                assert_eq!(tag, expected);
            }
            other => panic!("{src}: {other:?}"),
        }
    }
}

#[test]
fn parse_tag_u3_three_named_params() {
    let p = parse("I[U3(theta=0.34*pi, phi=0.21*pi, lambda=0.46*pi)] 5").unwrap();
    match &p.instructions[0] {
        Instruction::Gate(GateOp { tag, .. }) => {
            assert_eq!(tag, "U3(theta=0.34*pi, phi=0.21*pi, lambda=0.46*pi)")
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_loss_tag_with_args() {
    let p = parse("I_ERROR[loss](1.0) 0").unwrap();
    match &p.instructions[0] {
        Instruction::Noise(NoiseOp {
            name, tag, args, ..
        }) => {
            assert_eq!(*name, NoiseName::IError);
            assert_eq!(tag, "loss");
            approx_eq(args[0], 1.0);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_correlated_loss_three_args() {
    let p = parse("I_ERROR[correlated_loss](0.1, 0.2, 0.3) 0 1").unwrap();
    match &p.instructions[0] {
        Instruction::Noise(NoiseOp {
            tag, args, targets, ..
        }) => {
            assert_eq!(tag, "correlated_loss");
            assert_eq!(args.len(), 3);
            assert_eq!(targets, &[0, 1]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn commas_are_tag_content() {
    let p = parse("S[T,debug] 0").unwrap();
    match &p.instructions[0] {
        Instruction::Gate(GateOp { tag, .. }) => assert_eq!(tag, "T,debug"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn structured_looking_tag_remains_opaque() {
    let p = parse("I[R_X(0.5)] 0").unwrap();
    match &p.instructions[0] {
        Instruction::Gate(GateOp { tag, .. }) => assert_eq!(tag, "R_X(0.5)"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn arbitrary_stim_tag_is_one_opaque_string() {
    let p = parse("H[x^0*y^0,horizontal,(note),a[b,λ] 0").unwrap();
    assert_eq!(p.to_string(), "H[x^0*y^0,horizontal,(note),a[b,λ] 0\n");
}

#[test]
fn annotation_and_repeat_tags_are_preserved() {
    let p = parse("REPEAT[loop note] 2 { TICK[SE reset] }").unwrap();
    assert_eq!(
        p.to_string(),
        "REPEAT[loop note] 2 {\n    TICK[SE reset]\n}\n"
    );
}

#[test]
fn stim_tag_escapes_decode_and_print_canonically() {
    let p = parse(r"H[test \B\C\r\n] 0").unwrap();
    assert_eq!(p.to_string(), "H[test \\B\\C\\r\\n] 0\n");
}

#[test]
fn opaque_tags_round_trip_on_every_instruction_family() {
    let encoded = r"meta, [λ\B\C\n\r";
    let decoded = "meta, [λ\\]\n\r";
    let src = format!(
        "H[{encoded}] 0\n\
         X_ERROR[{encoded}](0.1) 0\n\
         M[{encoded}] 0\n\
         TICK[{encoded}]\n\
         MPP[{encoded}] X0\n\
         MPAD[{encoded}] 0\n\
         REPEAT[{encoded}] 1 {{ H 0 }}"
    );

    let parsed = parse(&src).unwrap();
    let tags = |program: &Program| {
        program
            .instructions
            .iter()
            .map(|instruction| match instruction {
                Instruction::Gate(op) => &op.tag,
                Instruction::Noise(op) => &op.tag,
                Instruction::Measure(op) => &op.tag,
                Instruction::Annotation(op) => &op.tag,
                Instruction::Mpp(op) => &op.tag,
                Instruction::MPad { tag, .. } | Instruction::Repeat { tag, .. } => tag,
            })
            .cloned()
            .collect::<Vec<_>>()
    };
    assert_eq!(tags(&parsed), vec![decoded.to_string(); 7]);

    let printed = parsed.to_string();
    assert_eq!(printed.matches(&format!("[{encoded}]")).count(), 7);
    assert_eq!(
        tags(&parse(&printed).unwrap()),
        vec![decoded.to_string(); 7]
    );
}
