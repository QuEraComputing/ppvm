// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

use stim_parser::prelude::*;

fn approx_eq(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} vs {b}");
}

#[test]
fn parse_tag_bare_ident() {
    let p = parse("S[T] 0").unwrap();
    match &p.instructions[0] {
        RawInstruction::Gate { name, tags, .. } => {
            assert_eq!(*name, GateName::S);
            assert_eq!(tags.len(), 1);
            assert_eq!(tags[0].name, "T");
            assert!(tags[0].params.is_empty());
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_tag_named_param_pi_expr() {
    let p = parse("I[R_X(theta=0.5*pi)] 0").unwrap();
    match &p.instructions[0] {
        RawInstruction::Gate { name, tags, .. } => {
            assert_eq!(*name, GateName::Identity);
            assert_eq!(tags.len(), 1);
            assert_eq!(tags[0].name, "R_X");
            match &tags[0].params[..] {
                [TagParam::Named { key, value }] => {
                    assert_eq!(key, "theta");
                    approx_eq(*value, 0.5 * std::f64::consts::PI);
                }
                other => panic!("{other:?}"),
            }
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_tag_u3_three_named_params() {
    let p = parse("I[U3(theta=0.34*pi, phi=0.21*pi, lambda=0.46*pi)] 5").unwrap();
    match &p.instructions[0] {
        RawInstruction::Gate { tags, .. } => {
            let params = &tags[0].params;
            assert_eq!(params.len(), 3);
            // Order is preserved.
            for (param, expected_key) in params.iter().zip(["theta", "phi", "lambda"]) {
                let TagParam::Named { key, .. } = param else {
                    panic!("expected named, got {param:?}");
                };
                assert_eq!(key, expected_key);
            }
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_loss_tag_with_args() {
    let p = parse("I_ERROR[loss](1.0) 0").unwrap();
    match &p.instructions[0] {
        RawInstruction::Noise {
            name, tags, args, ..
        } => {
            assert_eq!(*name, NoiseName::IError);
            assert_eq!(
                tags,
                &[Tag {
                    name: "loss".into(),
                    params: vec![]
                }]
            );
            approx_eq(args[0], 1.0);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_correlated_loss_three_args() {
    let p = parse("I_ERROR[correlated_loss](0.1, 0.2, 0.3) 0 1").unwrap();
    match &p.instructions[0] {
        RawInstruction::Noise {
            tags,
            args,
            targets,
            ..
        } => {
            assert_eq!(tags[0].name, "correlated_loss");
            assert_eq!(args.len(), 3);
            assert_eq!(targets, &[0, 1]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_multi_tag() {
    // Stim supports comma-separated multiple tags inside `[…]`.
    let p = parse("S[T,debug] 0").unwrap();
    match &p.instructions[0] {
        RawInstruction::Gate { tags, .. } => {
            assert_eq!(tags.len(), 2);
            assert_eq!(tags[0].name, "T");
            assert_eq!(tags[1].name, "debug");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_tag_positional_floats() {
    // `[R_X(0.5)]` — tag param without a key is positional.
    let p = parse("I[R_X(0.5)] 0").unwrap();
    match &p.instructions[0] {
        RawInstruction::Gate { tags, .. } => match &tags[0].params[..] {
            [TagParam::Positional(v)] => approx_eq(*v, 0.5),
            other => panic!("{other:?}"),
        },
        other => panic!("{other:?}"),
    }
}
