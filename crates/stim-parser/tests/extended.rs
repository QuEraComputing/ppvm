use stim_parser::ast::{AnnotationKind, GateName, MeasureName, NoiseName};
use stim_parser::extended::{
    Axis, ExtendedInstruction, ExtendedParseError, ExtendedProgram, RawPassthrough, parse_extended,
};

fn parse_ok(src: &str) -> ExtendedProgram {
    parse_extended(src).expect("parse_extended")
}

fn parse_err(src: &str) -> ExtendedParseError {
    parse_extended(src).expect_err("must reject")
}

fn approx_eq(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} vs {b}");
}

#[test]
fn vanilla_h_passes_through() {
    let p = parse_ok("H 0\n");
    assert_eq!(p.instructions.len(), 1);
    match &p.instructions[0] {
        ExtendedInstruction::Raw(RawPassthrough::Gate {
            name,
            tags,
            targets,
            line,
            ..
        }) => {
            assert_eq!(*name, GateName::H);
            assert!(tags.is_empty());
            assert_eq!(targets, &vec![0]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn vanilla_measure_passes_through() {
    let p = parse_ok("M 0 1\n");
    match &p.instructions[0] {
        ExtendedInstruction::Raw(RawPassthrough::Measure { name, targets, .. }) => {
            assert_eq!(*name, MeasureName::M);
            assert_eq!(targets, &vec![0, 1]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn vanilla_depolarize1_noise_passes_through() {
    let p = parse_ok("DEPOLARIZE1(0.01) 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Raw(RawPassthrough::Noise { name, args, .. }) => {
            assert_eq!(*name, NoiseName::Depolarize1);
            assert_eq!(args, &vec![0.01]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn vanilla_annotation_passes_through() {
    let p = parse_ok("TICK\n");
    match &p.instructions[0] {
        ExtendedInstruction::Raw(RawPassthrough::Annotation { kind, .. }) => {
            assert_eq!(*kind, AnnotationKind::Tick);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn vanilla_mpad_passes_through_as_bool_bits() {
    let p = parse_ok("MPAD 0 1\n");
    match &p.instructions[0] {
        ExtendedInstruction::MPad { bits, prob, .. } => {
            assert_eq!(bits, &vec![false, true]);
            assert!(prob.is_none());
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn mpad_non_bit_target_errors_in_extended_parser() {
    let err = parse_err("MPAD 0 2 1\n");
    match err {
        ExtendedParseError::InvalidMPadBit { line, index, value } => {
            assert_eq!(line, 1);
            assert_eq!(index, 1);
            assert_eq!(value, 2);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn repeat_recurses_into_body() {
    let p = parse_ok("REPEAT 3 {\n    H 0\n}\n");
    match &p.instructions[0] {
        ExtendedInstruction::Repeat { count, body, .. } => {
            assert_eq!(*count, 3);
            assert_eq!(body.len(), 1);
            assert!(matches!(
                &body[0],
                ExtendedInstruction::Raw(RawPassthrough::Gate {
                    name: GateName::H,
                    ..
                })
            ));
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn repeat_promotes_extended_rotation_in_body() {
    let p = parse_ok("REPEAT 2 { I[R_X(theta=0.25)] 0 }\n");
    match &p.instructions[0] {
        ExtendedInstruction::Repeat { count, body, .. } => {
            assert_eq!(*count, 2);
            assert_eq!(body.len(), 1);
            match &body[0] {
                ExtendedInstruction::Rotation {
                    axis,
                    theta,
                    targets,
                    ..
                } => {
                    assert!(matches!(axis, Axis::X));
                    approx_eq(*theta, 0.25);
                    assert_eq!(targets, &vec![0]);
                }
                other => panic!("{other:?}"),
            }
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn repeat_invalid_extended_tag_in_body_errors() {
    let err = parse_err("REPEAT 2 { I[R_X] 0 }\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn lenient_unknown_tag_on_h_passes_through() {
    let p = parse_ok("H[unrelated] 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Raw(RawPassthrough::Gate { name, tags, .. }) => {
            assert_eq!(*name, GateName::H);
            assert_eq!(tags.len(), 1);
            assert_eq!(tags[0].name, "unrelated");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_error_propagates() {
    let err = parse_extended("FROBNICATE 0\n").unwrap_err();
    assert!(matches!(err, ExtendedParseError::Parse(_)));
}

#[test]
fn axis_enum_has_xyz() {
    let _x = Axis::X;
    let _y = Axis::Y;
    let _z = Axis::Z;
}

#[test]
fn s_t_promotes_to_t() {
    let p = parse_ok("S[T] 0 1\n");
    match &p.instructions[0] {
        ExtendedInstruction::T { targets, line } => {
            assert_eq!(targets, &vec![0, 1]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn s_dag_t_promotes_to_t_dag() {
    let p = parse_ok("S_DAG[T] 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::TDag { targets, line } => {
            assert_eq!(targets, &vec![0]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn s_with_no_tag_is_vanilla_gate() {
    let p = parse_ok("S 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Raw(RawPassthrough::Gate { name, tags, .. }) => {
            assert_eq!(*name, GateName::S);
            assert!(tags.is_empty());
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn s_dag_with_no_tag_is_vanilla_gate() {
    let p = parse_ok("S_DAG 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Raw(RawPassthrough::Gate { name, tags, .. }) => {
            assert_eq!(*name, GateName::SDag);
            assert!(tags.is_empty());
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn s_with_unknown_tag_errors() {
    let err = parse_err("S[X] 0\n");
    match err {
        ExtendedParseError::InvalidTag {
            tag,
            instruction,
            line,
            ..
        } => {
            assert_eq!(tag, "X");
            assert_eq!(instruction, "S");
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn s_dag_with_unknown_tag_errors() {
    let err = parse_err("S_DAG[X] 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn s_with_multiple_tags_errors() {
    let err = parse_err("S[T, X] 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn s_t_with_params_errors() {
    let err = parse_err("S[T(0.5)] 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn s_dag_t_with_params_errors() {
    let err = parse_err("S_DAG[T(0.5)] 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn i_r_x_promotes_to_rotation_x() {
    let p = parse_ok("I[R_X(theta=0.5*pi)] 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Rotation {
            axis,
            theta,
            targets,
            line,
        } => {
            assert!(matches!(axis, Axis::X));
            approx_eq(*theta, 0.5 * std::f64::consts::PI);
            assert_eq!(targets, &vec![0]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_r_y_promotes_to_rotation_y() {
    let p = parse_ok("I[R_Y(theta=0.25)] 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Rotation { axis, theta, .. } => {
            assert!(matches!(axis, Axis::Y));
            approx_eq(*theta, 0.25);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_r_z_promotes_to_rotation_z() {
    let p = parse_ok("I[R_Z(theta=0.1)] 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Rotation { axis, theta, .. } => {
            assert!(matches!(axis, Axis::Z));
            approx_eq(*theta, 0.1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_u3_promotes_to_u3() {
    let p = parse_ok("I[U3(theta=0.1, phi=0.2, lambda=0.3)] 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::U3 {
            theta,
            phi,
            lambda,
            targets,
            line,
        } => {
            approx_eq(*theta, 0.1);
            approx_eq(*phi, 0.2);
            approx_eq(*lambda, 0.3);
            assert_eq!(targets, &vec![0]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_u3_missing_phi_errors() {
    let err = parse_err("I[U3(theta=0.1, lambda=0.2)] 0\n");
    match err {
        ExtendedParseError::InvalidTag {
            tag, instruction, ..
        } => {
            assert_eq!(tag, "U3");
            assert_eq!(instruction, "I");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_with_no_tag_is_vanilla_identity() {
    let p = parse_ok("I 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Raw(RawPassthrough::Gate { name, tags, .. }) => {
            assert_eq!(*name, GateName::Identity);
            assert!(tags.is_empty());
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_r_x_missing_theta_errors() {
    let err = parse_err("I[R_X] 0\n");
    match err {
        ExtendedParseError::InvalidTag {
            tag,
            instruction,
            line,
            ..
        } => {
            assert_eq!(tag, "R_X");
            assert_eq!(instruction, "I");
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_r_x_extra_named_param_errors() {
    let err = parse_err("I[R_X(theta=0.1, bogus=2.0)] 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn i_r_x_duplicate_theta_errors() {
    let err = parse_err("I[R_X(theta=0.1, theta=0.2)] 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn i_r_x_positional_theta_errors() {
    let err = parse_err("I[R_X(0.1)] 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn i_u3_extra_named_param_errors() {
    let err = parse_err("I[U3(theta=0.1, phi=0.2, lambda=0.3, bogus=0.4)] 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn i_u3_duplicate_theta_errors() {
    let err = parse_err("I[U3(theta=0.1, phi=0.2, lambda=0.3, theta=0.4)] 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn i_with_unknown_tag_errors() {
    let err = parse_err("I[wat] 0\n");
    match err {
        ExtendedParseError::InvalidTag {
            tag,
            instruction,
            line,
            ..
        } => {
            assert_eq!(tag, "wat");
            assert_eq!(instruction, "I");
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_with_multiple_rotation_tags_errors() {
    let err = parse_err("I[R_X(theta=0.1), R_Y(theta=0.2)] 0\n");
    match err {
        ExtendedParseError::InvalidTag {
            tag,
            instruction,
            line,
            ..
        } => {
            assert_eq!(tag, "R_X");
            assert_eq!(instruction, "I");
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_error_loss_promotes_to_loss() {
    let p = parse_ok("I_ERROR[loss](0.01) 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Loss { p, targets, line } => {
            approx_eq(*p, 0.01);
            assert_eq!(targets, &vec![0]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_error_correlated_loss_one_arg_expands_and_pairs_targets() {
    let p = parse_ok("I_ERROR[correlated_loss](0.5) 0 1 2 3\n");
    match &p.instructions[0] {
        ExtendedInstruction::CorrelatedLoss { ps, targets, line } => {
            approx_eq(ps[0], 0.5);
            approx_eq(ps[1], 0.0);
            approx_eq(ps[2], 0.0);
            assert_eq!(targets, &vec![(0, 1), (2, 3)]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_error_correlated_loss_three_args() {
    let p = parse_ok("I_ERROR[correlated_loss](0.1, 0.2, 0.3) 0 1\n");
    match &p.instructions[0] {
        ExtendedInstruction::CorrelatedLoss { ps, .. } => {
            approx_eq(ps[0], 0.1);
            approx_eq(ps[1], 0.2);
            approx_eq(ps[2], 0.3);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_error_correlated_loss_single_target_errors() {
    let err = parse_err("I_ERROR[correlated_loss](0.5) 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn i_error_correlated_loss_odd_targets_errors() {
    let err = parse_err("I_ERROR[correlated_loss](0.5) 0 1 2\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn i_error_correlated_loss_two_args_errors() {
    let err = parse_err("I_ERROR[correlated_loss](0.1, 0.2) 0 1\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn i_error_with_no_tag_errors() {
    let err = parse_err("I_ERROR(0.1) 0\n");
    match err {
        ExtendedParseError::InvalidTag { instruction, .. } => {
            assert_eq!(instruction, "I_ERROR");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_error_loss_wrong_arg_count_errors() {
    let err = parse_err("I_ERROR[loss](0.1, 0.2) 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn i_error_unknown_tag_errors() {
    let err = parse_err("I_ERROR[bogus](0.1) 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn i_error_multiple_tags_errors() {
    let err = parse_err("I_ERROR[loss, bogus](0.1) 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn prelude_exposes_parse_extended_and_types() {
    mod prelude_scope {
        use stim_parser::prelude::*;

        pub fn check() {
            let p: ExtendedProgram = parse_extended("H 0\n").unwrap();
            assert_eq!(p.instructions.len(), 1);
            // ExtendedParseError, ExtendedInstruction, Axis must also be in scope.
            fn _is_axis(_: Axis) {}
            fn _is_inst(_: &ExtendedInstruction) {}
            fn _is_err(_: ExtendedParseError) {}
        }
    }

    prelude_scope::check();
}

// ----------------------------------------------------------------
// measurement_count
// ----------------------------------------------------------------

#[test]
fn measurement_count_counts_m_mz_mr() {
    let p = parse_ok("X 0\nM 0 1 2\nMR 5");
    assert_eq!(p.measurement_count(), 4);
}

#[test]
fn measurement_count_includes_repeat_multiplier() {
    let p = parse_ok("REPEAT 10 {\n    X 0\n    M 0 1\n}");
    assert_eq!(p.measurement_count(), 20);
}

#[test]
fn measurement_count_mpad_inside_repeat_block_multiplies() {
    let p = parse_ok("REPEAT 3 {\n    MPAD 1\n}");
    assert_eq!(p.measurement_count(), 3);
}

#[test]
fn measurement_count_nested_repeats_multiply_measure_and_mpad() {
    let p = parse_ok("REPEAT 2 {\n    M 0\n    REPEAT 3 {\n        MPAD 0 1\n    }\n}");
    assert_eq!(p.measurement_count(), 14);
}
