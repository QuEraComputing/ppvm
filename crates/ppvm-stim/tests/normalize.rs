use ppvm_stim::{
    parse, normalize, GateKind, Instruction, MeasureKind, NoiseKind,
    NormalizeError, TableauProgram,
};

fn norm(src: &str) -> TableauProgram {
    let prog = parse(src).expect("parse");
    normalize::to_tableau(&prog).expect("normalize")
}

fn norm_err(src: &str) -> NormalizeError {
    let prog = parse(src).expect("parse");
    normalize::to_tableau(&prog).expect_err("must reject")
}

fn approx_eq(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} vs {b}");
}

#[test]
fn h_maps_to_gate_h() {
    let p = norm("H 0");
    assert!(matches!(p.instructions[0], Instruction::Gate { kind: GateKind::H, .. }));
}

#[test]
fn s_t_tag_maps_to_gate_t() {
    let p = norm("S[T] 0");
    let Instruction::Gate { kind, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(*kind, GateKind::T);
}

#[test]
fn s_dag_t_tag_maps_to_gate_t_dag() {
    let p = norm("S_DAG[T] 0");
    let Instruction::Gate { kind, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(*kind, GateKind::TDag);
}

#[test]
fn i_rx_tag_maps_to_rx() {
    let p = norm("I[R_X(theta=0.5*pi)] 0");
    let Instruction::Gate { kind, .. } = &p.instructions[0] else { panic!() };
    if let GateKind::RX { theta } = kind {
        approx_eq(*theta, 0.5 * std::f64::consts::PI);
    } else {
        panic!("expected RX, got {kind:?}");
    }
}

#[test]
fn i_u3_tag_maps_to_u3() {
    let p = norm("I[U3(theta=0.1, phi=0.2, lambda=0.3)] 0");
    let Instruction::Gate { kind, .. } = &p.instructions[0] else { panic!() };
    if let GateKind::U3 { theta, phi, lambda } = kind {
        approx_eq(*theta, 0.1);
        approx_eq(*phi, 0.2);
        approx_eq(*lambda, 0.3);
    } else {
        panic!("{kind:?}");
    }
}

#[test]
fn i_error_loss_tag_maps_to_loss() {
    let p = norm("I_ERROR[loss](1.0) 0");
    let Instruction::Noise { kind, args, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(*kind, NoiseKind::Loss);
    approx_eq(args[0], 1.0);
}

#[test]
fn i_error_correlated_loss_one_arg_expands_to_three() {
    let p = norm("I_ERROR[correlated_loss](0.5) 0 1");
    let Instruction::Noise { kind, args, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(*kind, NoiseKind::CorrelatedLoss);
    assert_eq!(args.len(), 3);
    approx_eq(args[0], 0.5);
    approx_eq(args[1], 0.0);
    approx_eq(args[2], 0.0);
}

#[test]
fn i_error_correlated_loss_three_args_passthrough() {
    let p = norm("I_ERROR[correlated_loss](0.1, 0.2, 0.3) 0 1");
    let Instruction::Noise { args, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(args.len(), 3);
    approx_eq(args[0], 0.1);
    approx_eq(args[1], 0.2);
    approx_eq(args[2], 0.3);
}

#[test]
fn cnot_alias_maps_to_cx() {
    let p = norm("CNOT 0 1");
    assert!(matches!(p.instructions[0], Instruction::Gate { kind: GateKind::CX, .. }));
}

#[test]
fn h_xz_alias_maps_to_h() {
    let p = norm("H_XZ 0");
    assert!(matches!(p.instructions[0], Instruction::Gate { kind: GateKind::H, .. }));
}

#[test]
fn sqrt_z_alias_maps_to_s() {
    let p = norm("SQRT_Z 0");
    assert!(matches!(p.instructions[0], Instruction::Gate { kind: GateKind::S, .. }));
}

#[test]
fn r_and_rz_both_map_to_reset() {
    let p = norm("R 0\nRZ 1");
    assert!(matches!(p.instructions[0], Instruction::Gate { kind: GateKind::Reset, .. }));
    assert!(matches!(p.instructions[1], Instruction::Gate { kind: GateKind::Reset, .. }));
}

#[test]
fn x_error_y_error_z_error_supported() {
    for (src, expected) in [
        ("X_ERROR(0.1) 0", NoiseKind::XError),
        ("Y_ERROR(0.1) 0", NoiseKind::YError),
        ("Z_ERROR(0.1) 0", NoiseKind::ZError),
    ] {
        let p = norm(src);
        let Instruction::Noise { kind, .. } = &p.instructions[0] else { panic!() };
        assert_eq!(*kind, expected, "src = {src}");
    }
}

#[test]
fn measurements_m_mz_map_to_m() {
    let p = norm("M 0\nMZ 1");
    assert!(matches!(p.instructions[0], Instruction::Measure { kind: MeasureKind::M, .. }));
    assert!(matches!(p.instructions[1], Instruction::Measure { kind: MeasureKind::M, .. }));
}

#[test]
fn measurement_mr_maps_to_mr() {
    let p = norm("MR 0");
    assert!(matches!(p.instructions[0], Instruction::Measure { kind: MeasureKind::MR, .. }));
}

#[test]
fn annotations_become_no_op_annotations() {
    let p = norm("DETECTOR\nTICK\nQUBIT_COORDS(0,0) 0");
    assert!(matches!(p.instructions[0], Instruction::Annotation { .. }));
    assert!(matches!(p.instructions[1], Instruction::Annotation { .. }));
    assert!(matches!(p.instructions[2], Instruction::Annotation { .. }));
}

#[test]
fn expected_measurement_count_counts_m_mz_mr() {
    let p = norm("X 0\nM 0 1 2\nMR 5");
    assert_eq!(p.expected_measurement_count, 4);
}

#[test]
fn expected_measurement_count_includes_repeat_multiplier() {
    let p = norm("REPEAT 10 {\n    X 0\n    M 0 1\n}");
    assert_eq!(p.expected_measurement_count, 20);
}

#[test]
fn unsupported_swap_rejected() {
    let err = norm_err("SWAP 0 1");
    match err {
        NormalizeError::Unsupported { name, line } => {
            assert_eq!(name, "SWAP");
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn unsupported_mx_rejected() {
    let err = norm_err("MX 0");
    assert!(matches!(err, NormalizeError::Unsupported { .. }));
}

#[test]
fn unsupported_heralded_erase_rejected() {
    let err = norm_err("HERALDED_ERASE(0.1) 0");
    assert!(matches!(err, NormalizeError::Unsupported { .. }));
}

#[test]
fn malformed_rx_tag_missing_theta_rejected() {
    let err = norm_err("I[R_X] 0");
    match err {
        NormalizeError::InvalidTag { tag, instruction, .. } => {
            assert_eq!(tag, "R_X");
            assert_eq!(instruction, "I");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn untagged_i_error_rejected_as_invalid_tag() {
    let err = norm_err("I_ERROR(0.1) 0");
    assert!(matches!(err, NormalizeError::InvalidTag { .. }));
}

#[test]
fn measure_noise_arg_passes_through_normalize() {
    let p = norm("MZ(0.01) 0");
    let Instruction::Measure { kind, noise, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(*kind, MeasureKind::M);
    assert!((*noise - 0.01).abs() < 1e-12);
}

#[test]
fn measure_no_noise_arg_defaults_to_zero() {
    let p = norm("M 0");
    let Instruction::Measure { noise, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(*noise, 0.0);
}

#[test]
fn mr_noise_passes_through() {
    let p = norm("MR(0.5) 0");
    let Instruction::Measure { kind, noise, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(*kind, MeasureKind::MR);
    assert!((*noise - 0.5).abs() < 1e-12);
}
