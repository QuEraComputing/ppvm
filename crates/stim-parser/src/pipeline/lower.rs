// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Stage 3: lowering. Promotes the vanilla [`Program`]'s tag-based PPVM
//! extensions to first-class [`ExtendedInstruction`] variants. Unlike the
//! reference crate's first-error `Result` returns, every recoverable error is
//! reported to a [`DiagnosticSink`]; the sink's returned [`Flow`] decides
//! whether to abort the whole stage or to skip the offending instruction and
//! keep lowering.

use std::sync::Arc;

use chumsky::error::Rich;
use chumsky::extra;
use chumsky::prelude::*;

use crate::ast::extended::{ExtendedInstruction, ExtendedProgram};
use crate::ast::shared::{Axis, GateOp, NoiseOp, Target};
use crate::ast::vanilla::{Instruction, Program};
use crate::diagnostics::{Aborted, DiagnosticSink, Span};
use crate::instructions::{GateName, NoiseName};

use super::emit_skip;

type TagExtra<'src> = extra::Err<Rich<'src, char>>;

#[derive(Debug, Clone, PartialEq)]
struct StructuredTag {
    name: String,
    params: Vec<StructuredTagParam>,
}

#[derive(Debug, Clone, PartialEq)]
enum StructuredTagParam {
    Positional(f64),
    Named {
        key: String,
        value: f64,
        had_pi: bool,
    },
}

fn structured_tag<'src>() -> impl Parser<'src, &'src str, StructuredTag, TagExtra<'src>> + Clone {
    use crate::syntax::{ident, inline_pad, pi_expr, pi_expr_flagged};

    let named = ident()
        .then_ignore(inline_pad())
        .then_ignore(just('='))
        .then_ignore(inline_pad())
        .then(pi_expr_flagged())
        .map(|(key, (value, had_pi))| StructuredTagParam::Named { key, value, had_pi });
    let positional = pi_expr().map(StructuredTagParam::Positional);
    let params = choice((named, positional))
        .separated_by(inline_pad().then(just(',')).then(inline_pad()))
        .allow_trailing()
        .collect::<Vec<_>>()
        .delimited_by(just('(').then(inline_pad()), inline_pad().then(just(')')));

    inline_pad()
        .ignore_then(ident().then(params.or_not()))
        .then_ignore(inline_pad())
        .map(|(name, params)| StructuredTag {
            name,
            params: params.unwrap_or_default(),
        })
}

fn parse_structured_tag(src: &str) -> Option<StructuredTag> {
    structured_tag()
        .then_ignore(end())
        .parse(src)
        .into_result()
        .ok()
}

/// Lower a vanilla [`Program`] into an [`ExtendedProgram`], forwarding every
/// recoverable error to the sink.
pub(crate) fn lower(
    program: Program,
    sink: &mut dyn DiagnosticSink,
) -> Result<ExtendedProgram, Aborted> {
    let line_map = Arc::clone(&program.line_map);
    let instructions = lower_slice(program.instructions, sink)?;
    Ok(ExtendedProgram {
        instructions,
        line_map,
    })
}

/// Lower a (possibly nested) list of instructions. Lowered instructions are
/// pushed; recoverable errors are skipped; an aborting sink short-circuits the
/// whole walk.
fn lower_slice(
    src: Vec<Instruction>,
    sink: &mut dyn DiagnosticSink,
) -> Result<Vec<ExtendedInstruction>, Aborted> {
    let mut out = Vec::with_capacity(src.len());
    for instr in src {
        if let Some(lowered) = lower_one(instr, sink)? {
            out.push(lowered);
        }
    }
    Ok(out)
}

/// Lower a single instruction.
///
/// - `Ok(Some(i))` — a lowered instruction.
/// - `Ok(None)` — a recoverable error was emitted; skip this instruction.
/// - `Err(Aborted)` — the sink demanded the stage abort.
fn lower_one(
    instr: Instruction,
    sink: &mut dyn DiagnosticSink,
) -> Result<Option<ExtendedInstruction>, Aborted> {
    match instr {
        Instruction::Gate(op) => lower_gate(op, sink),
        Instruction::Noise(op) => lower_noise(op, sink),
        // Pass-through families move the shared `*Op` struct straight through.
        Instruction::Measure(op) => Ok(Some(ExtendedInstruction::Measure(op))),
        Instruction::Annotation(op) => Ok(Some(ExtendedInstruction::Annotation(op))),
        Instruction::Mpp(op) => Ok(Some(ExtendedInstruction::Mpp(op))),
        Instruction::MPad {
            tag,
            prob,
            bits,
            span,
        } => {
            let Some(bits) = convert_mpad_bits(&bits, span, sink)? else {
                return Ok(None);
            };
            Ok(Some(ExtendedInstruction::MPad {
                tag,
                prob,
                bits,
                span,
            }))
        }
        Instruction::Repeat {
            tag,
            count,
            body,
            span,
        } => {
            let body = lower_slice(body, sink)?;
            Ok(Some(ExtendedInstruction::Repeat {
                tag,
                count,
                body,
                span,
            }))
        }
    }
}

fn lower_gate(
    op: GateOp,
    sink: &mut dyn DiagnosticSink,
) -> Result<Option<ExtendedInstruction>, Aborted> {
    use GateName::*;

    let GateOp {
        name,
        tag,
        args,
        targets,
        span,
    } = op;

    match name {
        // Tags on native T/T_DAG do not affect simulation.
        T | TDag => {
            let Some(targets) = qubit_targets(targets, name.canonical_name(), span, sink)? else {
                return Ok(None);
            };
            Ok(Some(if matches!(name, GateName::T) {
                ExtendedInstruction::T { targets, span }
            } else {
                ExtendedInstruction::TDag { targets, span }
            }))
        }
        S | SDag => {
            if tag == "T" {
                let Some(targets) = qubit_targets(targets, name.canonical_name(), span, sink)?
                else {
                    return Ok(None);
                };
                return Ok(Some(if matches!(name, S) {
                    ExtendedInstruction::T { targets, span }
                } else {
                    ExtendedInstruction::TDag { targets, span }
                }));
            }
            Ok(Some(ExtendedInstruction::Gate(GateOp {
                name,
                tag,
                args,
                targets,
                span,
            })))
        }
        Identity if tag.is_empty() => Ok(Some(ExtendedInstruction::Gate(GateOp {
            name,
            tag,
            args,
            targets,
            span,
        }))),
        Identity => {
            let Some(targets) = qubit_targets(targets, "I", span, sink)? else {
                return Ok(None);
            };
            interpret_identity_tag(&tag, targets, span, sink)
        }
        _ => Ok(Some(ExtendedInstruction::Gate(GateOp {
            name,
            tag,
            args,
            targets,
            span,
        }))),
    }
}

fn interpret_identity_tag(
    tag_text: &str,
    targets: Vec<usize>,
    span: Span,
    sink: &mut dyn DiagnosticSink,
) -> Result<Option<ExtendedInstruction>, Aborted> {
    let Some(tag) = parse_structured_tag(tag_text) else {
        return invalid_tag(tag_text, "I", span, "expected a PPVM modifier", sink);
    };
    let axis = match tag.name.as_str() {
        "R_X" => Some(Axis::X),
        "R_Y" => Some(Axis::Y),
        "R_Z" => Some(Axis::Z),
        _ => None,
    };

    if let Some(axis) = axis {
        let Some([theta]) = exact_named_params(&tag, ["theta"], "I", span, sink)? else {
            return Ok(None);
        };
        return Ok(Some(ExtendedInstruction::Rotation {
            axis,
            theta,
            targets,
            span,
        }));
    }

    if tag.name == "U3" {
        let Some([theta, phi, lambda]) =
            exact_named_params(&tag, ["theta", "phi", "lambda"], "I", span, sink)?
        else {
            return Ok(None);
        };
        return Ok(Some(ExtendedInstruction::U3 {
            theta,
            phi,
            lambda,
            targets,
            span,
        }));
    }

    invalid_tag(
        &tag.name,
        "I",
        span,
        "unrecognized tag (expected R_X / R_Y / R_Z / U3)",
        sink,
    )
}

fn lower_noise(
    op: NoiseOp,
    sink: &mut dyn DiagnosticSink,
) -> Result<Option<ExtendedInstruction>, Aborted> {
    use NoiseName::*;

    let NoiseOp {
        name,
        tag,
        args,
        targets,
        span,
    } = op;

    match (name, tag.as_str()) {
        (IError, "loss") => {
            if args.len() != 1 {
                return invalid_tag(
                    "loss",
                    name.canonical_name(),
                    span,
                    format!("[loss] expects 1 arg, got {}", args.len()),
                    sink,
                );
            }
            Ok(Some(ExtendedInstruction::Loss {
                p: args[0],
                targets,
                span,
            }))
        }
        (IError, "correlated_loss") => {
            if targets.is_empty() || !targets.len().is_multiple_of(2) {
                return invalid_tag(
                    "correlated_loss",
                    name.canonical_name(),
                    span,
                    format!(
                        "[correlated_loss] expects a nonzero even target count, got {}",
                        targets.len()
                    ),
                    sink,
                );
            }
            let ps = match args.len() {
                1 => [args[0], 0.0, 0.0],
                3 => [args[0], args[1], args[2]],
                n => {
                    return invalid_tag(
                        "correlated_loss",
                        name.canonical_name(),
                        span,
                        format!("[correlated_loss] expects 1 or 3 args, got {n}"),
                        sink,
                    );
                }
            };
            Ok(Some(ExtendedInstruction::CorrelatedLoss {
                ps,
                targets: pair_targets(&targets),
                span,
            }))
        }
        (IError, "leakage") => {
            if args.len() != 2 {
                return invalid_tag(
                    "leakage",
                    name.canonical_name(),
                    span,
                    format!("[leakage] expects 2 args, got {}", args.len()),
                    sink,
                );
            }
            Ok(Some(ExtendedInstruction::Leakage {
                p0: args[0],
                p1: args[1],
                targets,
                span,
            }))
        }
        (IError, "") => invalid_tag(
            "<missing>",
            name.canonical_name(),
            span,
            "I_ERROR requires a [loss], [correlated_loss], or [leakage] tag",
            sink,
        ),
        (IError, other) => invalid_tag(
            other,
            name.canonical_name(),
            span,
            "expected [loss], [correlated_loss], or [leakage]",
            sink,
        ),
        _ => Ok(Some(ExtendedInstruction::Noise(NoiseOp {
            name,
            tag,
            args,
            targets,
            span,
        }))),
    }
}

/// Emit an `invalid-tag` diagnostic and translate the sink's decision into the
/// skip/abort return shape.
fn invalid_tag<T>(
    tag_name: &str,
    instruction: &str,
    span: Span,
    message: impl Into<String>,
    sink: &mut dyn DiagnosticSink,
) -> Result<Option<T>, Aborted> {
    let message = message.into();
    emit_skip(
        sink,
        span,
        "invalid-tag",
        format!("invalid tag {tag_name:?} on {instruction}: {message}"),
    )
}

/// Convert MPAD bit literals (`0`/`1`) to booleans. A bit outside `{0, 1}` is
/// reported as `invalid-mpad-bit`.
fn convert_mpad_bits(
    bits: &[usize],
    span: Span,
    sink: &mut dyn DiagnosticSink,
) -> Result<Option<Vec<bool>>, Aborted> {
    let mut out = Vec::with_capacity(bits.len());
    for (index, value) in bits.iter().copied().enumerate() {
        match value {
            0 => out.push(false),
            1 => out.push(true),
            _ => {
                return emit_skip(
                    sink,
                    span,
                    "invalid-mpad-bit",
                    format!("MPAD bit {index} must be 0 or 1, got {value}"),
                );
            }
        }
    }
    Ok(Some(out))
}

/// Lower gate targets to bare qubit indices for the extended-dialect sugar
/// variants (`T`, rotations, `U3`). Those gates only ever take qubit targets;
/// only the controlled Clifford gates carry record controls (and pass through
/// unchanged). The grammar still accepts `rec[-k]` on any gate, so reject a
/// record target here rather than panicking.
///
/// - `Ok(Some(qubits))` — every target was a plain qubit.
/// - `Ok(None)` — a record target was reported and the sink chose to continue.
/// - `Err(Aborted)` — the sink demanded the stage abort.
fn qubit_targets(
    targets: Vec<Target>,
    instruction: &str,
    span: Span,
    sink: &mut dyn DiagnosticSink,
) -> Result<Option<Vec<usize>>, Aborted> {
    let mut out = Vec::with_capacity(targets.len());
    for t in targets {
        match t.as_qubit() {
            Some(q) => out.push(q),
            None => {
                return emit_skip(
                    sink,
                    span,
                    "record-target-not-allowed",
                    format!("record target rec[-k] not allowed on {instruction}"),
                );
            }
        }
    }
    Ok(Some(out))
}

fn pair_targets(targets: &[usize]) -> Vec<(usize, usize)> {
    targets
        .chunks_exact(2)
        .map(|pair| (pair[0], pair[1]))
        .collect()
}

/// Validate that a tag carries exactly the `required` named parameters — no
/// positional params, no unexpected/duplicate/missing keys — and return their
/// values in `required` order. Any violation is reported as `invalid-tag`.
fn exact_named_params<const N: usize>(
    tag: &StructuredTag,
    required: [&str; N],
    instruction: &str,
    span: Span,
    sink: &mut dyn DiagnosticSink,
) -> Result<Option<[f64; N]>, Aborted> {
    let mut values = [0.0; N];
    let mut seen = [false; N];

    for param in &tag.params {
        match param {
            StructuredTagParam::Positional(_) => {
                return invalid_tag(
                    &tag.name,
                    instruction,
                    span,
                    "tag parameters must be named",
                    sink,
                );
            }
            StructuredTagParam::Named { key, value, had_pi } => {
                let Some(index) = required.iter().position(|required_key| key == required_key)
                else {
                    return invalid_tag(
                        &tag.name,
                        instruction,
                        span,
                        format!("unexpected named parameter '{key}'"),
                        sink,
                    );
                };
                if seen[index] {
                    return invalid_tag(
                        &tag.name,
                        instruction,
                        span,
                        format!("duplicate named parameter '{key}'"),
                        sink,
                    );
                }
                // Rotation/U3 angles are in half-turns: require the `<n>*pi` form,
                // mirroring tsim, so a bare number can't be mistaken for radians.
                if !had_pi {
                    return invalid_tag(
                        &tag.name,
                        instruction,
                        span,
                        format!("parameter '{key}' must be written as <n>*pi (half-turns)"),
                        sink,
                    );
                }
                seen[index] = true;
                values[index] = *value;
            }
        }
    }

    for (index, key) in required.iter().enumerate() {
        if !seen[index] {
            return invalid_tag(
                &tag.name,
                instruction,
                span,
                format!("missing required named parameter '{key}'"),
                sink,
            );
        }
    }

    Ok(Some(values))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ast::shared::Axis;

    /// Drive the whole pipeline (parse → validate → lower) with a `FailFast`
    /// sink, returning the lowered program or the diagnostics that aborted it.
    fn lower_extended(src: &str) -> Result<ExtendedProgram, Vec<crate::diagnostics::Diagnostic>> {
        use crate::diagnostics::FailFast;
        use crate::pipeline::Pipeline;
        let mut sink = FailFast::new();
        Pipeline::new(src)
            .parse(&mut sink)
            .and_then(|p| p.validate(&mut sink))
            .and_then(|p| p.lower(&mut sink))
            .map(|p| p.finish())
            .map_err(|_| sink.into_items())
    }

    #[test]
    fn s_t_tag_lowers_to_t() {
        let prog = lower_extended("S[T] 0").expect("lower");
        assert!(matches!(
            &prog.instructions[0],
            ExtendedInstruction::T { targets, .. } if targets == &vec![0]
        ));
    }

    #[test]
    fn native_t_lowers_to_t() {
        let prog = lower_extended("T 0").expect("lower");
        assert!(matches!(
            &prog.instructions[0],
            ExtendedInstruction::T { targets, .. } if targets == &vec![0]
        ));
    }

    #[test]
    fn s_dag_t_tag_lowers_to_tdag() {
        let prog = lower_extended("S_DAG[T] 0").expect("lower");
        assert!(matches!(
            &prog.instructions[0],
            ExtendedInstruction::TDag { targets, .. } if targets == &vec![0]
        ));
    }

    #[test]
    fn native_t_dag_lowers_to_tdag() {
        let prog = lower_extended("T_DAG 0").expect("lower");
        assert!(matches!(
            &prog.instructions[0],
            ExtendedInstruction::TDag { .. }
        ));
    }

    #[test]
    fn s_without_tag_passes_through() {
        let prog = lower_extended("S 0").expect("lower");
        assert!(matches!(
            &prog.instructions[0],
            ExtendedInstruction::Gate(op) if op.name == GateName::S
        ));
    }

    #[test]
    fn identity_rotation_x_lowers() {
        let prog = lower_extended("I[R_X(theta=0.5*pi)] 0").expect("lower");
        match &prog.instructions[0] {
            ExtendedInstruction::Rotation {
                axis,
                theta,
                targets,
                ..
            } => {
                assert_eq!(*axis, Axis::X);
                assert!((*theta - 0.5 * std::f64::consts::PI).abs() < 1e-12);
                assert_eq!(targets, &vec![0]);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn identity_u3_lowers() {
        let prog =
            lower_extended("I[U3(theta=0.5*pi, phi=1.0*pi, lambda=1.5*pi)] 0").expect("lower");
        match &prog.instructions[0] {
            ExtendedInstruction::U3 {
                theta,
                phi,
                lambda,
                targets,
                ..
            } => {
                let pi = std::f64::consts::PI;
                assert!((*theta - 0.5 * pi).abs() < 1e-12);
                assert!((*phi - 1.0 * pi).abs() < 1e-12);
                assert!((*lambda - 1.5 * pi).abs() < 1e-12);
                assert_eq!(targets, &vec![0]);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn bare_t_ignores_tag() {
        let prog = lower_extended("T[foo] 0").expect("lower");
        assert!(matches!(
            &prog.instructions[0],
            ExtendedInstruction::T { targets, .. } if targets == &vec![0]
        ));
    }

    #[test]
    fn bare_t_dag_ignores_tag() {
        let prog = lower_extended("T_DAG[foo] 0").expect("lower");
        assert!(matches!(
            &prog.instructions[0],
            ExtendedInstruction::TDag { targets, .. } if targets == &vec![0]
        ));
    }

    #[test]
    fn rotation_tag_without_pi_is_rejected() {
        // Mirror tsim: rotation tag angles must be written as <n>*pi (half-turns).
        let err = lower_extended("I[R_Z(theta=0.5)] 0").unwrap_err();
        assert_eq!(err.last().unwrap().code, Some("invalid-tag"));
    }

    #[test]
    fn u3_tag_without_pi_is_rejected() {
        let err = lower_extended("I[U3(theta=0.5, phi=1.0, lambda=1.5)] 0").unwrap_err();
        assert_eq!(err.last().unwrap().code, Some("invalid-tag"));
    }

    #[test]
    fn rotation_tag_with_pi_is_accepted() {
        let prog = lower_extended("I[R_Z(theta=0.5*pi)] 0").expect("lower");
        match &prog.instructions[0] {
            ExtendedInstruction::Rotation { axis, theta, .. } => {
                assert_eq!(*axis, Axis::Z);
                assert!((*theta - 0.5 * std::f64::consts::PI).abs() < 1e-12);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn rotation_tag_accepts_alternate_pi_syntax() {
        let half_pi = 0.5 * std::f64::consts::PI;
        for src in [
            "I[R_Z(theta=0.5*pi)] 0",
            "I[R_Z(theta=0.5 * pi)] 0",
            "I[R_Z(theta=0.5pi)] 0",
        ] {
            let prog = lower_extended(src).unwrap_or_else(|e| panic!("{src}: {e:?}"));
            match &prog.instructions[0] {
                ExtendedInstruction::Rotation { axis, theta, .. } => {
                    assert_eq!(*axis, Axis::Z);
                    assert!((*theta - half_pi).abs() < 1e-12, "{src}");
                }
                other => panic!("{src}: {other:?}"),
            }
        }
    }

    #[test]
    fn i_error_loss_lowers() {
        let prog = lower_extended("I_ERROR[loss](0.01) 0").expect("lower");
        match &prog.instructions[0] {
            ExtendedInstruction::Loss { p, targets, .. } => {
                assert_eq!(*p, 0.01);
                assert_eq!(targets, &vec![0]);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn i_error_correlated_loss_lowers() {
        let prog = lower_extended("I_ERROR[correlated_loss](0.1,0.05,0.05) 0 1").expect("lower");
        match &prog.instructions[0] {
            ExtendedInstruction::CorrelatedLoss { ps, targets, .. } => {
                assert_eq!(*ps, [0.1, 0.05, 0.05]);
                assert_eq!(targets, &vec![(0, 1)]);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn i_error_correlated_loss_single_arg_pads_with_zero() {
        let prog = lower_extended("I_ERROR[correlated_loss](0.1) 0 1").expect("lower");
        match &prog.instructions[0] {
            ExtendedInstruction::CorrelatedLoss { ps, .. } => {
                assert_eq!(*ps, [0.1, 0.0, 0.0]);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn i_error_leakage_lowers() {
        let prog = lower_extended("I_ERROR[leakage](0.1, 0.2) 0").expect("lower");
        match &prog.instructions[0] {
            ExtendedInstruction::Leakage {
                p0, p1, targets, ..
            } => {
                assert_eq!(*p0, 0.1);
                assert_eq!(*p1, 0.2);
                assert_eq!(targets, &vec![0]);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn h_passes_through_as_gate() {
        let prog = lower_extended("H 0").expect("lower");
        assert!(matches!(
            &prog.instructions[0],
            ExtendedInstruction::Gate(op) if op.name == GateName::H
        ));
    }

    #[test]
    fn measure_passes_through() {
        let prog = lower_extended("M 0").expect("lower");
        assert!(matches!(
            &prog.instructions[0],
            ExtendedInstruction::Measure(_)
        ));
    }

    #[test]
    fn mpad_bits_convert_to_bools() {
        let prog = lower_extended("MPAD 0 1 0").expect("lower");
        match &prog.instructions[0] {
            ExtendedInstruction::MPad { bits, .. } => {
                assert_eq!(bits, &vec![false, true, false]);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn mpad_bad_bit_is_invalid_mpad_bit() {
        let err = lower_extended("MPAD 0 2 0").unwrap_err();
        assert_eq!(err.last().unwrap().code, Some("invalid-mpad-bit"));
    }

    #[test]
    fn record_target_on_sugar_gate_is_rejected() {
        let err = lower_extended("T rec[-1]").unwrap_err();
        assert_eq!(err.last().unwrap().code, Some("record-target-not-allowed"));
    }

    #[test]
    fn i_error_without_tag_is_invalid_tag() {
        let err = lower_extended("I_ERROR(0.01) 0").unwrap_err();
        assert_eq!(err.last().unwrap().code, Some("invalid-tag"));
    }

    #[test]
    fn unexpected_named_param_is_invalid_tag() {
        let err = lower_extended("I[R_X(phi=0.5)] 0").unwrap_err();
        assert_eq!(err.last().unwrap().code, Some("invalid-tag"));
    }

    #[test]
    fn missing_named_param_is_invalid_tag() {
        let err = lower_extended("I[U3(theta=0.5, phi=1.0)] 0").unwrap_err();
        assert_eq!(err.last().unwrap().code, Some("invalid-tag"));
    }

    #[test]
    fn s_with_unknown_tag_passes_through() {
        let prog = lower_extended("S[BOGUS] 0").expect("lower");
        assert!(matches!(
            &prog.instructions[0],
            ExtendedInstruction::Gate(op) if op.name == GateName::S && op.tag == "BOGUS"
        ));
    }

    #[test]
    fn s_non_exact_t_tag_passes_through() {
        let prog = lower_extended("S[T(theta=0.5)] 0").expect("lower");
        assert!(matches!(
            &prog.instructions[0],
            ExtendedInstruction::Gate(op) if op.name == GateName::S && op.tag == "T(theta=0.5)"
        ));
    }
}
