// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Stage 2: table-driven validation. Walks the `RawSyntaxTree` produced by the
//! grammar and emits a typed vanilla [`Program`]. Unlike the reference crate's
//! first-error `Result` returns, every recoverable error is reported to a
//! [`DiagnosticSink`]; the sink's returned [`Flow`] decides whether to abort the
//! whole stage or to skip the offending instruction and keep validating.

use std::sync::Arc;

use crate::ast::shared::{
    AnnotationOp, GateOp, MeasureOp, MppOp, NoiseOp, PauliAxis, PauliFactor, Target,
};
use crate::ast::vanilla::{Instruction, Program};
use crate::diagnostics::{Aborted, DiagnosticSink, LineMap, Span};
use crate::instructions::{ArgCount, EntryKind, MeasureName, TableEntry, TargetArity, lookup};
use crate::syntax::{RawSyntaxNode, RawSyntaxTree, RawTarget};

use super::emit_skip;

/// Walk the raw syntactic tree and build the validated [`Program`].
pub(crate) fn validate(
    tree: RawSyntaxTree,
    line_map: &Arc<LineMap>,
    sink: &mut dyn DiagnosticSink,
) -> Result<Program, Aborted> {
    let instructions = validate_slice(tree, line_map, sink)?;
    Ok(Program {
        instructions,
        line_map: Arc::clone(line_map),
    })
}

/// Validate a (possibly nested) list of raw nodes into instructions. Valid
/// instructions are pushed; recoverable errors are skipped; an aborting sink
/// short-circuits the whole walk.
fn validate_slice(
    nodes: Vec<RawSyntaxNode>,
    line_map: &Arc<LineMap>,
    sink: &mut dyn DiagnosticSink,
) -> Result<Vec<Instruction>, Aborted> {
    let mut out = Vec::with_capacity(nodes.len());
    for node in nodes {
        if let Some(instruction) = validate_node(node, line_map, sink)? {
            out.push(instruction);
        }
    }
    Ok(out)
}

/// Validate a single raw node.
///
/// - `Ok(Some(i))` — a valid instruction.
/// - `Ok(None)` — a recoverable error was emitted; skip this instruction.
/// - `Err(Aborted)` — the sink demanded the stage abort.
fn validate_node(
    node: RawSyntaxNode,
    line_map: &Arc<LineMap>,
    sink: &mut dyn DiagnosticSink,
) -> Result<Option<Instruction>, Aborted> {
    match node {
        RawSyntaxNode::Instruction {
            name,
            tag,
            args,
            targets,
            span,
        } => {
            let span: Span = span.into();
            let Some(entry) = lookup(&name) else {
                return emit_skip(
                    sink,
                    span,
                    "unknown-instruction",
                    format!("unknown instruction '{name}'"),
                );
            };
            let arg_rule = entry.args;
            let target_rule = entry.targets;
            let canonical = entry.canonical;

            // MPP carries Pauli-product targets (`X0*Y1*Z2`), not qubit indices,
            // so it parses on a dedicated path rather than the qubit/record loop.
            if matches!(entry.kind, EntryKind::Measure(MeasureName::MPP)) {
                let Some(products) = parse_mpp_products(&targets, sink)? else {
                    return Ok(None);
                };
                if let Some(expected) = check_arg_count(arg_rule, args.len()) {
                    return emit_skip(
                        sink,
                        span,
                        "arg-count",
                        format!("'{canonical}' expected {expected} args, got {}", args.len()),
                    );
                }
                if products.is_empty() {
                    return emit_skip(
                        sink,
                        span,
                        "target-count",
                        format!("'{canonical}' expected target count divisible by 1, got 0"),
                    );
                }
                return Ok(Some(Instruction::Mpp(MppOp {
                    tag,
                    args,
                    products,
                    span,
                })));
            }

            let is_annotation = matches!(entry.kind, EntryKind::Annotation(_));
            // Only gates accept measurement-record controls (`rec[-k]`), as in
            // `CX rec[-1] 1`. For every other instruction a `rec[...]` target is
            // either tolerated-and-dropped (annotations) or an error.
            let is_gate = matches!(entry.kind, EntryKind::Gate(_));
            let mut parsed_targets: Vec<Target> = Vec::with_capacity(targets.len());
            for t in &targets {
                if let Ok(n) = t.text.parse::<usize>() {
                    parsed_targets.push(Target::Qubit(n));
                    continue;
                }
                if is_gate && let Some(k) = parse_rec(&t.text) {
                    parsed_targets.push(Target::Rec(k));
                    continue;
                }
                if is_annotation {
                    // Annotations (DETECTOR, OBSERVABLE_INCLUDE, …) tolerate
                    // non-numeric targets like `rec[-1]` by dropping them.
                    continue;
                }
                let target_span: Span = t.span.into();
                return emit_skip(
                    sink,
                    target_span,
                    "invalid-target",
                    format!("invalid target {:?}", t.text),
                );
            }

            if let Some(expected) = check_arg_count(arg_rule, args.len()) {
                return emit_skip(
                    sink,
                    span,
                    "arg-count",
                    format!("'{canonical}' expected {expected} args, got {}", args.len()),
                );
            }

            let divisor = match target_rule {
                TargetArity::Any => None,
                TargetArity::AtLeastOne => Some(1),
                TargetArity::Pairs => Some(2),
                TargetArity::Quadruples => Some(4),
            };
            if let Some(d) = divisor {
                let n = parsed_targets.len();
                if n == 0 || !n.is_multiple_of(d) {
                    return emit_skip(
                        sink,
                        span,
                        "target-count",
                        format!("'{canonical}' expected target count divisible by {d}, got {n}"),
                    );
                }
            }

            Ok(Some(build_instruction(
                entry,
                tag,
                args,
                parsed_targets,
                span,
            )))
        }
        RawSyntaxNode::Repeat {
            tag,
            count,
            body,
            span,
        } => {
            let span: Span = span.into();
            let body = validate_slice(body, line_map, sink)?;
            Ok(Some(Instruction::Repeat {
                tag,
                count,
                body,
                span,
            }))
        }
    }
}

/// Validate an instruction's argument count against its table rule. Returns
/// `None` when the count is acceptable, or `Some(expected)` (the count the
/// instruction wanted) when it is not.
fn check_arg_count(arg_rule: ArgCount, found: usize) -> Option<usize> {
    match arg_rule {
        ArgCount::Deferred | ArgCount::Any => None,
        ArgCount::None if found != 0 => Some(0),
        ArgCount::Exact(n) if found != n => Some(n),
        ArgCount::Optional(n) if found != 0 && found != n => Some(n),
        _ => None,
    }
}

/// Parse the space-separated Pauli-product targets of an `MPP` instruction,
/// e.g. `X0*Y1*Z2 Z3*Z4` into two products. Each product is a `*`-joined run
/// of single-qubit Pauli factors (`<X|Y|Z><qubit>`).
///
/// - `Ok(Some(products))` — every target parsed (possibly an empty list when
///   there were no targets; the caller maps that to a `target-count` error).
/// - `Ok(None)` — a bad factor was reported and the sink chose to continue;
///   skip this instruction.
/// - `Err(Aborted)` — the sink demanded the stage abort.
fn parse_mpp_products(
    targets: &[RawTarget],
    sink: &mut dyn DiagnosticSink,
) -> Result<Option<Vec<Vec<PauliFactor>>>, Aborted> {
    let mut products = Vec::with_capacity(targets.len());
    for t in targets {
        let mut factors = Vec::new();
        for factor in t.text.split('*') {
            let mut chars = factor.chars();
            let axis = match chars.next() {
                Some('X') => PauliAxis::X,
                Some('Y') => PauliAxis::Y,
                Some('Z') => PauliAxis::Z,
                _ => return invalid_mpp_target(t, sink),
            };
            let Ok(qubit) = chars.as_str().parse::<usize>() else {
                return invalid_mpp_target(t, sink);
            };
            factors.push(PauliFactor { axis, qubit });
        }
        products.push(factors);
    }
    Ok(Some(products))
}

/// Report a bad MPP target and translate the sink's decision into the
/// `parse_mpp_products` return shape.
fn invalid_mpp_target(
    t: &RawTarget,
    sink: &mut dyn DiagnosticSink,
) -> Result<Option<Vec<Vec<PauliFactor>>>, Aborted> {
    let span: Span = t.span.into();
    emit_skip(
        sink,
        span,
        "invalid-mpp-target",
        format!("invalid MPP target {:?}", t.text),
    )
}

/// Parse a measurement-record lookback target `rec[-k]` into its lookback
/// distance `k >= 1`. Returns `None` for anything that isn't a well-formed,
/// strictly-negative record reference (so `rec[0]`, `rec[1]`, `rec[]` and
/// non-`rec` text all fall through to the caller's error handling).
fn parse_rec(text: &str) -> Option<usize> {
    let inner = text.strip_prefix("rec[")?.strip_suffix(']')?;
    let magnitude = inner.strip_prefix('-')?;
    match magnitude.parse::<usize>() {
        Ok(k) if k >= 1 => Some(k),
        _ => None,
    }
}

/// Convert validated targets to bare qubit indices. Only gate targets may be
/// `rec[...]`; every other instruction's targets were validated as qubits in
/// [`validate_node`], so this never drops information for them.
fn qubit_indices(targets: Vec<Target>) -> Vec<usize> {
    targets
        .into_iter()
        .map(|t| {
            t.as_qubit()
                .expect("non-gate targets are validated as plain qubits")
        })
        .collect()
}

fn build_instruction(
    entry: TableEntry,
    tag: String,
    args: Vec<f64>,
    targets: Vec<Target>,
    span: Span,
) -> Instruction {
    match entry.kind {
        EntryKind::Gate(name) => Instruction::Gate(GateOp {
            name,
            tag,
            args,
            targets,
            span,
        }),
        EntryKind::Noise(name) => Instruction::Noise(NoiseOp {
            name,
            tag,
            args,
            targets: qubit_indices(targets),
            span,
        }),
        EntryKind::Measure(name) => Instruction::Measure(MeasureOp {
            name,
            tag,
            args,
            targets: qubit_indices(targets),
            span,
        }),
        EntryKind::Annotation(kind) => Instruction::Annotation(AnnotationOp {
            kind,
            tag,
            args,
            targets: qubit_indices(targets),
            span,
        }),
        EntryKind::MPad => Instruction::MPad {
            tag,
            prob: args.into_iter().next(),
            bits: qubit_indices(targets),
            span,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ast::shared::Target;
    use crate::ast::vanilla::Instruction;
    use crate::diagnostics::{Collect, Diagnostic, FailFast, LineMap};
    use crate::instructions::{GateName, MeasureName, NoiseName};
    use crate::syntax::raw::{RawSyntaxNode, RawTarget};
    use chumsky::span::SimpleSpan;
    use std::sync::Arc;

    fn lm() -> Arc<LineMap> {
        Arc::new(LineMap::new("H 0\nM 0"))
    }

    fn instr(
        name: &str,
        args: Vec<f64>,
        targets: Vec<&str>,
        span: (usize, usize),
    ) -> RawSyntaxNode {
        RawSyntaxNode::Instruction {
            name: name.to_string(),
            tag: String::new(),
            args,
            targets: targets
                .into_iter()
                .map(|t| RawTarget {
                    text: t.to_string(),
                    span: SimpleSpan::from(span.0..span.1),
                })
                .collect(),
            span: SimpleSpan::from(span.0..span.1),
        }
    }

    fn instr_with_target_spans(
        name: &str,
        args: Vec<f64>,
        targets: Vec<(&str, (usize, usize))>,
        span: (usize, usize),
    ) -> RawSyntaxNode {
        RawSyntaxNode::Instruction {
            name: name.to_string(),
            tag: String::new(),
            args,
            targets: targets
                .into_iter()
                .map(|(text, span)| RawTarget {
                    text: text.to_string(),
                    span: SimpleSpan::from(span.0..span.1),
                })
                .collect(),
            span: SimpleSpan::from(span.0..span.1),
        }
    }

    fn instr_with_tag(name: &str, tag: &str) -> RawSyntaxNode {
        RawSyntaxNode::Instruction {
            name: name.to_string(),
            tag: tag.to_string(),
            args: vec![],
            targets: vec![RawTarget {
                text: "0".to_string(),
                span: SimpleSpan::from(2..3),
            }],
            span: SimpleSpan::from(0..1),
        }
    }

    /// Validate with a `FailFast` sink, asserting no error was emitted, and
    /// return the produced program.
    fn ok_program(nodes: Vec<RawSyntaxNode>, line_map: &Arc<LineMap>) -> Program {
        let mut sink = FailFast::new();
        let program = validate(nodes, line_map, &mut sink).expect("validation should succeed");
        assert!(!sink.saw_error(), "expected no diagnostics");
        program
    }

    /// Validate with a `Collect` sink (never aborts) and return the collected
    /// diagnostics for inspection. Callers must pass nodes that are ALL
    /// expected to be invalid — every node is emitted-and-skipped, so the
    /// resulting program is empty (asserted here as a guard).
    fn collect_errors(nodes: Vec<RawSyntaxNode>, line_map: &Arc<LineMap>) -> Vec<Diagnostic> {
        let mut sink = Collect::new();
        let program = validate(nodes, line_map, &mut sink).expect("Collect never aborts");
        assert!(
            program.instructions.is_empty(),
            "collect_errors expects every input node to be invalid"
        );
        sink.into_items()
    }

    #[test]
    fn validates_simple_gate() {
        let prog = ok_program(vec![instr("H", vec![], vec!["0"], (0, 1))], &lm());
        assert!(matches!(
            &prog.instructions[0],
            Instruction::Gate(GateOp {
                name: GateName::H,
                ..
            })
        ));
        // Span resolves to line 1.
        match &prog.instructions[0] {
            Instruction::Gate(op) => assert_eq!(op.span.line(&prog.line_map), 1),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn validates_measure() {
        let prog = ok_program(vec![instr("M", vec![], vec!["0"], (4, 5))], &lm());
        match &prog.instructions[0] {
            Instruction::Measure(op) => {
                assert_eq!(op.name, MeasureName::M);
                assert_eq!(op.span.line(&prog.line_map), 2);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn unknown_instruction_errors() {
        let items = collect_errors(vec![instr("FROBNICATE", vec![], vec!["0"], (0, 10))], &lm());
        assert_eq!(items[0].code, Some("unknown-instruction"));
        assert_eq!(items[0].message, "unknown instruction 'FROBNICATE'");
        assert_eq!(items[0].span.line(&lm()), 1);
    }

    #[test]
    fn arg_count_errors() {
        // DEPOLARIZE1 expects exactly 1 arg.
        let items = collect_errors(
            vec![instr("DEPOLARIZE1", vec![0.1, 0.2], vec!["0"], (0, 11))],
            &lm(),
        );
        assert_eq!(items[0].code, Some("arg-count"));
    }

    #[test]
    fn target_pair_errors() {
        let items = collect_errors(
            vec![instr("CX", vec![], vec!["0", "1", "2"], (0, 2))],
            &lm(),
        );
        assert_eq!(items[0].code, Some("target-count"));
        assert_eq!(
            items[0].message,
            "'CX' expected target count divisible by 2, got 3"
        );
    }

    #[test]
    fn invalid_target_for_gate_is_invalid_target() {
        let items = collect_errors(vec![instr("H", vec![], vec!["abc"], (0, 1))], &lm());
        assert_eq!(items[0].code, Some("invalid-target"));
    }

    #[test]
    fn annotation_tolerates_non_numeric_targets() {
        // DETECTOR is an annotation: rec[-1] should be silently dropped.
        let prog = ok_program(
            vec![instr("DETECTOR", vec![], vec!["rec[-1]"], (0, 8))],
            &lm(),
        );
        assert!(matches!(&prog.instructions[0], Instruction::Annotation(_)));
    }

    #[test]
    fn invalid_target_takes_precedence_over_arg_count_error() {
        // The bad target is reported before the arg-count check runs.
        let items = collect_errors(
            vec![instr("DEPOLARIZE1", vec![0.1, 0.2], vec!["abc"], (0, 11))],
            &lm(),
        );
        assert_eq!(items[0].code, Some("invalid-target"));
    }

    #[test]
    fn invalid_target_uses_target_span_for_line_col() {
        let line_map = Arc::new(LineMap::new("X 0\nH 0 abc"));
        let nodes = vec![instr_with_target_spans(
            "H",
            vec![],
            vec![("0", (6, 7)), ("abc", (8, 11))],
            (4, 5),
        )];
        let mut sink = Collect::new();
        validate(nodes, &line_map, &mut sink).expect("Collect never aborts");
        let items = sink.into_items();
        assert_eq!(items[0].code, Some("invalid-target"));
        // The diagnostic's span is the bad target's span (8..11) → line 2, col 5.
        assert_eq!(items[0].span.line_col(&line_map), (2, 5));
        assert_eq!(items[0].message, "invalid target \"abc\"");
    }

    #[test]
    fn repeat_body_is_validated_recursively() {
        let nodes = vec![RawSyntaxNode::Repeat {
            tag: String::new(),
            count: 3,
            body: vec![instr("H", vec![], vec!["0"], (11, 12))],
            span: SimpleSpan::from(0..6),
        }];
        let line_map = Arc::new(LineMap::new("REPEAT 3 { H 0 }"));
        let prog = ok_program(nodes, &line_map);
        match &prog.instructions[0] {
            Instruction::Repeat {
                count, body, span, ..
            } => {
                assert_eq!(*count, 3);
                assert_eq!(span.line(&prog.line_map), 1);
                assert!(matches!(
                    &body[0],
                    Instruction::Gate(GateOp {
                        name: GateName::H,
                        targets,
                        ..
                    }) if targets == &vec![Target::Qubit(0)]
                ));
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn tag_passes_through_validator() {
        let nodes = vec![instr_with_tag("H", "R(0.5, theta=0.25)")];
        let prog = ok_program(nodes, &lm());
        match &prog.instructions[0] {
            Instruction::Gate(GateOp { tag, .. }) => assert_eq!(tag, "R(0.5, theta=0.25)"),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn at_least_one_target_errors_on_zero_targets() {
        let items = collect_errors(vec![instr("H", vec![], vec![], (0, 1))], &lm());
        assert_eq!(items[0].code, Some("target-count"));
        assert_eq!(
            items[0].message,
            "'H' expected target count divisible by 1, got 0"
        );
    }

    #[test]
    fn optional_measure_arg_accepts_zero_or_one_arg() {
        let no_arg = ok_program(vec![instr("M", vec![], vec!["0"], (0, 1))], &lm());
        assert!(matches!(&no_arg.instructions[0], Instruction::Measure(_)));

        let one_arg = ok_program(vec![instr("M", vec![0.25], vec!["0"], (0, 1))], &lm());
        match &one_arg.instructions[0] {
            Instruction::Measure(op) => assert_eq!(op.args, vec![0.25]),
            other => panic!("{other:?}"),
        }

        let items = collect_errors(vec![instr("M", vec![0.1, 0.2], vec!["0"], (0, 1))], &lm());
        assert_eq!(items[0].code, Some("arg-count"));
        assert_eq!(items[0].message, "'M' expected 1 args, got 2");
    }

    #[test]
    fn i_error_arg_count_is_deferred() {
        // I_ERROR's arg count is Deferred, so 2 args is accepted at this stage.
        let prog = ok_program(
            vec![instr("I_ERROR", vec![0.1, 0.2], vec!["0"], (0, 7))],
            &lm(),
        );
        match &prog.instructions[0] {
            Instruction::Noise(op) => {
                assert_eq!(op.name, NoiseName::IError);
                assert_eq!(op.args, vec![0.1, 0.2]);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn mpp_happy_path() {
        // MPP X0*Y1*Z2 Z3 → two products: [X0,Y1,Z2] and [Z3].
        let line_map = Arc::new(LineMap::new("MPP X0*Y1*Z2 Z3"));
        let nodes = vec![instr_with_target_spans(
            "MPP",
            vec![],
            vec![("X0*Y1*Z2", (4, 12)), ("Z3", (13, 15))],
            (0, 3),
        )];
        let prog = ok_program(nodes, &line_map);
        match &prog.instructions[0] {
            Instruction::Mpp(op) => {
                assert_eq!(op.products.len(), 2);
                assert_eq!(
                    op.products[0],
                    vec![
                        PauliFactor {
                            axis: PauliAxis::X,
                            qubit: 0
                        },
                        PauliFactor {
                            axis: PauliAxis::Y,
                            qubit: 1
                        },
                        PauliFactor {
                            axis: PauliAxis::Z,
                            qubit: 2
                        },
                    ]
                );
                assert_eq!(
                    op.products[1],
                    vec![PauliFactor {
                        axis: PauliAxis::Z,
                        qubit: 3
                    }]
                );
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn mpp_bad_factor_errors() {
        // `W0` has an invalid Pauli axis.
        let line_map = Arc::new(LineMap::new("MPP W0"));
        let nodes = vec![instr_with_target_spans(
            "MPP",
            vec![],
            vec![("W0", (4, 6))],
            (0, 3),
        )];
        let mut sink = Collect::new();
        validate(nodes, &line_map, &mut sink).expect("Collect never aborts");
        let items = sink.into_items();
        assert_eq!(items[0].code, Some("invalid-mpp-target"));
        assert_eq!(items[0].message, "invalid MPP target \"W0\"");
        assert_eq!(items[0].span.line_col(&line_map), (1, 5));
    }

    #[test]
    fn mpp_empty_targets_is_target_count_error() {
        let nodes = vec![instr("MPP", vec![], vec![], (0, 3))];
        let items = collect_errors(nodes, &lm());
        assert_eq!(items[0].code, Some("target-count"));
        assert_eq!(
            items[0].message,
            "'MPP' expected target count divisible by 1, got 0"
        );
    }

    #[test]
    fn fail_fast_aborts_and_recoverable_skip_continues() {
        // FailFast on an error aborts the stage.
        let mut sink = FailFast::new();
        let res = validate(
            vec![instr("FROBNICATE", vec![], vec!["0"], (0, 10))],
            &lm(),
            &mut sink,
        );
        assert!(res.is_err());
        assert!(sink.saw_error());

        // Collect skips the bad instruction but keeps the good one.
        let mut sink = Collect::new();
        let prog = validate(
            vec![
                instr("FROBNICATE", vec![], vec!["0"], (0, 10)),
                instr("H", vec![], vec!["0"], (0, 1)),
            ],
            &lm(),
            &mut sink,
        )
        .expect("Collect never aborts");
        assert_eq!(prog.instructions.len(), 1);
        assert!(matches!(
            &prog.instructions[0],
            Instruction::Gate(GateOp {
                name: GateName::H,
                ..
            })
        ));
        assert_eq!(sink.into_items().len(), 1);
    }
}
