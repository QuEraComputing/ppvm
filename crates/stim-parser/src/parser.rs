use std::sync::Arc;

use crate::ast::{ArgCount, ParseError, Program, RawInstruction, SyntaxError, Tag, TargetArity};
use crate::grammar;
use crate::line_map::LineMap;
use crate::table::{EntryKind, TableEntry, lookup};

fn syntax(
    line: usize,
    col: usize,
    message: impl Into<String>,
    line_map: &Arc<LineMap>,
) -> ParseError {
    ParseError::Syntax(SyntaxError::synth(line, col, message, Arc::clone(line_map)))
}

use chumsky::span::SimpleSpan;

/// Raw syntactic tree produced by the chumsky grammar before
/// table-driven validation. `pub(crate)` because it is plumbing
/// between `grammar.rs` and the validator post-pass; not part of the
/// public API.
#[derive(Debug, Clone)]
pub(crate) enum RawSyntaxNode {
    Instruction {
        name: String,
        tags: Vec<Tag>,
        args: Vec<f64>,
        targets: Vec<RawTarget>,
        span: SimpleSpan<usize>,
    },
    Repeat {
        count: u64,
        body: Vec<RawSyntaxNode>,
        span: SimpleSpan<usize>,
    },
}

#[derive(Debug, Clone)]
pub(crate) struct RawTarget {
    pub text: String,
    pub span: SimpleSpan<usize>,
}

/// Parse Stim source into a [`Program`].
pub fn parse(src: &str) -> Result<Program, ParseError> {
    use chumsky::Parser;
    let line_map = Arc::new(LineMap::new(src));
    let parse_result = grammar::program_parser().parse(src);
    let nodes = match parse_result.into_result() {
        Ok(nodes) => nodes,
        Err(errors) => {
            let first = errors.into_iter().next().expect("non-empty error list");
            return Err(ParseError::Syntax(SyntaxError::from_rich(
                first,
                Arc::clone(&line_map),
            )));
        }
    };
    let instructions = validate_program(nodes, &line_map)?;
    Ok(Program { instructions })
}

/// Walk the raw syntactic tree and emit validated instructions.
fn validate_program(
    nodes: Vec<RawSyntaxNode>,
    line_map: &Arc<LineMap>,
) -> Result<Vec<RawInstruction>, ParseError> {
    let mut out = Vec::with_capacity(nodes.len());
    for node in nodes {
        out.push(validate_node(node, line_map)?);
    }
    Ok(out)
}

fn validate_node(
    node: RawSyntaxNode,
    line_map: &Arc<LineMap>,
) -> Result<RawInstruction, ParseError> {
    match node {
        RawSyntaxNode::Instruction {
            name,
            tags,
            args,
            targets,
            span,
        } => {
            let line = line_map.line_of(span.start);
            let entry = lookup(&name).ok_or(ParseError::UnknownInstruction {
                name: name.clone(),
                line,
            })?;
            let arg_rule = entry.args;
            let target_rule = entry.targets;
            let canonical = entry.canonical();

            let is_annotation = matches!(entry.kind, EntryKind::Annotation(_));
            let mut numeric_targets = Vec::with_capacity(targets.len());
            for t in &targets {
                match t.text.parse::<usize>() {
                    Ok(n) => numeric_targets.push(n),
                    Err(_) if is_annotation => {}
                    Err(_) => {
                        let (l, c) = line_map.line_col(t.span.start);
                        return Err(syntax(
                            l,
                            c,
                            format!("invalid target {:?}", t.text),
                            line_map,
                        ));
                    }
                }
            }

            let skip_arg_validation = matches!(arg_rule, ArgCount::Deferred | ArgCount::Any);
            if !skip_arg_validation {
                match arg_rule {
                    ArgCount::Any => {}
                    ArgCount::None if !args.is_empty() => {
                        return Err(ParseError::ArgCount {
                            name: canonical.to_string(),
                            expected: 0,
                            found: args.len(),
                            line,
                        });
                    }
                    ArgCount::Exact(n) if args.len() != n => {
                        return Err(ParseError::ArgCount {
                            name: canonical.to_string(),
                            expected: n,
                            found: args.len(),
                            line,
                        });
                    }
                    ArgCount::Optional(n) if !args.is_empty() && args.len() != n => {
                        return Err(ParseError::ArgCount {
                            name: canonical.to_string(),
                            expected: n,
                            found: args.len(),
                            line,
                        });
                    }
                    _ => {}
                }
            }

            let divisor = match target_rule {
                TargetArity::Any => None,
                TargetArity::AtLeastOne => Some(1),
                TargetArity::Pairs => Some(2),
                TargetArity::Quadruples => Some(4),
            };
            if let Some(d) = divisor {
                let n = numeric_targets.len();
                if n == 0 || !n.is_multiple_of(d) {
                    return Err(ParseError::TargetCount {
                        name: canonical.to_string(),
                        divisor: d,
                        found: n,
                        line,
                    });
                }
            }

            Ok(build_instruction(entry, tags, args, numeric_targets, line))
        }
        RawSyntaxNode::Repeat { count, body, span } => {
            let line = line_map.line_of(span.start);
            let body = validate_program(body, line_map)?;
            Ok(RawInstruction::Repeat { count, body, line })
        }
    }
}

fn build_instruction(
    entry: TableEntry,
    tags: Vec<Tag>,
    args: Vec<f64>,
    targets: Vec<usize>,
    line: usize,
) -> RawInstruction {
    match entry.kind {
        EntryKind::Gate(name) => RawInstruction::Gate {
            name,
            tags,
            args,
            targets,
            line,
        },
        EntryKind::Noise(name) => RawInstruction::Noise {
            name,
            tags,
            args,
            targets,
            line,
        },
        EntryKind::Measure(name) => RawInstruction::Measure {
            name,
            tags,
            args,
            targets,
            line,
        },
        EntryKind::Annotation(kind) => RawInstruction::Annotation {
            kind,
            args,
            targets,
            line,
        },
        EntryKind::MPad => RawInstruction::MPad {
            tags,
            prob: args.into_iter().next(),
            bits: targets,
            line,
        },
    }
}

#[cfg(test)]
mod validate_tests {
    use super::*;
    use crate::ast::{GateName, MeasureName, NoiseName, RawInstruction, TagParam};
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
            tags: vec![],
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
            tags: vec![],
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

    fn instr_with_tags(name: &str, tags: Vec<Tag>) -> RawSyntaxNode {
        RawSyntaxNode::Instruction {
            name: name.to_string(),
            tags,
            args: vec![],
            targets: vec![RawTarget {
                text: "0".to_string(),
                span: SimpleSpan::from(2..3),
            }],
            span: SimpleSpan::from(0..1),
        }
    }

    #[test]
    fn validates_simple_gate() {
        let nodes = vec![instr("H", vec![], vec!["0"], (0, 1))];
        let result = validate_program(nodes, &lm()).unwrap();
        assert!(matches!(
            &result[0],
            RawInstruction::Gate {
                name: GateName::H,
                line: 1,
                ..
            }
        ));
    }

    #[test]
    fn validates_measure() {
        let nodes = vec![instr("M", vec![], vec!["0"], (4, 5))];
        let result = validate_program(nodes, &lm()).unwrap();
        assert!(matches!(
            &result[0],
            RawInstruction::Measure {
                name: MeasureName::M,
                line: 2,
                ..
            }
        ));
    }

    #[test]
    fn unknown_instruction_errors() {
        let nodes = vec![instr("FROBNICATE", vec![], vec!["0"], (0, 10))];
        let err = validate_program(nodes, &lm()).unwrap_err();
        match err {
            ParseError::UnknownInstruction { name, line } => {
                assert_eq!(name, "FROBNICATE");
                assert_eq!(line, 1);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn arg_count_errors() {
        // DEPOLARIZE1 expects exactly 1 arg.
        let nodes = vec![instr("DEPOLARIZE1", vec![0.1, 0.2], vec!["0"], (0, 11))];
        let err = validate_program(nodes, &lm()).unwrap_err();
        assert!(matches!(err, ParseError::ArgCount { .. }));
    }

    #[test]
    fn target_pair_errors() {
        let nodes = vec![instr("CX", vec![], vec!["0", "1", "2"], (0, 2))];
        let err = validate_program(nodes, &lm()).unwrap_err();
        assert!(matches!(
            err,
            ParseError::TargetCount {
                divisor: 2,
                found: 3,
                ..
            }
        ));
    }

    #[test]
    fn invalid_target_for_gate_is_syntax_error() {
        let nodes = vec![instr("H", vec![], vec!["abc"], (0, 1))];
        let err = validate_program(nodes, &lm()).unwrap_err();
        assert!(matches!(err, ParseError::Syntax(_)));
    }

    #[test]
    fn annotation_tolerates_non_numeric_targets() {
        // DETECTOR is an annotation: rec[-1] should be silently dropped.
        let nodes = vec![instr("DETECTOR", vec![], vec!["rec[-1]"], (0, 8))];
        let result = validate_program(nodes, &lm()).unwrap();
        assert!(matches!(&result[0], RawInstruction::Annotation { .. }));
    }

    #[test]
    fn invalid_target_takes_precedence_over_arg_count_error() {
        let nodes = vec![instr("DEPOLARIZE1", vec![0.1, 0.2], vec!["abc"], (0, 11))];
        let err = validate_program(nodes, &lm()).unwrap_err();
        assert!(matches!(err, ParseError::Syntax(_)));
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
        let err = validate_program(nodes, &line_map).unwrap_err();
        assert_eq!(
            err.to_string(),
            "syntax error at line 2, col 5: invalid target \"abc\""
        );
    }

    #[test]
    fn repeat_body_is_validated_recursively() {
        let nodes = vec![RawSyntaxNode::Repeat {
            count: 3,
            body: vec![instr("H", vec![], vec!["0"], (11, 12))],
            span: SimpleSpan::from(0..6),
        }];
        let result = validate_program(nodes, &Arc::new(LineMap::new("REPEAT 3 { H 0 }"))).unwrap();
        match &result[0] {
            RawInstruction::Repeat { count, body, line } => {
                assert_eq!(*count, 3);
                assert_eq!(*line, 1);
                assert!(matches!(
                    &body[0],
                    RawInstruction::Gate {
                        name: GateName::H,
                        targets,
                        ..
                    } if targets == &vec![0]
                ));
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn tags_pass_through_validator() {
        let nodes = vec![instr_with_tags(
            "H",
            vec![Tag {
                name: "R".to_string(),
                params: vec![
                    TagParam::Positional(0.5),
                    TagParam::Named {
                        key: "theta".to_string(),
                        value: 0.25,
                    },
                ],
            }],
        )];
        let result = validate_program(nodes, &lm()).unwrap();
        match &result[0] {
            RawInstruction::Gate { tags, .. } => {
                assert_eq!(tags[0].name, "R");
                assert!(matches!(tags[0].params[0], TagParam::Positional(0.5)));
                assert!(matches!(
                    &tags[0].params[1],
                    TagParam::Named { key, value } if key == "theta" && *value == 0.25
                ));
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn at_least_one_target_errors_on_zero_targets() {
        let nodes = vec![instr("H", vec![], vec![], (0, 1))];
        let err = validate_program(nodes, &lm()).unwrap_err();
        assert!(matches!(
            err,
            ParseError::TargetCount {
                divisor: 1,
                found: 0,
                ..
            }
        ));
    }

    #[test]
    fn optional_measure_arg_accepts_zero_or_one_arg() {
        let no_arg = validate_program(vec![instr("M", vec![], vec!["0"], (0, 1))], &lm()).unwrap();
        assert!(matches!(&no_arg[0], RawInstruction::Measure { .. }));

        let one_arg =
            validate_program(vec![instr("M", vec![0.25], vec!["0"], (0, 1))], &lm()).unwrap();
        match &one_arg[0] {
            RawInstruction::Measure { args, .. } => assert_eq!(args, &vec![0.25]),
            other => panic!("{other:?}"),
        }

        let err = validate_program(vec![instr("M", vec![0.1, 0.2], vec!["0"], (0, 1))], &lm())
            .unwrap_err();
        assert!(matches!(err, ParseError::ArgCount { expected: 1, .. }));
    }

    #[test]
    fn i_error_arg_count_is_deferred_to_extended_parser() {
        let nodes = vec![instr("I_ERROR", vec![0.1, 0.2], vec!["0"], (0, 7))];
        let result = validate_program(nodes, &lm()).unwrap();
        assert!(matches!(
            &result[0],
            RawInstruction::Noise {
                name: NoiseName::IError,
                args,
                ..
            } if args == &vec![0.1, 0.2]
        ));
    }
}
