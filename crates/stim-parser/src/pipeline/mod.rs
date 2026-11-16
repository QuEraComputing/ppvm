// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Typestate lowering pipeline:
//!
//! ```text
//! Source ─parse─▶ Parsed ─validate─▶ Validated ─fuse─▶ Fused<Program>
//!                                        │
//!                                        └─lower─▶ Lowered ─fuse─▶ Fused<ExtendedProgram>
//! ```
//!
//! Each transition consumes `self` and returns the next state type, so
//! illegal orderings do not compile.

mod fuse;
mod lower;
mod validate;

use std::sync::Arc;

use chumsky::Parser;

use crate::ast::{ExtendedProgram, Program};
use crate::diagnostics::{Aborted, Diagnostic, DiagnosticSink, Flow, LineMap, Severity, Span};
use crate::syntax::{RawSyntaxTree, program_parser};

/// Emit a diagnostic to the sink and report the sink's continuation decision.
pub(crate) fn emit(
    sink: &mut dyn DiagnosticSink,
    span: Span,
    code: &'static str,
    message: String,
) -> Flow {
    sink.emit(Diagnostic::error(span, code, message))
}

/// Emit a recoverable error, then either abort the stage or signal "skip this item".
pub(crate) fn emit_skip<T>(
    sink: &mut dyn DiagnosticSink,
    span: Span,
    code: &'static str,
    message: String,
) -> Result<Option<T>, Aborted> {
    if emit(sink, span, code, message) == Flow::Abort {
        Err(Aborted)
    } else {
        Ok(None)
    }
}

/// Forwards to the wrapped sink, remembering whether an error went through.
struct TrackErrors<'s> {
    sink: &'s mut dyn DiagnosticSink,
    seen: bool,
}

impl<'s> TrackErrors<'s> {
    fn new(sink: &'s mut dyn DiagnosticSink) -> Self {
        TrackErrors { sink, seen: false }
    }
}

impl DiagnosticSink for TrackErrors<'_> {
    fn emit(&mut self, diagnostic: Diagnostic) -> Flow {
        self.seen |= diagnostic.severity == Severity::Error;
        self.sink.emit(diagnostic)
    }
}

pub struct Pipeline<S> {
    state: S,
}

pub struct Source<'a> {
    src: &'a str,
}
pub struct Parsed {
    pub(crate) tree: RawSyntaxTree,
    pub(crate) line_map: Arc<LineMap>,
}
// Once a Program exists it owns the LineMap, so the later states need only
// hold the program, and whether any stage so far reported an error.
pub struct Validated {
    pub(crate) program: Program,
    pub(crate) had_errors: bool,
}
pub struct Lowered {
    pub(crate) program: ExtendedProgram,
    pub(crate) had_errors: bool,
}
pub struct Fused<P> {
    pub(crate) program: P,
}

impl<'a> Pipeline<Source<'a>> {
    pub fn new(src: &'a str) -> Self {
        Pipeline {
            state: Source { src },
        }
    }

    /// Stage 1: pure syntax. Forwards every chumsky error into the sink.
    pub fn parse(self, sink: &mut dyn DiagnosticSink) -> Result<Pipeline<Parsed>, Aborted> {
        let src = self.state.src;
        let line_map = Arc::new(LineMap::new(src));
        let result = program_parser().parse(src);
        match result.into_result() {
            Ok(tree) => Ok(Pipeline {
                state: Parsed { tree, line_map },
            }),
            Err(errors) => {
                for err in errors {
                    let span: Span = (*err.span()).into();
                    let flow = sink.emit(Diagnostic::error(span, "syntax", err.to_string()));
                    if flow == Flow::Abort {
                        return Err(Aborted);
                    }
                }
                // All syntax errors forwarded; with a non-aborting sink we still
                // cannot produce a tree, so abort the stage.
                Err(Aborted)
            }
        }
    }
}

impl Pipeline<Parsed> {
    /// Stage 2: table-driven validation. Builds a typed vanilla [`Program`],
    /// forwarding every recoverable error to the sink.
    pub fn validate(self, sink: &mut dyn DiagnosticSink) -> Result<Pipeline<Validated>, Aborted> {
        let Parsed { tree, line_map } = self.state;
        let mut sink = TrackErrors::new(sink);
        let program = validate::validate(tree, &line_map, &mut sink)?;
        Ok(Pipeline {
            state: Validated {
                program,
                had_errors: sink.seen,
            },
        })
    }
}

impl Pipeline<Validated> {
    /// Stage 3: lowering. Promotes tag-based PPVM extensions to first-class
    /// [`ExtendedInstruction`] variants, forwarding every recoverable error to
    /// the sink.
    pub fn lower(self, sink: &mut dyn DiagnosticSink) -> Result<Pipeline<Lowered>, Aborted> {
        let Validated {
            program,
            had_errors,
        } = self.state;
        let mut sink = TrackErrors::new(sink);
        let program = lower::lower(program, &mut sink)?;
        Ok(Pipeline {
            state: Lowered {
                program,
                had_errors: had_errors || sink.seen,
            },
        })
    }

    /// Stage 4: fusion. Merges consecutive compatible instructions, as Stim
    /// does (`H 0` / `H 1` becomes `H 0 1`). Skipped if any stage reported an
    /// error, since the program is then incomplete.
    pub fn fuse(self) -> Pipeline<Fused<Program>> {
        let Validated {
            mut program,
            had_errors,
        } = self.state;
        if !had_errors {
            program.instructions = fuse::fuse(program.instructions);
        }
        Pipeline {
            state: Fused { program },
        }
    }

    pub fn finish(self) -> Program {
        self.state.program
    }
}

impl Pipeline<Lowered> {
    /// Stage 4: fusion. Merges consecutive compatible instructions, as Stim
    /// does (`H 0` / `H 1` becomes `H 0 1`). Skipped if any stage reported an
    /// error, since the program is then incomplete.
    pub fn fuse(self) -> Pipeline<Fused<ExtendedProgram>> {
        let Lowered {
            mut program,
            had_errors,
        } = self.state;
        if !had_errors {
            program.instructions = fuse::fuse(program.instructions);
        }
        Pipeline {
            state: Fused { program },
        }
    }

    pub fn finish(self) -> ExtendedProgram {
        self.state.program
    }
}

impl<P> Pipeline<Fused<P>> {
    pub fn finish(self) -> P {
        self.state.program
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagnostics::FailFast;

    #[test]
    fn parse_transition_produces_parsed_state() {
        let mut sink = FailFast::new();
        let parsed = Pipeline::new("H 0\n").parse(&mut sink);
        assert!(parsed.is_ok());
    }

    #[test]
    fn parse_transition_emits_diagnostic_and_aborts_on_syntax_error() {
        let mut sink = FailFast::new();
        let res = Pipeline::new("REPEAT 2 {\nH 0\n").parse(&mut sink); // unclosed
        assert!(res.is_err());
        assert!(sink.saw_error());
    }

    #[test]
    fn parse_then_validate_produces_program() {
        use crate::ast::Instruction;
        use crate::instructions::GateName;

        let mut sink = FailFast::new();
        let program = Pipeline::new("H 0\nCX 0 1\n")
            .parse(&mut sink)
            .expect("parse")
            .validate(&mut sink)
            .expect("validate")
            .finish();
        assert_eq!(program.instructions.len(), 2);
        assert!(matches!(
            &program.instructions[0],
            Instruction::Gate(op) if op.name == GateName::H
        ));
    }

    #[test]
    fn validate_transition_emits_and_aborts_on_failfast() {
        let mut sink = FailFast::new();
        let parsed = Pipeline::new("FROBNICATE 0\n")
            .parse(&mut sink)
            .expect("parse");
        let res = parsed.validate(&mut sink);
        assert!(res.is_err());
        assert!(sink.saw_error());
    }
}
