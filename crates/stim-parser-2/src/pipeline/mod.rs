// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Typestate lowering pipeline: Source → Parsed → Validated → Lowered.
//! Each transition consumes `self` and returns the next state type, so
//! illegal orderings do not compile.

mod lower;
mod validate;

use std::sync::Arc;

use chumsky::Parser;

use crate::ast::{ExtendedProgram, Program};
use crate::diagnostics::{Aborted, Diagnostic, DiagnosticSink, Flow, LineMap, Span};
use crate::syntax::{RawSyntaxTree, program_parser};

pub struct Pipeline<S> {
    state: S,
}

pub struct Source<'a> {
    src: &'a str,
}
#[expect(
    dead_code,
    reason = "fields consumed by validate transition in Task 12"
)]
pub struct Parsed {
    pub(crate) tree: RawSyntaxTree,
    pub(crate) line_map: Arc<LineMap>,
}
// Once a Program exists it owns the LineMap, so the later states need only
// hold the program.
pub struct Validated {
    pub(crate) program: Program,
}
pub struct Lowered {
    pub(crate) program: ExtendedProgram,
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

impl Pipeline<Validated> {
    pub fn finish(self) -> Program {
        self.state.program
    }
}

impl Pipeline<Lowered> {
    pub fn finish(self) -> ExtendedProgram {
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
}
