// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

pub mod ast;
pub mod diagnostics;
pub mod instructions;
pub mod pipeline;
pub mod print;
pub(crate) mod syntax;

use std::sync::Arc;

use crate::ast::{ExtendedProgram, Program};
use crate::diagnostics::{Diagnostics, FailFast, LineMap};
use crate::pipeline::Pipeline;
use crate::syntax::run_on_parser_stack;

/// Build a [`Diagnostics`] from a fail-fast sink that aborted a parse stage.
fn fail(sink: FailFast, src: &str) -> Diagnostics {
    Diagnostics::new(sink.into_items(), Arc::new(LineMap::new(src)))
}

/// Parse Stim source into the vanilla [`Program`] AST. Uses a fail-fast
/// policy; the returned [`Diagnostics`] holds the first error.
pub fn parse(src: &str) -> Result<Program, Diagnostics> {
    run_on_parser_stack(|| {
        let mut sink = FailFast::new();
        let result = Pipeline::new(src)
            .parse(&mut sink)
            .and_then(|p| p.validate(&mut sink));
        match result {
            Ok(p) => Ok(p.finish()),
            Err(_) => Err(fail(sink, src)),
        }
    })
}

/// Parse Stim source into the extended-dialect [`ExtendedProgram`] AST.
pub fn parse_extended(src: &str) -> Result<ExtendedProgram, Diagnostics> {
    run_on_parser_stack(|| {
        let mut sink = FailFast::new();
        let result = Pipeline::new(src)
            .parse(&mut sink)
            .and_then(|p| p.validate(&mut sink))
            .and_then(|p| p.lower(&mut sink));
        match result {
            Ok(p) => Ok(p.finish()),
            Err(_) => Err(fail(sink, src)),
        }
    })
}

pub mod prelude {
    pub use crate::ast::{
        AnnotationOp, Axis, ExtendedInstruction, ExtendedProgram, GateOp, Instruction, MeasureOp,
        MppOp, NoiseOp, PauliAxis, PauliFactor, Program, Tag, TagParam, Target,
    };
    pub use crate::diagnostics::{
        Diagnostic, DiagnosticSink, Diagnostics, Flow, LineMap, Severity, Span,
    };
    pub use crate::instructions::{AnnotationKind, GateName, MeasureName, NoiseName};
    pub use crate::pipeline::Pipeline;
    pub use crate::print::{PrintOptions, StimPrint};
    pub use crate::{parse, parse_extended};
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_returns_program() {
        let prog = parse("H 0\nM 0\n").unwrap();
        assert_eq!(prog.instructions.len(), 2);
    }

    #[test]
    fn parse_extended_returns_extended_program() {
        let prog = parse_extended("S[T] 0\n").unwrap();
        assert_eq!(prog.instructions.len(), 1);
    }

    #[test]
    fn parse_error_renders_line_col() {
        let err = parse("REPEAT 2 {\nH 0\n").unwrap_err();
        assert!(err.to_string().starts_with("error at line"));
    }
}
