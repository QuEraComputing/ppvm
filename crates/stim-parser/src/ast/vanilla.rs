// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Vanilla Stim AST. Tags are preserved verbatim; the parser does not
//! resolve the Stim dialect — that is the consumer's responsibility.

use std::sync::Arc;

use crate::ast::shared::{AnnotationOp, GateOp, MeasureOp, MppOp, NoiseOp, Tag};
use crate::diagnostics::{LineMap, Span};

#[derive(Debug, Clone, PartialEq)]
pub enum Instruction {
    Gate(GateOp),
    Noise(NoiseOp),
    Measure(MeasureOp),
    Annotation(AnnotationOp),
    Mpp(MppOp),
    MPad {
        tags: Vec<Tag>,
        prob: Option<f64>,
        bits: Vec<usize>,
        span: Span,
    },
    Repeat {
        count: u64,
        body: Vec<Instruction>,
        span: Span,
    },
}

/// Vanilla program. Owns the `LineMap` so consumers can resolve any node's
/// `span` to `line:col`. Equality is by `instructions` only — the line map
/// is positional metadata, not identity.
#[derive(Debug, Clone)]
pub struct Program {
    pub instructions: Vec<Instruction>,
    pub line_map: Arc<LineMap>,
}

impl PartialEq for Program {
    fn eq(&self, other: &Self) -> bool {
        self.instructions == other.instructions
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ast::shared::GateOp;
    use crate::diagnostics::{LineMap, Span};
    use crate::instructions::GateName;
    use std::sync::Arc;

    #[test]
    fn program_holds_instructions() {
        let p = Program {
            instructions: vec![Instruction::Gate(GateOp {
                name: GateName::H,
                tags: vec![],
                args: vec![],
                targets: vec![],
                span: Span::new(0, 1),
            })],
            line_map: Arc::new(LineMap::new("H 0")),
        };
        assert_eq!(p.instructions.len(), 1);
    }

    #[test]
    fn program_eq_ignores_line_map() {
        let g = || {
            Instruction::Gate(GateOp {
                name: GateName::H,
                tags: vec![],
                args: vec![],
                targets: vec![],
                span: Span::new(0, 1),
            })
        };
        let a = Program {
            instructions: vec![g()],
            line_map: Arc::new(LineMap::new("H 0")),
        };
        let b = Program {
            instructions: vec![g()],
            line_map: Arc::new(LineMap::new("\n\nH 0")),
        };
        assert_eq!(a, b); // equality is by instructions only
    }
}
