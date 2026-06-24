// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Typed AST for Stim with PPVM tag-based extensions promoted to
//! first-class instruction variants.

use std::sync::Arc;

use crate::ast::shared::{AnnotationOp, Axis, GateOp, MeasureOp, MppOp, NoiseOp, Tag};
use crate::diagnostics::{LineMap, Span};

#[derive(Debug, Clone, PartialEq)]
pub enum ExtendedInstruction {
    // Pass-through families — the SAME structs as the vanilla AST.
    Gate(GateOp),
    Noise(NoiseOp),
    Measure(MeasureOp),
    Annotation(AnnotationOp),
    Mpp(MppOp),

    // Promoted sugar.
    T {
        targets: Vec<usize>,
        span: Span,
    },
    TDag {
        targets: Vec<usize>,
        span: Span,
    },
    Rotation {
        axis: Axis,
        theta: f64,
        targets: Vec<usize>,
        span: Span,
    },
    U3 {
        theta: f64,
        phi: f64,
        lambda: f64,
        targets: Vec<usize>,
        span: Span,
    },
    Loss {
        p: f64,
        targets: Vec<usize>,
        span: Span,
    },
    CorrelatedLoss {
        ps: [f64; 3],
        targets: Vec<(usize, usize)>,
        span: Span,
    },
    MPad {
        tags: Vec<Tag>,
        prob: Option<f64>,
        bits: Vec<bool>,
        span: Span,
    },
    Repeat {
        count: u64,
        body: Vec<ExtendedInstruction>,
        span: Span,
    },
}

#[derive(Debug, Clone)]
pub struct ExtendedProgram {
    pub instructions: Vec<ExtendedInstruction>,
    pub line_map: Arc<LineMap>,
}

impl PartialEq for ExtendedProgram {
    fn eq(&self, other: &Self) -> bool {
        self.instructions == other.instructions
    }
}

impl ExtendedProgram {
    /// Total recorded bits the program produces, accounting for REPEAT
    /// factors. Pure AST property; backend-agnostic.
    pub fn measurement_count(&self) -> usize {
        count_in_slice(&self.instructions, 1)
    }
}

fn count_in_slice(instructions: &[ExtendedInstruction], factor: u64) -> usize {
    let mut total = 0usize;
    let factor_usize = usize::try_from(factor).unwrap_or(usize::MAX);
    for instr in instructions {
        match instr {
            ExtendedInstruction::Measure(op) => {
                total = total.saturating_add(op.targets.len().saturating_mul(factor_usize));
            }
            ExtendedInstruction::MPad { bits, .. } => {
                total = total.saturating_add(bits.len().saturating_mul(factor_usize));
            }
            ExtendedInstruction::Mpp(op) => {
                total = total.saturating_add(op.products.len().saturating_mul(factor_usize));
            }
            ExtendedInstruction::Repeat { count, body, .. } => {
                total = total.saturating_add(count_in_slice(body, factor.saturating_mul(*count)));
            }
            ExtendedInstruction::Gate(_)
            | ExtendedInstruction::Noise(_)
            | ExtendedInstruction::Annotation(_)
            | ExtendedInstruction::T { .. }
            | ExtendedInstruction::TDag { .. }
            | ExtendedInstruction::Rotation { .. }
            | ExtendedInstruction::U3 { .. }
            | ExtendedInstruction::Loss { .. }
            | ExtendedInstruction::CorrelatedLoss { .. } => {}
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ast::shared::{GateOp, MeasureOp};
    use crate::diagnostics::{LineMap, Span};
    use crate::instructions::{GateName, MeasureName};
    use std::sync::Arc;

    fn span() -> Span {
        Span::new(0, 1)
    }

    #[test]
    fn measurement_count_scales_with_repeat() {
        let m = ExtendedInstruction::Measure(MeasureOp {
            name: MeasureName::M,
            tags: vec![],
            args: vec![],
            targets: vec![0, 1],
            span: span(),
        });
        let prog = ExtendedProgram {
            instructions: vec![ExtendedInstruction::Repeat {
                count: 3,
                body: vec![m],
                span: span(),
            }],
            line_map: Arc::new(LineMap::new("")),
        };
        assert_eq!(prog.measurement_count(), 6);
    }

    #[test]
    fn gate_op_is_shared_with_vanilla() {
        let _ = ExtendedInstruction::Gate(GateOp {
            name: GateName::H,
            tags: vec![],
            args: vec![],
            targets: vec![],
            span: span(),
        });
    }
}
