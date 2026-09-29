// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Typed AST for Stim with PPVM tag-based extensions promoted to
//! first-class instruction variants.

use std::sync::Arc;

use crate::ast::shared::{AnnotationOp, Axis, GateOp, MeasureOp, MppOp, NoiseOp};
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
    Leakage {
        p0: f64,
        p1: f64,
        targets: Vec<usize>,
        span: Span,
    },
    MPad {
        tag: String,
        prob: Option<f64>,
        bits: Vec<bool>,
        span: Span,
    },
    Repeat {
        tag: String,
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

    /// Number of qubits the program operates on: one past the highest qubit
    /// index referenced by any executable instruction, or `0` if it touches no
    /// qubits. Annotations (`DETECTOR`, `QUBIT_COORDS`, …) are ignored — their
    /// operands are measurement-record lookbacks or coordinates, not executable
    /// qubits. Pure AST property; backend-agnostic, mirrors [`measurement_count`].
    pub fn num_qubits(&self) -> usize {
        max_qubit_in_slice(&self.instructions).map_or(0, |m| m + 1)
    }
}

/// Highest qubit index referenced by any executable instruction in `slice`,
/// recursing into `REPEAT` bodies. `None` if nothing touches a qubit.
fn max_qubit_in_slice(instructions: &[ExtendedInstruction]) -> Option<usize> {
    let mut max: Option<usize> = None;
    for instr in instructions {
        // `Option<usize>` orders `None` below every `Some`, so `max.max(local)`
        // tracks the running maximum and treats "no qubit" as absent.
        let local = match instr {
            ExtendedInstruction::Gate(op) => op.targets.iter().filter_map(|t| t.as_qubit()).max(),
            ExtendedInstruction::Noise(op) => op.targets.iter().copied().max(),
            ExtendedInstruction::Measure(op) => op.targets.iter().copied().max(),
            ExtendedInstruction::Mpp(op) => op.products.iter().flatten().map(|f| f.qubit).max(),
            ExtendedInstruction::T { targets, .. }
            | ExtendedInstruction::TDag { targets, .. }
            | ExtendedInstruction::Rotation { targets, .. }
            | ExtendedInstruction::U3 { targets, .. }
            | ExtendedInstruction::Loss { targets, .. }
            | ExtendedInstruction::Leakage { targets, .. } => targets.iter().copied().max(),
            ExtendedInstruction::CorrelatedLoss { targets, .. } => {
                targets.iter().flat_map(|&(a, b)| [a, b]).max()
            }
            ExtendedInstruction::Repeat { body, .. } => max_qubit_in_slice(body),
            ExtendedInstruction::Annotation(_) | ExtendedInstruction::MPad { .. } => None,
        };
        max = max.max(local);
    }
    max
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
            | ExtendedInstruction::CorrelatedLoss { .. }
            | ExtendedInstruction::Leakage { .. } => {}
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ast::shared::{GateOp, MeasureOp, Target};
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
            tag: String::new(),
            args: vec![],
            targets: vec![0, 1],
            span: span(),
        });
        let prog = ExtendedProgram {
            instructions: vec![ExtendedInstruction::Repeat {
                tag: String::new(),
                count: 3,
                body: vec![m],
                span: span(),
            }],
            line_map: Arc::new(LineMap::new("")),
        };
        assert_eq!(prog.measurement_count(), 6);
    }

    #[test]
    fn num_qubits_is_one_past_highest_index() {
        // Gate on qubits {0, 4}, measure {2} -> 5 qubits (indices 0..=4).
        let prog = ExtendedProgram {
            instructions: vec![
                ExtendedInstruction::Gate(GateOp {
                    name: GateName::H,
                    tag: String::new(),
                    args: vec![],
                    targets: vec![Target::Qubit(0), Target::Qubit(4)],
                    span: span(),
                }),
                ExtendedInstruction::Measure(MeasureOp {
                    name: MeasureName::M,
                    tag: String::new(),
                    args: vec![],
                    targets: vec![2],
                    span: span(),
                }),
            ],
            line_map: Arc::new(LineMap::new("")),
        };
        assert_eq!(prog.num_qubits(), 5);
    }

    #[test]
    fn num_qubits_recurses_into_repeat() {
        let prog = ExtendedProgram {
            instructions: vec![ExtendedInstruction::Repeat {
                tag: String::new(),
                count: 3,
                body: vec![ExtendedInstruction::Measure(MeasureOp {
                    name: MeasureName::M,
                    tag: String::new(),
                    args: vec![],
                    targets: vec![7],
                    span: span(),
                })],
                span: span(),
            }],
            line_map: Arc::new(LineMap::new("")),
        };
        assert_eq!(prog.num_qubits(), 8);
    }

    #[test]
    fn num_qubits_is_zero_for_no_qubit_program() {
        // Empty program, and an annotation-only program, both touch no qubits.
        let empty = ExtendedProgram {
            instructions: vec![],
            line_map: Arc::new(LineMap::new("")),
        };
        assert_eq!(empty.num_qubits(), 0);
    }

    #[test]
    fn gate_op_is_shared_with_vanilla() {
        let _ = ExtendedInstruction::Gate(GateOp {
            name: GateName::H,
            tag: String::new(),
            args: vec![],
            targets: vec![],
            span: span(),
        });
    }
}
