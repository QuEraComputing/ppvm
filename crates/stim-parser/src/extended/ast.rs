// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

//! Typed AST for Stim with PPVM tag-based extensions promoted to
//! first-class instruction variants.

use crate::ast::{AnnotationKind, GateName, MeasureName, NoiseName, RawInstruction, Tag};

#[derive(Debug, Clone, PartialEq)]
pub struct ExtendedProgram {
    pub instructions: Vec<ExtendedInstruction>,
}

/// Subset of [`RawInstruction`] that passes through the extended-dialect
/// interpreter unchanged. Excludes `MPad` and `Repeat`, which are always
/// lowered to the typed [`ExtendedInstruction::MPad`] /
/// [`ExtendedInstruction::Repeat`] variants — so by construction
/// `ExtendedInstruction::Raw(_)` can never wrap one of them.
#[derive(Debug, Clone, PartialEq)]
pub enum RawPassthrough {
    Gate {
        name: GateName,
        tags: Vec<Tag>,
        args: Vec<f64>,
        targets: Vec<usize>,
        line: usize,
    },
    Noise {
        name: NoiseName,
        tags: Vec<Tag>,
        args: Vec<f64>,
        targets: Vec<usize>,
        line: usize,
    },
    Measure {
        name: MeasureName,
        tags: Vec<Tag>,
        args: Vec<f64>,
        targets: Vec<usize>,
        line: usize,
    },
    Annotation {
        kind: AnnotationKind,
        args: Vec<f64>,
        targets: Vec<usize>,
        line: usize,
    },
}

impl RawPassthrough {
    /// Lift back into [`RawInstruction`]. Used by the printer to share
    /// formatting code with the vanilla AST.
    pub fn into_raw(self) -> RawInstruction {
        match self {
            RawPassthrough::Gate {
                name,
                tags,
                args,
                targets,
                line,
            } => RawInstruction::Gate {
                name,
                tags,
                args,
                targets,
                line,
            },
            RawPassthrough::Noise {
                name,
                tags,
                args,
                targets,
                line,
            } => RawInstruction::Noise {
                name,
                tags,
                args,
                targets,
                line,
            },
            RawPassthrough::Measure {
                name,
                tags,
                args,
                targets,
                line,
            } => RawInstruction::Measure {
                name,
                tags,
                args,
                targets,
                line,
            },
            RawPassthrough::Annotation {
                kind,
                args,
                targets,
                line,
            } => RawInstruction::Annotation {
                kind,
                args,
                targets,
                line,
            },
        }
    }

    pub fn line(&self) -> usize {
        match self {
            RawPassthrough::Gate { line, .. }
            | RawPassthrough::Noise { line, .. }
            | RawPassthrough::Measure { line, .. }
            | RawPassthrough::Annotation { line, .. } => *line,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ExtendedInstruction {
    /// Pass-through from vanilla Stim — covers `Gate`, `Noise`, `Measure`,
    /// `Annotation`. `MPad` and `Repeat` get their own typed variants on
    /// [`ExtendedInstruction`] because their bits/body shapes diverge
    /// between dialects.
    Raw(RawPassthrough),

    // --- Extended-dialect sugar variants ---
    T {
        targets: Vec<usize>,
        line: usize,
    },
    TDag {
        targets: Vec<usize>,
        line: usize,
    },
    Rotation {
        axis: Axis,
        theta: f64,
        targets: Vec<usize>,
        line: usize,
    },
    U3 {
        theta: f64,
        phi: f64,
        lambda: f64,
        targets: Vec<usize>,
        line: usize,
    },
    Loss {
        p: f64,
        targets: Vec<usize>,
        line: usize,
    },
    CorrelatedLoss {
        ps: [f64; 3],
        targets: Vec<(usize, usize)>,
        line: usize,
    },
    /// Extended-dialect MPad with validated bits (each ∈ {0, 1}).
    MPad {
        tags: Vec<Tag>,
        prob: Option<f64>,
        bits: Vec<bool>,
        line: usize,
    },
    Repeat {
        count: u64,
        body: Vec<ExtendedInstruction>,
        line: usize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    X,
    Y,
    Z,
}

impl ExtendedProgram {
    /// Total number of recorded bits the program will produce, accounting for
    /// `REPEAT` factors. Pure AST property; backend-agnostic.
    pub fn measurement_count(&self) -> usize {
        count_in_slice(&self.instructions, 1)
    }
}

fn count_in_slice(instructions: &[ExtendedInstruction], factor: u64) -> usize {
    let mut total = 0usize;
    let factor_usize = usize::try_from(factor).unwrap_or(usize::MAX);
    for instr in instructions {
        match instr {
            ExtendedInstruction::Raw(RawPassthrough::Measure { targets, .. }) => {
                total = total.saturating_add(targets.len().saturating_mul(factor_usize));
            }
            ExtendedInstruction::MPad { bits, .. } => {
                total = total.saturating_add(bits.len().saturating_mul(factor_usize));
            }
            ExtendedInstruction::Repeat { count, body, .. } => {
                total = total.saturating_add(count_in_slice(body, factor.saturating_mul(*count)));
            }
            ExtendedInstruction::Raw(_)
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
