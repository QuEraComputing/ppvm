//! Typed AST for Stim with PPVM tag-based extensions promoted to
//! first-class instruction variants.

use crate::ast::{RawInstruction, Tag};

#[derive(Debug, Clone, PartialEq)]
pub struct ExtendedProgram {
    pub instructions: Vec<ExtendedInstruction>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ExtendedInstruction {
    /// Pass-through from vanilla Stim — covers `Gate`, `Noise`, `Measure`,
    /// `Annotation`. `MPad` and `Repeat` are NOT in `Raw` because their
    /// bits/body shapes diverge between dialects.
    Raw(RawInstruction),

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
    for instr in instructions {
        match instr {
            ExtendedInstruction::Raw(RawInstruction::Measure { targets, .. }) => {
                total = total.saturating_add(targets.len().saturating_mul(factor as usize));
            }
            ExtendedInstruction::MPad { bits, .. } => {
                total = total.saturating_add(bits.len().saturating_mul(factor as usize));
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
