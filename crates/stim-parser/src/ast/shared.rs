// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Shared leaf types and per-family payload structs reused by both the vanilla
//! `Instruction` AST and the `ExtendedInstruction` AST.

use crate::diagnostics::Span;
use crate::instructions::{AnnotationKind, GateName, MeasureName, NoiseName};

/// A gate operand. Most operands are plain qubit indices, but the control
/// of a classically-controlled gate — e.g. the `rec[-1]` in `CX rec[-1] 1` —
/// is a measurement-record lookback rather than a qubit. Mirrors the qubit /
/// record distinction Stim draws in its `GateTarget`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Target {
    /// A qubit index.
    Qubit(usize),
    /// A measurement-record lookback `rec[-k]`; stores `k >= 1` (the number
    /// of measurements to count back from the most recent).
    Rec(usize),
}

impl Target {
    /// The qubit index, if this is a plain qubit target.
    pub fn as_qubit(self) -> Option<usize> {
        match self {
            Target::Qubit(q) => Some(q),
            Target::Rec(_) => None,
        }
    }
}

/// A qubit target compares equal to its bare index; a record target never
/// does. Lets callers (and tests) match qubit targets against plain `usize`s.
impl PartialEq<usize> for Target {
    fn eq(&self, other: &usize) -> bool {
        matches!(self, Target::Qubit(q) if q == other)
    }
}

/// The Pauli basis of a single-qubit factor in an `MPP` product.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PauliAxis {
    X,
    Y,
    Z,
}

impl PauliAxis {
    /// The single-character Stim spelling (`X`/`Y`/`Z`).
    pub fn as_char(self) -> char {
        match self {
            PauliAxis::X => 'X',
            PauliAxis::Y => 'Y',
            PauliAxis::Z => 'Z',
        }
    }
}

/// One single-qubit Pauli factor of an `MPP` product, e.g. the `Y3` in
/// `MPP X0*Y3*Z7`: the measured `axis` on its support `qubit`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PauliFactor {
    pub axis: PauliAxis,
    pub qubit: usize,
}

/// The rotation axis for an extended-dialect `R_X` / `R_Y` / `R_Z` rotation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    X,
    Y,
    Z,
}

// ---------------------------------------------------------------------------
// Shared per-family payload structs
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
pub struct GateOp {
    pub name: GateName,
    pub tag: String,
    pub args: Vec<f64>,
    pub targets: Vec<Target>,
    pub span: Span,
}

#[derive(Debug, Clone, PartialEq)]
pub struct NoiseOp {
    pub name: NoiseName,
    pub tag: String,
    pub args: Vec<f64>,
    pub targets: Vec<usize>,
    pub span: Span,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MeasureOp {
    pub name: MeasureName,
    pub tag: String,
    pub args: Vec<f64>,
    pub targets: Vec<usize>,
    pub span: Span,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AnnotationOp {
    pub kind: AnnotationKind,
    pub tag: String,
    pub args: Vec<f64>,
    pub targets: Vec<usize>,
    pub span: Span,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MppOp {
    pub tag: String,
    pub args: Vec<f64>,
    pub products: Vec<Vec<PauliFactor>>,
    pub span: Span,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn target_equals_bare_qubit_index() {
        assert_eq!(Target::Qubit(3), 3usize);
        assert_ne!(Target::Rec(1), 0usize);
        assert_eq!(Target::Qubit(3).as_qubit(), Some(3));
        assert_eq!(Target::Rec(1).as_qubit(), None);
    }

    #[test]
    fn pauli_axis_char() {
        assert_eq!(PauliAxis::X.as_char(), 'X');
        assert_eq!(PauliAxis::Z.as_char(), 'Z');
    }
}
