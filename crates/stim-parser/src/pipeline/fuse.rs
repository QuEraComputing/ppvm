// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Stage 4: fusion. Like Stim's `Circuit` parser, merges consecutive
//! instructions with the same name, tag and arguments (`H 0` / `H 1` becomes
//! `H 0 1`). Runs after every stage that can reject input, so those checks
//! see one source line at a time.

use std::mem;

use crate::ast::extended::ExtendedInstruction;
use crate::ast::shared::{GateOp, MeasureOp, MppOp, NoiseOp};
use crate::ast::vanilla::Instruction;
use crate::diagnostics::Span;

/// Fuse a list of instructions, recursing into `REPEAT` bodies. Never fuses
/// across a `skipped` node.
pub(crate) fn fuse<I: Fusable>(instructions: Vec<I>, mut skipped: Vec<Span>) -> Vec<I> {
    // validate's and lower's records are each in source order, but concatenated.
    skipped.sort_unstable_by_key(|s| s.start);
    fuse_slice(instructions, &skipped)
}

fn fuse_slice<I: Fusable>(instructions: Vec<I>, skipped: &[Span]) -> Vec<I> {
    let mut out: Vec<I> = Vec::with_capacity(instructions.len());
    for mut instr in instructions {
        if let Some(body) = instr.body_mut() {
            *body = fuse_slice(mem::take(body), skipped);
        }
        let next = *instr.span_mut();
        if let Some(prev) = out.last_mut()
            && prev.can_fuse(&instr)
            && !skipped_between(skipped, *prev.span_mut(), next)
        {
            prev.append_operands(instr);
            prev.span_mut().end = next.end;
        } else {
            out.push(instr);
        }
    }
    out
}

/// Was any node skipped between `prev` and `next`? `skipped` is sorted by start.
fn skipped_between(skipped: &[Span], prev: Span, next: Span) -> bool {
    let i = skipped.partition_point(|s| s.start <= prev.start);
    skipped.get(i).is_some_and(|s| s.start < next.start)
}

/// An instruction type the fusion pass can work on.
pub(crate) trait Fusable: Sized {
    fn span_mut(&mut self) -> &mut Span;
    /// The body of a `REPEAT` block, fused independently.
    fn body_mut(&mut self) -> Option<&mut Vec<Self>>;
    /// Can `next` be fused onto the end of `self`?
    fn can_fuse(&self, next: &Self) -> bool;
    /// Append `next`'s targets onto `self`. Precondition: `self.can_fuse(&next)`.
    fn append_operands(&mut self, next: Self);
}

// The shared op structs fuse on the same name , tag and args.
// TODO: spelling aliases such as `CX`/`CNOT`

fn gates_fuse(a: &GateOp, b: &GateOp) -> bool {
    a.name == b.name && a.tag == b.tag && a.args == b.args
}

fn noises_fuse(a: &NoiseOp, b: &NoiseOp) -> bool {
    a.name.is_fusable() && a.name == b.name && a.tag == b.tag && a.args == b.args
}

fn measures_fuse(a: &MeasureOp, b: &MeasureOp) -> bool {
    a.name == b.name && a.tag == b.tag && a.args == b.args
}

fn mpps_fuse(a: &MppOp, b: &MppOp) -> bool {
    a.tag == b.tag && a.args == b.args
}

/// Never fused: annotations (`TICK`, `DETECTOR`, `OBSERVABLE_INCLUDE`,
/// coordinates), `REPEAT` blocks, and `CORRELATED_ERROR` /
/// `ELSE_CORRELATED_ERROR` (see `NoiseName::is_fusable`).
impl Fusable for Instruction {
    fn span_mut(&mut self) -> &mut Span {
        match self {
            Instruction::Gate(op) => &mut op.span,
            Instruction::Noise(op) => &mut op.span,
            Instruction::Measure(op) => &mut op.span,
            Instruction::Annotation(op) => &mut op.span,
            Instruction::Mpp(op) => &mut op.span,
            Instruction::MPad { span, .. } | Instruction::Repeat { span, .. } => span,
        }
    }

    fn body_mut(&mut self) -> Option<&mut Vec<Self>> {
        match self {
            Instruction::Repeat { body, .. } => Some(body),
            _ => None,
        }
    }

    fn can_fuse(&self, next: &Self) -> bool {
        use Instruction as I;
        match (self, next) {
            (I::Gate(a), I::Gate(b)) => gates_fuse(a, b),
            (I::Noise(a), I::Noise(b)) => noises_fuse(a, b),
            (I::Measure(a), I::Measure(b)) => measures_fuse(a, b),
            (I::Mpp(a), I::Mpp(b)) => mpps_fuse(a, b),
            (
                I::MPad { tag, prob, .. },
                I::MPad {
                    tag: t, prob: p, ..
                },
            ) => tag == t && prob == p,
            _ => false,
        }
    }

    fn append_operands(&mut self, next: Self) {
        use Instruction as I;
        match (self, next) {
            (I::Gate(a), I::Gate(b)) => a.targets.extend(b.targets),
            (I::Noise(a), I::Noise(b)) => a.targets.extend(b.targets),
            (I::Measure(a), I::Measure(b)) => a.targets.extend(b.targets),
            (I::Mpp(a), I::Mpp(b)) => a.products.extend(b.products),
            (I::MPad { bits, .. }, I::MPad { bits: b, .. }) => bits.extend(b),
            _ => unreachable!("append_operands called on a pair that can_fuse rejects"),
        }
    }
}

/// Fuses by operation, not spelling: `S[T] 0` / `T 1` becomes one `T`. Never
/// fused: annotations and `REPEAT` blocks.
impl Fusable for ExtendedInstruction {
    fn span_mut(&mut self) -> &mut Span {
        use ExtendedInstruction as E;
        match self {
            E::Gate(op) => &mut op.span,
            E::Noise(op) => &mut op.span,
            E::Measure(op) => &mut op.span,
            E::Annotation(op) => &mut op.span,
            E::Mpp(op) => &mut op.span,
            E::T { span, .. }
            | E::TDag { span, .. }
            | E::Rotation { span, .. }
            | E::U3 { span, .. }
            | E::Loss { span, .. }
            | E::CorrelatedLoss { span, .. }
            | E::Leakage { span, .. }
            | E::ResetLoss { span, .. }
            | E::ResetLeakage { span, .. }
            | E::ResetState { span, .. }
            | E::MPad { span, .. }
            | E::Repeat { span, .. } => span,
        }
    }

    fn body_mut(&mut self) -> Option<&mut Vec<Self>> {
        match self {
            ExtendedInstruction::Repeat { body, .. } => Some(body),
            _ => None,
        }
    }

    fn can_fuse(&self, next: &Self) -> bool {
        use ExtendedInstruction as E;
        match (self, next) {
            (E::Gate(a), E::Gate(b)) => gates_fuse(a, b),
            (E::Noise(a), E::Noise(b)) => noises_fuse(a, b),
            (E::Measure(a), E::Measure(b)) => measures_fuse(a, b),
            (E::Mpp(a), E::Mpp(b)) => mpps_fuse(a, b),
            (E::T { .. }, E::T { .. })
            | (E::TDag { .. }, E::TDag { .. })
            | (E::ResetLoss { .. }, E::ResetLoss { .. })
            | (E::ResetLeakage { .. }, E::ResetLeakage { .. })
            | (E::ResetState { .. }, E::ResetState { .. }) => true,
            (
                E::Rotation { axis, theta, .. },
                E::Rotation {
                    axis: a, theta: t, ..
                },
            ) => axis == a && theta == t,
            (
                E::U3 {
                    theta, phi, lambda, ..
                },
                E::U3 {
                    theta: t,
                    phi: p,
                    lambda: l,
                    ..
                },
            ) => theta == t && phi == p && lambda == l,
            (E::Loss { p, .. }, E::Loss { p: q, .. }) => p == q,
            (E::Leakage { p0, p1, .. }, E::Leakage { p0: q0, p1: q1, .. }) => p0 == q0 && p1 == q1,
            (E::CorrelatedLoss { ps, .. }, E::CorrelatedLoss { ps: qs, .. }) => ps == qs,
            (
                E::MPad { tag, prob, .. },
                E::MPad {
                    tag: t, prob: p, ..
                },
            ) => tag == t && prob == p,
            _ => false,
        }
    }

    fn append_operands(&mut self, next: Self) {
        use ExtendedInstruction as E;
        match (self, next) {
            (E::Gate(a), E::Gate(b)) => a.targets.extend(b.targets),
            (E::Noise(a), E::Noise(b)) => a.targets.extend(b.targets),
            (E::Measure(a), E::Measure(b)) => a.targets.extend(b.targets),
            (E::Mpp(a), E::Mpp(b)) => a.products.extend(b.products),
            (E::T { targets, .. }, E::T { targets: t, .. })
            | (E::TDag { targets, .. }, E::TDag { targets: t, .. })
            | (E::Rotation { targets, .. }, E::Rotation { targets: t, .. })
            | (E::U3 { targets, .. }, E::U3 { targets: t, .. })
            | (E::Loss { targets, .. }, E::Loss { targets: t, .. })
            | (E::Leakage { targets, .. }, E::Leakage { targets: t, .. })
            | (E::ResetLoss { targets, .. }, E::ResetLoss { targets: t, .. })
            | (E::ResetLeakage { targets, .. }, E::ResetLeakage { targets: t, .. })
            | (E::ResetState { targets, .. }, E::ResetState { targets: t, .. }) => {
                targets.extend(t)
            }
            (E::CorrelatedLoss { targets, .. }, E::CorrelatedLoss { targets: t, .. }) => {
                targets.extend(t)
            }
            (E::MPad { bits, .. }, E::MPad { bits: b, .. }) => bits.extend(b),
            _ => unreachable!("append_operands called on a pair that can_fuse rejects"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn skipped_between_only_counts_spans_between_the_two() {
        let skipped = [Span::new(5, 9), Span::new(20, 24)];
        assert!(skipped_between(
            &skipped,
            Span::new(0, 1),
            Span::new(10, 11)
        ));
        assert!(!skipped_between(
            &skipped,
            Span::new(10, 11),
            Span::new(15, 16)
        ));
        assert!(skipped_between(
            &skipped,
            Span::new(10, 11),
            Span::new(30, 31)
        ));
        assert!(!skipped_between(
            &skipped,
            Span::new(25, 26),
            Span::new(30, 31)
        ));
    }
}
