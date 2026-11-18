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
use crate::instructions::NoiseName;

/// Fuse a list of instructions, recursing into `REPEAT` bodies.
pub(crate) fn fuse<I: Fusable>(instructions: Vec<I>) -> Vec<I> {
    let mut out: Vec<I> = Vec::with_capacity(instructions.len());
    for mut instr in instructions {
        if let Some(body) = instr.body_mut() {
            *body = fuse(mem::take(body));
        }
        let Some(prev) = out.last_mut() else {
            out.push(instr);
            continue;
        };
        if let Err(instr) = prev.try_fuse(instr) {
            out.push(instr);
        }
    }
    out
}

/// An instruction type the fusion pass can work on.
pub(crate) trait Fusable: Sized {
    /// The body of a `REPEAT` block, fused independently.
    fn body_mut(&mut self) -> Option<&mut Vec<Self>>;
    /// Append `next`'s targets onto `self`, or hand `next` back if the two
    /// can't be fused.
    fn try_fuse(&mut self, next: Self) -> Result<(), Self>;
}

// The shared op structs fuse on the same name, tag and args.

fn gates_fuse(a: &GateOp, b: &GateOp) -> bool {
    a.name == b.name && a.tag == b.tag && a.args == b.args
}

fn noises_fuse(a: &NoiseOp, b: &NoiseOp) -> bool {
    // A correlated error is one joint event over all its targets.
    let correlated = matches!(
        a.name,
        NoiseName::CorrelatedError | NoiseName::ElseCorrelatedError
    );
    !correlated && a.name == b.name && a.tag == b.tag && a.args == b.args
}

fn measures_fuse(a: &MeasureOp, b: &MeasureOp) -> bool {
    a.name == b.name && a.tag == b.tag && a.args == b.args
}

fn mpps_fuse(a: &MppOp, b: &MppOp) -> bool {
    a.tag == b.tag && a.args == b.args
}

/// Never fused: annotations (`TICK`, `DETECTOR`, `OBSERVABLE_INCLUDE`,
/// coordinates), `REPEAT` blocks, and `CORRELATED_ERROR` /
/// `ELSE_CORRELATED_ERROR`.
impl Fusable for Instruction {
    fn body_mut(&mut self) -> Option<&mut Vec<Self>> {
        match self {
            Instruction::Repeat { body, .. } => Some(body),
            _ => None,
        }
    }

    fn try_fuse(&mut self, next: Self) -> Result<(), Self> {
        use Instruction as I;
        match (self, next) {
            (I::Gate(a), I::Gate(b)) if gates_fuse(a, &b) => a.targets.extend(b.targets),
            (I::Noise(a), I::Noise(b)) if noises_fuse(a, &b) => a.targets.extend(b.targets),
            (I::Measure(a), I::Measure(b)) if measures_fuse(a, &b) => a.targets.extend(b.targets),
            (I::Mpp(a), I::Mpp(b)) if mpps_fuse(a, &b) => a.products.extend(b.products),
            (
                I::MPad {
                    tag, prob, bits, ..
                },
                I::MPad {
                    tag: t,
                    prob: p,
                    bits: b,
                    ..
                },
            ) if *tag == t && *prob == p => bits.extend(b),
            (_, next) => return Err(next),
        }
        Ok(())
    }
}

/// Fuses by operation, not spelling: `S[T] 0` / `T 1` becomes one `T`. Never
/// fused: annotations and `REPEAT` blocks.
impl Fusable for ExtendedInstruction {
    fn body_mut(&mut self) -> Option<&mut Vec<Self>> {
        match self {
            ExtendedInstruction::Repeat { body, .. } => Some(body),
            _ => None,
        }
    }

    fn try_fuse(&mut self, next: Self) -> Result<(), Self> {
        use ExtendedInstruction as E;
        match (self, next) {
            (E::Gate(a), E::Gate(b)) if gates_fuse(a, &b) => a.targets.extend(b.targets),
            (E::Noise(a), E::Noise(b)) if noises_fuse(a, &b) => a.targets.extend(b.targets),
            (E::Measure(a), E::Measure(b)) if measures_fuse(a, &b) => a.targets.extend(b.targets),
            (E::Mpp(a), E::Mpp(b)) if mpps_fuse(a, &b) => a.products.extend(b.products),
            (E::T { targets, .. }, E::T { targets: t, .. })
            | (E::TDag { targets, .. }, E::TDag { targets: t, .. })
            | (E::ResetLoss { targets, .. }, E::ResetLoss { targets: t, .. })
            | (E::ResetLeakage { targets, .. }, E::ResetLeakage { targets: t, .. })
            | (E::ResetState { targets, .. }, E::ResetState { targets: t, .. }) => {
                targets.extend(t)
            }
            (
                E::Rotation {
                    axis,
                    theta,
                    targets,
                    ..
                },
                E::Rotation {
                    axis: a,
                    theta: th,
                    targets: t,
                    ..
                },
            ) if *axis == a && *theta == th => targets.extend(t),
            (
                E::U3 {
                    theta,
                    phi,
                    lambda,
                    targets,
                    ..
                },
                E::U3 {
                    theta: th,
                    phi: ph,
                    lambda: l,
                    targets: t,
                    ..
                },
            ) if *theta == th && *phi == ph && *lambda == l => targets.extend(t),
            (
                E::Loss { p, targets, .. },
                E::Loss {
                    p: q, targets: t, ..
                },
            ) if *p == q => targets.extend(t),
            (
                E::Leakage {
                    p0, p1, targets, ..
                },
                E::Leakage {
                    p0: q0,
                    p1: q1,
                    targets: t,
                    ..
                },
            ) if *p0 == q0 && *p1 == q1 => targets.extend(t),
            (
                E::CorrelatedLoss { ps, targets, .. },
                E::CorrelatedLoss {
                    ps: qs, targets: t, ..
                },
            ) if *ps == qs => targets.extend(t),
            (
                E::MPad {
                    tag, prob, bits, ..
                },
                E::MPad {
                    tag: t,
                    prob: p,
                    bits: b,
                    ..
                },
            ) if *tag == t && *prob == p => bits.extend(b),
            (_, next) => return Err(next),
        }
        Ok(())
    }
}
