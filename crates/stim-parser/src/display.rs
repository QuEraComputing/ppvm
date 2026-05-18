// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

//! `Display` impls for both Stim AST layers.
//!
//! Output is canonical Stim: 4-space indentation inside `REPEAT` blocks,
//! one instruction per line, `[tags](args) targets…` ordering matching
//! the grammar. The printer is deterministic — `parse → print → parse →
//! print` reaches a fixpoint after one round (see `tests/roundtrip.rs`).

use std::fmt;

use crate::ast::{Program, RawInstruction, Tag, TagParam};
use crate::extended::ast::{Axis, ExtendedInstruction, ExtendedProgram, RawPassthrough};

const INDENT: &str = "    ";

impl fmt::Display for Program {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_raw_slice(&self.instructions, f, 0)
    }
}

impl fmt::Display for ExtendedProgram {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_ext_slice(&self.instructions, f, 0)
    }
}

fn fmt_raw_slice(
    instrs: &[RawInstruction],
    f: &mut fmt::Formatter<'_>,
    depth: usize,
) -> fmt::Result {
    for i in instrs {
        fmt_raw(i, f, depth)?;
    }
    Ok(())
}

fn fmt_ext_slice(
    instrs: &[ExtendedInstruction],
    f: &mut fmt::Formatter<'_>,
    depth: usize,
) -> fmt::Result {
    for i in instrs {
        fmt_ext(i, f, depth)?;
    }
    Ok(())
}

fn write_indent(f: &mut fmt::Formatter<'_>, depth: usize) -> fmt::Result {
    for _ in 0..depth {
        f.write_str(INDENT)?;
    }
    Ok(())
}

fn fmt_raw(i: &RawInstruction, f: &mut fmt::Formatter<'_>, depth: usize) -> fmt::Result {
    write_indent(f, depth)?;
    match i {
        RawInstruction::Gate {
            name,
            tags,
            args,
            targets,
            ..
        } => {
            f.write_str(name.canonical_name())?;
            write_tags(f, tags)?;
            write_args(f, args)?;
            write_usize_targets(f, targets)?;
        }
        RawInstruction::Noise {
            name,
            tags,
            args,
            targets,
            ..
        } => {
            f.write_str(name.canonical_name())?;
            write_tags(f, tags)?;
            write_args(f, args)?;
            write_usize_targets(f, targets)?;
        }
        RawInstruction::Measure {
            name,
            tags,
            args,
            targets,
            ..
        } => {
            f.write_str(name.canonical_name())?;
            write_tags(f, tags)?;
            write_args(f, args)?;
            write_usize_targets(f, targets)?;
        }
        RawInstruction::Annotation {
            kind,
            args,
            targets,
            ..
        } => {
            f.write_str(kind.canonical_name())?;
            write_args(f, args)?;
            write_usize_targets(f, targets)?;
        }
        RawInstruction::MPad {
            tags, prob, bits, ..
        } => {
            f.write_str("MPAD")?;
            write_tags(f, tags)?;
            if let Some(p) = prob {
                write!(f, "({})", FloatLit(*p))?;
            }
            write_usize_targets(f, bits)?;
        }
        RawInstruction::Repeat { count, body, .. } => {
            writeln!(f, "REPEAT {count} {{")?;
            fmt_raw_slice(body, f, depth + 1)?;
            write_indent(f, depth)?;
            f.write_str("}")?;
        }
    }
    writeln!(f)
}

fn fmt_raw_passthrough(
    i: &RawPassthrough,
    f: &mut fmt::Formatter<'_>,
    depth: usize,
) -> fmt::Result {
    write_indent(f, depth)?;
    match i {
        RawPassthrough::Gate {
            name,
            tags,
            args,
            targets,
            ..
        } => {
            f.write_str(name.canonical_name())?;
            write_tags(f, tags)?;
            write_args(f, args)?;
            write_usize_targets(f, targets)?;
        }
        RawPassthrough::Noise {
            name,
            tags,
            args,
            targets,
            ..
        } => {
            f.write_str(name.canonical_name())?;
            write_tags(f, tags)?;
            write_args(f, args)?;
            write_usize_targets(f, targets)?;
        }
        RawPassthrough::Measure {
            name,
            tags,
            args,
            targets,
            ..
        } => {
            f.write_str(name.canonical_name())?;
            write_tags(f, tags)?;
            write_args(f, args)?;
            write_usize_targets(f, targets)?;
        }
        RawPassthrough::Annotation {
            kind,
            args,
            targets,
            ..
        } => {
            f.write_str(kind.canonical_name())?;
            write_args(f, args)?;
            write_usize_targets(f, targets)?;
        }
    }
    writeln!(f)
}

fn fmt_ext(i: &ExtendedInstruction, f: &mut fmt::Formatter<'_>, depth: usize) -> fmt::Result {
    match i {
        ExtendedInstruction::Raw(r) => return fmt_raw_passthrough(r, f, depth),
        ExtendedInstruction::Repeat { count, body, .. } => {
            write_indent(f, depth)?;
            writeln!(f, "REPEAT {count} {{")?;
            fmt_ext_slice(body, f, depth + 1)?;
            write_indent(f, depth)?;
            return writeln!(f, "}}");
        }
        _ => {}
    }

    write_indent(f, depth)?;
    match i {
        ExtendedInstruction::Raw(_) | ExtendedInstruction::Repeat { .. } => unreachable!(),
        ExtendedInstruction::T { targets, .. } => {
            f.write_str("S[T]")?;
            write_usize_targets(f, targets)?;
        }
        ExtendedInstruction::TDag { targets, .. } => {
            f.write_str("S_DAG[T]")?;
            write_usize_targets(f, targets)?;
        }
        ExtendedInstruction::Rotation {
            axis,
            theta,
            targets,
            ..
        } => {
            let axis_tag = match axis {
                Axis::X => "R_X",
                Axis::Y => "R_Y",
                Axis::Z => "R_Z",
            };
            write!(f, "I[{}(theta={})]", axis_tag, FloatLit(*theta))?;
            write_usize_targets(f, targets)?;
        }
        ExtendedInstruction::U3 {
            theta,
            phi,
            lambda,
            targets,
            ..
        } => {
            write!(
                f,
                "I[U3(theta={}, phi={}, lambda={})]",
                FloatLit(*theta),
                FloatLit(*phi),
                FloatLit(*lambda),
            )?;
            write_usize_targets(f, targets)?;
        }
        ExtendedInstruction::Loss { p, targets, .. } => {
            write!(f, "I_ERROR[loss]({})", FloatLit(*p))?;
            write_usize_targets(f, targets)?;
        }
        ExtendedInstruction::CorrelatedLoss { ps, targets, .. } => {
            write!(
                f,
                "I_ERROR[correlated_loss]({}, {}, {})",
                FloatLit(ps[0]),
                FloatLit(ps[1]),
                FloatLit(ps[2]),
            )?;
            for &(a, b) in targets {
                write!(f, " {a} {b}")?;
            }
        }
        ExtendedInstruction::MPad {
            tags, prob, bits, ..
        } => {
            f.write_str("MPAD")?;
            write_tags(f, tags)?;
            if let Some(p) = prob {
                write!(f, "({})", FloatLit(*p))?;
            }
            for &bit in bits {
                write!(f, " {}", u8::from(bit))?;
            }
        }
    }
    writeln!(f)
}

fn write_tags(f: &mut fmt::Formatter<'_>, tags: &[Tag]) -> fmt::Result {
    if tags.is_empty() {
        return Ok(());
    }
    f.write_str("[")?;
    for (i, tag) in tags.iter().enumerate() {
        if i > 0 {
            f.write_str(", ")?;
        }
        f.write_str(&tag.name)?;
        if !tag.params.is_empty() {
            f.write_str("(")?;
            for (j, p) in tag.params.iter().enumerate() {
                if j > 0 {
                    f.write_str(", ")?;
                }
                match p {
                    TagParam::Positional(v) => write!(f, "{}", FloatLit(*v))?,
                    TagParam::Named { key, value } => {
                        write!(f, "{key}={}", FloatLit(*value))?;
                    }
                }
            }
            f.write_str(")")?;
        }
    }
    f.write_str("]")
}

fn write_args(f: &mut fmt::Formatter<'_>, args: &[f64]) -> fmt::Result {
    if args.is_empty() {
        return Ok(());
    }
    f.write_str("(")?;
    for (i, a) in args.iter().enumerate() {
        if i > 0 {
            f.write_str(", ")?;
        }
        write!(f, "{}", FloatLit(*a))?;
    }
    f.write_str(")")
}

fn write_usize_targets(f: &mut fmt::Formatter<'_>, targets: &[usize]) -> fmt::Result {
    for t in targets {
        write!(f, " {t}")?;
    }
    Ok(())
}

/// f64 formatter that always emits a decimal point. The grammar's
/// `signed_float` accepts bare integers too (`42` parses), so this is
/// purely a readability choice — printing `1.0` instead of `1` keeps the
/// canonical output looking like floating-point everywhere args are
/// expected.
struct FloatLit(f64);

impl fmt::Display for FloatLit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = format!("{}", self.0);
        if s.contains('.') || s.contains('e') || s.contains('E') || s.contains("inf") || s == "NaN"
        {
            f.write_str(&s)
        } else {
            write!(f, "{s}.0")
        }
    }
}
