// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! First-class canonical printer for both AST layers. Output is canonical
//! Stim — 4-space REPEAT indentation, `[tags](args) targets` ordering — and
//! round-trips: parse → print → parse is a fixpoint.

use std::borrow::Cow;
use std::fmt;

use crate::ast::shared::{
    AnnotationOp, Axis, GateOp, MeasureOp, MppOp, NoiseOp, PauliFactor, Tag, TagParam, Target,
};
use crate::ast::{ExtendedInstruction, ExtendedProgram, Instruction, Program};

pub struct PrintOptions {
    pub indent: Cow<'static, str>,
}

impl Default for PrintOptions {
    fn default() -> Self {
        PrintOptions {
            indent: Cow::Borrowed("    "),
        }
    }
}

pub trait StimPrint {
    fn print(&self, out: &mut dyn fmt::Write, opts: &PrintOptions, depth: usize) -> fmt::Result;

    fn to_stim(&self) -> String {
        self.to_stim_with(&PrintOptions::default())
    }

    fn to_stim_with(&self, opts: &PrintOptions) -> String {
        let mut s = String::new();
        let _ = self.print(&mut s, opts, 0);
        s
    }
}

// ---------------------------------------------------------------------------
// Low-level writers shared by every `StimPrint` impl. They produce canonical
// Stim byte-for-byte and write to any `&mut dyn fmt::Write`.
// ---------------------------------------------------------------------------

fn write_indent(out: &mut dyn fmt::Write, opts: &PrintOptions, depth: usize) -> fmt::Result {
    for _ in 0..depth {
        out.write_str(&opts.indent)?;
    }
    Ok(())
}

fn write_tags(out: &mut dyn fmt::Write, tags: &[Tag]) -> fmt::Result {
    if tags.is_empty() {
        return Ok(());
    }
    out.write_str("[")?;
    for (i, tag) in tags.iter().enumerate() {
        if i > 0 {
            out.write_str(", ")?;
        }
        out.write_str(&tag.name)?;
        if !tag.params.is_empty() {
            out.write_str("(")?;
            for (j, p) in tag.params.iter().enumerate() {
                if j > 0 {
                    out.write_str(", ")?;
                }
                match p {
                    TagParam::Positional(v) => write!(out, "{}", FloatLit(*v))?,
                    TagParam::Named { key, value } => {
                        write!(out, "{key}={}", FloatLit(*value))?;
                    }
                }
            }
            out.write_str(")")?;
        }
    }
    out.write_str("]")
}

fn write_args(out: &mut dyn fmt::Write, args: &[f64]) -> fmt::Result {
    if args.is_empty() {
        return Ok(());
    }
    out.write_str("(")?;
    for (i, a) in args.iter().enumerate() {
        if i > 0 {
            out.write_str(", ")?;
        }
        write!(out, "{}", FloatLit(*a))?;
    }
    out.write_str(")")
}

fn write_usize_targets(out: &mut dyn fmt::Write, targets: &[usize]) -> fmt::Result {
    for t in targets {
        write!(out, " {t}")?;
    }
    Ok(())
}

/// Print gate targets, rendering measurement-record controls as `rec[-k]`
/// so the output round-trips back through the parser.
fn write_targets(out: &mut dyn fmt::Write, targets: &[Target]) -> fmt::Result {
    for t in targets {
        match t {
            Target::Qubit(q) => write!(out, " {q}")?,
            Target::Rec(k) => write!(out, " rec[-{k}]")?,
        }
    }
    Ok(())
}

/// Print a `REPEAT count { … }` block, recursively printing the body one
/// indent level deeper and closing the brace at the block's own depth. The
/// caller is responsible for the trailing newline after the closing brace.
fn write_repeat_block<T: StimPrint>(
    out: &mut dyn fmt::Write,
    opts: &PrintOptions,
    depth: usize,
    count: u64,
    body: &[T],
) -> fmt::Result {
    writeln!(out, "REPEAT {count} {{")?;
    for instr in body {
        instr.print(out, opts, depth + 1)?;
    }
    write_indent(out, opts, depth)?;
    out.write_str("}")
}

/// Print `MPP` products as space-separated, `*`-joined Pauli factors
/// (`X0*Y1*Z2`) so the output round-trips back through the parser.
fn write_mpp_products(out: &mut dyn fmt::Write, products: &[Vec<PauliFactor>]) -> fmt::Result {
    for product in products {
        out.write_str(" ")?;
        for (i, factor) in product.iter().enumerate() {
            if i > 0 {
                out.write_str("*")?;
            }
            write!(out, "{}{}", factor.axis.as_char(), factor.qubit)?;
        }
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

// ---------------------------------------------------------------------------
// StimPrint for shared *Op structs
// ---------------------------------------------------------------------------

impl StimPrint for GateOp {
    fn print(&self, out: &mut dyn fmt::Write, _opts: &PrintOptions, _depth: usize) -> fmt::Result {
        out.write_str(self.name.canonical_name())?;
        write_tags(out, &self.tags)?;
        write_args(out, &self.args)?;
        write_targets(out, &self.targets)
    }
}

impl StimPrint for NoiseOp {
    fn print(&self, out: &mut dyn fmt::Write, _opts: &PrintOptions, _depth: usize) -> fmt::Result {
        out.write_str(self.name.canonical_name())?;
        write_tags(out, &self.tags)?;
        write_args(out, &self.args)?;
        write_usize_targets(out, &self.targets)
    }
}

impl StimPrint for MeasureOp {
    fn print(&self, out: &mut dyn fmt::Write, _opts: &PrintOptions, _depth: usize) -> fmt::Result {
        out.write_str(self.name.canonical_name())?;
        write_tags(out, &self.tags)?;
        write_args(out, &self.args)?;
        write_usize_targets(out, &self.targets)
    }
}

impl StimPrint for AnnotationOp {
    fn print(&self, out: &mut dyn fmt::Write, _opts: &PrintOptions, _depth: usize) -> fmt::Result {
        out.write_str(self.kind.canonical_name())?;
        write_args(out, &self.args)?;
        write_usize_targets(out, &self.targets)
    }
}

impl StimPrint for MppOp {
    fn print(&self, out: &mut dyn fmt::Write, _opts: &PrintOptions, _depth: usize) -> fmt::Result {
        out.write_str(crate::instructions::MeasureName::MPP.canonical_name())?;
        write_tags(out, &self.tags)?;
        write_args(out, &self.args)?;
        write_mpp_products(out, &self.products)
    }
}

// ---------------------------------------------------------------------------
// StimPrint for vanilla Instruction / Program
// ---------------------------------------------------------------------------

impl StimPrint for Instruction {
    fn print(&self, out: &mut dyn fmt::Write, opts: &PrintOptions, depth: usize) -> fmt::Result {
        write_indent(out, opts, depth)?;
        match self {
            Instruction::Gate(op) => op.print(out, opts, depth)?,
            Instruction::Noise(op) => op.print(out, opts, depth)?,
            Instruction::Measure(op) => op.print(out, opts, depth)?,
            Instruction::Annotation(op) => op.print(out, opts, depth)?,
            Instruction::Mpp(op) => op.print(out, opts, depth)?,
            Instruction::MPad {
                tags, prob, bits, ..
            } => {
                out.write_str("MPAD")?;
                write_tags(out, tags)?;
                if let Some(p) = prob {
                    write!(out, "({})", FloatLit(*p))?;
                }
                write_usize_targets(out, bits)?;
            }
            Instruction::Repeat { count, body, .. } => {
                write_repeat_block(out, opts, depth, *count, body)?;
            }
        }
        writeln!(out)
    }
}

impl StimPrint for Program {
    fn print(&self, out: &mut dyn fmt::Write, opts: &PrintOptions, depth: usize) -> fmt::Result {
        for instr in &self.instructions {
            instr.print(out, opts, depth)?;
        }
        Ok(())
    }
}

impl fmt::Display for Program {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.print(f, &PrintOptions::default(), 0)
    }
}

// ---------------------------------------------------------------------------
// StimPrint for ExtendedInstruction / ExtendedProgram
// ---------------------------------------------------------------------------

impl StimPrint for ExtendedInstruction {
    fn print(&self, out: &mut dyn fmt::Write, opts: &PrintOptions, depth: usize) -> fmt::Result {
        write_indent(out, opts, depth)?;
        match self {
            ExtendedInstruction::Gate(op) => op.print(out, opts, depth)?,
            ExtendedInstruction::Noise(op) => op.print(out, opts, depth)?,
            ExtendedInstruction::Measure(op) => op.print(out, opts, depth)?,
            ExtendedInstruction::Annotation(op) => op.print(out, opts, depth)?,
            ExtendedInstruction::Mpp(op) => op.print(out, opts, depth)?,
            ExtendedInstruction::T { targets, .. } => {
                out.write_str("S[T]")?;
                write_usize_targets(out, targets)?;
            }
            ExtendedInstruction::TDag { targets, .. } => {
                out.write_str("S_DAG[T]")?;
                write_usize_targets(out, targets)?;
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
                write!(out, "I[{}(theta={})]", axis_tag, FloatLit(*theta))?;
                write_usize_targets(out, targets)?;
            }
            ExtendedInstruction::U3 {
                theta,
                phi,
                lambda,
                targets,
                ..
            } => {
                write!(
                    out,
                    "I[U3(theta={}, phi={}, lambda={})]",
                    FloatLit(*theta),
                    FloatLit(*phi),
                    FloatLit(*lambda),
                )?;
                write_usize_targets(out, targets)?;
            }
            ExtendedInstruction::Loss { p, targets, .. } => {
                write!(out, "I_ERROR[loss]({})", FloatLit(*p))?;
                write_usize_targets(out, targets)?;
            }
            ExtendedInstruction::CorrelatedLoss { ps, targets, .. } => {
                write!(
                    out,
                    "I_ERROR[correlated_loss]({}, {}, {})",
                    FloatLit(ps[0]),
                    FloatLit(ps[1]),
                    FloatLit(ps[2]),
                )?;
                for &(a, b) in targets {
                    write!(out, " {a} {b}")?;
                }
            }
            ExtendedInstruction::MPad {
                tags, prob, bits, ..
            } => {
                out.write_str("MPAD")?;
                write_tags(out, tags)?;
                if let Some(p) = prob {
                    write!(out, "({})", FloatLit(*p))?;
                }
                for &bit in bits {
                    write!(out, " {}", u8::from(bit))?;
                }
            }
            ExtendedInstruction::Repeat { count, body, .. } => {
                write_repeat_block(out, opts, depth, *count, body)?;
            }
        }
        writeln!(out)
    }
}

impl StimPrint for ExtendedProgram {
    fn print(&self, out: &mut dyn fmt::Write, opts: &PrintOptions, depth: usize) -> fmt::Result {
        for instr in &self.instructions {
            instr.print(out, opts, depth)?;
        }
        Ok(())
    }
}

impl fmt::Display for ExtendedProgram {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.print(f, &PrintOptions::default(), 0)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use crate::print::StimPrint;
    use crate::{parse, parse_extended};

    #[test]
    fn vanilla_printed_form_is_canonical_shape() {
        let src = "H 0  # trail\nCX  0   1\nDEPOLARIZE1(0.05) 0 1\nREPEAT 2 { X 0 }\n";
        let ast = parse(src).unwrap();
        let expected = "H 0\nCX 0 1\nDEPOLARIZE1(0.05) 0 1\nREPEAT 2 {\n    X 0\n}\n";
        assert_eq!(ast.to_stim(), expected);
        assert_eq!(format!("{ast}"), expected);
    }

    #[test]
    fn extended_printed_form_lowers_sugar_into_canonical_stim() {
        let src = "S[T] 0\nI[R_X(theta=0.25)] 1\nI_ERROR[loss](0.01) 2\n";
        let ast = parse_extended(src).unwrap();
        let expected = "S[T] 0\nI[R_X(theta=0.25)] 1\nI_ERROR[loss](0.01) 2\n";
        assert_eq!(ast.to_stim(), expected);
    }

    #[test]
    fn rec_and_mpp_targets_round_trip() {
        // rec[-k] feed-forward control and MPP Pauli products print canonically.
        let ast = parse("CX rec[-1] 0\nMPP X0*Y1*Z2\n").unwrap();
        assert_eq!(ast.to_stim(), "CX rec[-1] 0\nMPP X0*Y1*Z2\n");
    }
}
