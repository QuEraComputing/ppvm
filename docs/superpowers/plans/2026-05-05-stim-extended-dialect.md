# Stim Extended Dialect Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move recognition of PPVM tag-based extensions (`S[T]`, `I[R_X/R_Y/R_Z/U3]`, `I_ERROR[loss/correlated_loss]`) from `ppvm-stim/src/normalize.rs` into a new `stim-parser::extended` module that produces a typed `ExtendedProgram`, and migrate `ppvm-stim` to consume it.

**Architecture:** Two-layer parse: vanilla `parse()` returns `Program` (unchanged); new `parse_extended()` calls `parse()` then runs a post-pass `interpret(&Program) -> Result<ExtendedProgram, ExtendedParseError>` that promotes recognized tag-bearing instructions to first-class `ExtendedInstruction` variants. The four "strict" hosts (`S`, `S_DAG`, `I`, `I_ERROR`) reject unrecognized or malformed tags at parse time; all other instructions are lenient pass-through with their tags retained. `ppvm-stim/src/normalize.rs` switches from `&Program` to `&ExtendedProgram` and becomes a near-1:1 translation pass; the `find_tag` / `require_no_params` / `identity_to_kind` helpers and `NormalizeError::InvalidTag` are deleted.

**Tech Stack:** Rust 2024 edition, chumsky 0.12 (already wired), thiserror 2 (already a dep of `stim-parser`), workspace test runner via `cargo test`. Rustfmt before every commit.

**Reference spec:** [docs/superpowers/specs/2026-05-05-stim-extended-dialect-design.md](../specs/2026-05-05-stim-extended-dialect-design.md)

---

## Task 1: Vanilla pass-through scaffold

Create the AST, error type, interpret pass-through (no extensions yet), `parse_extended` entry point, and the integration-test file. Verify pass-through works for all non-extension instruction kinds. Strict-host recognition arrives in tasks 2–6.

**Files:**
- Create: `crates/stim-parser/src/extended/ast.rs`
- Modify: `crates/stim-parser/src/extended/parser.rs` (currently empty)
- Create: `crates/stim-parser/src/extended/interpret.rs`
- Modify: `crates/stim-parser/src/extended/mod.rs` (currently empty)
- Modify: `crates/stim-parser/src/lib.rs:1-15`
- Create: `crates/stim-parser/tests/extended.rs`

- [ ] **Step 1: Write the failing tests**

Create `crates/stim-parser/tests/extended.rs`:

```rust
use stim_parser::extended::{
    Axis, ExtendedInstruction, ExtendedParseError, ExtendedProgram, parse_extended,
};
use stim_parser::ast::{GateName, MeasureName, NoiseName, AnnotationKind};

fn parse_ok(src: &str) -> ExtendedProgram {
    parse_extended(src).expect("parse_extended")
}

#[test]
fn vanilla_h_passes_through() {
    let p = parse_ok("H 0\n");
    assert_eq!(p.instructions.len(), 1);
    match &p.instructions[0] {
        ExtendedInstruction::Gate { name, tags, targets, line, .. } => {
            assert_eq!(*name, GateName::H);
            assert!(tags.is_empty());
            assert_eq!(targets, &vec![0]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn vanilla_measure_passes_through() {
    let p = parse_ok("M 0 1\n");
    match &p.instructions[0] {
        ExtendedInstruction::Measure { name, targets, .. } => {
            assert_eq!(*name, MeasureName::M);
            assert_eq!(targets, &vec![0, 1]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn vanilla_depolarize1_noise_passes_through() {
    let p = parse_ok("DEPOLARIZE1(0.01) 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Noise { name, args, .. } => {
            assert_eq!(*name, NoiseName::Depolarize1);
            assert_eq!(args, &vec![0.01]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn vanilla_annotation_passes_through() {
    let p = parse_ok("TICK\n");
    match &p.instructions[0] {
        ExtendedInstruction::Annotation { kind, .. } => {
            assert_eq!(*kind, AnnotationKind::Tick);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn vanilla_mpad_passes_through() {
    let p = parse_ok("MPAD 0 1\n");
    match &p.instructions[0] {
        ExtendedInstruction::MPad { bits, prob, .. } => {
            assert_eq!(bits, &vec![0, 1]);
            assert!(prob.is_none());
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn repeat_recurses_into_body() {
    let p = parse_ok("REPEAT 3 {\n    H 0\n}\n");
    match &p.instructions[0] {
        ExtendedInstruction::Repeat { count, body, .. } => {
            assert_eq!(*count, 3);
            assert_eq!(body.len(), 1);
            assert!(matches!(
                &body[0],
                ExtendedInstruction::Gate { name: GateName::H, .. }
            ));
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn lenient_unknown_tag_on_h_passes_through() {
    let p = parse_ok("H[unrelated] 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Gate { name, tags, .. } => {
            assert_eq!(*name, GateName::H);
            assert_eq!(tags.len(), 1);
            assert_eq!(tags[0].name, "unrelated");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_error_propagates() {
    let err = parse_extended("FROBNICATE 0\n").unwrap_err();
    assert!(matches!(err, ExtendedParseError::Parse(_)));
}

#[test]
fn axis_enum_has_xyz() {
    // Compile-time check that Axis exists with the three variants.
    let _x = Axis::X;
    let _y = Axis::Y;
    let _z = Axis::Z;
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p stim-parser --test extended`
Expected: FAIL with errors about `extended` module / `parse_extended` / `ExtendedProgram` not found.

- [ ] **Step 3: Write the AST types**

Create `crates/stim-parser/src/extended/ast.rs`:

```rust
//! Typed AST for Stim with PPVM tag-based extensions promoted to
//! first-class instruction variants.

use crate::ast::{AnnotationKind, GateName, MeasureName, NoiseName, Tag};

#[derive(Debug, Clone, PartialEq)]
pub struct ExtendedProgram {
    pub instructions: Vec<ExtendedInstruction>,
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ExtendedInstruction {
    // Vanilla pass-through. Tags retained on lenient hosts.
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
    MPad {
        tags: Vec<Tag>,
        prob: Option<f64>,
        bits: Vec<usize>,
        line: usize,
    },
    Repeat {
        count: u64,
        body: Vec<ExtendedInstruction>,
        line: usize,
    },

    // Promoted from extension tags.
    T { targets: Vec<usize>, line: usize },
    TDag { targets: Vec<usize>, line: usize },
    Rotation { axis: Axis, theta: f64, targets: Vec<usize>, line: usize },
    U3 { theta: f64, phi: f64, lambda: f64, targets: Vec<usize>, line: usize },
    Loss { p: f64, targets: Vec<usize>, line: usize },
    CorrelatedLoss { ps: [f64; 3], targets: Vec<usize>, line: usize },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Axis { X, Y, Z }
```

- [ ] **Step 4: Write the error type and `parse_extended` entry point**

Create `crates/stim-parser/src/extended/parser.rs`:

```rust
//! Public entry point for the extended-dialect parse.

use crate::ast::ParseError;
use crate::extended::ast::ExtendedProgram;
use crate::extended::interpret::interpret;

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ExtendedParseError {
    #[error(transparent)]
    Parse(#[from] ParseError),
    #[error("invalid tag '{tag}' on '{instruction}' at line {line}: {message}")]
    InvalidTag {
        tag: String,
        instruction: String,
        line: usize,
        message: String,
    },
}

/// Parse Stim source and interpret PPVM tag-based extensions into a typed
/// [`ExtendedProgram`]. Strict for the four extension hosts (`S`, `S_DAG`,
/// `I`, `I_ERROR`); lenient pass-through for everything else.
pub fn parse_extended(src: &str) -> Result<ExtendedProgram, ExtendedParseError> {
    let prog = crate::parser::parse(src)?;
    interpret(&prog)
}
```

- [ ] **Step 5: Write the interpret pass with vanilla pass-through only**

Create `crates/stim-parser/src/extended/interpret.rs`:

```rust
//! Post-pass interpretation: walks a vanilla `Program` and produces an
//! `ExtendedProgram`. Strict for hosts `S`, `S_DAG`, `I`, `I_ERROR`;
//! lenient pass-through for everything else.

use crate::ast::{Program, RawInstruction};
use crate::extended::ast::{ExtendedInstruction, ExtendedProgram};
use crate::extended::parser::ExtendedParseError;

pub(crate) fn interpret(prog: &Program) -> Result<ExtendedProgram, ExtendedParseError> {
    let mut out = Vec::with_capacity(prog.instructions.len());
    interpret_slice(&prog.instructions, &mut out)?;
    Ok(ExtendedProgram { instructions: out })
}

fn interpret_slice(
    src: &[RawInstruction],
    out: &mut Vec<ExtendedInstruction>,
) -> Result<(), ExtendedParseError> {
    for raw in src {
        out.push(interpret_one(raw)?);
    }
    Ok(())
}

fn interpret_one(raw: &RawInstruction) -> Result<ExtendedInstruction, ExtendedParseError> {
    match raw {
        RawInstruction::Gate { name, tags, args, targets, line } => {
            // Strict-host recognition is added in later tasks; for now,
            // every gate is a vanilla pass-through.
            Ok(ExtendedInstruction::Gate {
                name: *name,
                tags: tags.clone(),
                args: args.clone(),
                targets: targets.clone(),
                line: *line,
            })
        }
        RawInstruction::Noise { name, tags, args, targets, line } => {
            Ok(ExtendedInstruction::Noise {
                name: *name,
                tags: tags.clone(),
                args: args.clone(),
                targets: targets.clone(),
                line: *line,
            })
        }
        RawInstruction::Measure { name, tags, args, targets, line } => {
            Ok(ExtendedInstruction::Measure {
                name: *name,
                tags: tags.clone(),
                args: args.clone(),
                targets: targets.clone(),
                line: *line,
            })
        }
        RawInstruction::Annotation { kind, args, targets, line } => {
            Ok(ExtendedInstruction::Annotation {
                kind: *kind,
                args: args.clone(),
                targets: targets.clone(),
                line: *line,
            })
        }
        RawInstruction::MPad { tags, prob, bits, line } => {
            Ok(ExtendedInstruction::MPad {
                tags: tags.clone(),
                prob: *prob,
                bits: bits.clone(),
                line: *line,
            })
        }
        RawInstruction::Repeat { count, body, line } => {
            let mut inner = Vec::with_capacity(body.len());
            interpret_slice(body, &mut inner)?;
            Ok(ExtendedInstruction::Repeat {
                count: *count,
                body: inner,
                line: *line,
            })
        }
    }
}
```

- [ ] **Step 6: Wire `extended/mod.rs` and lib.rs**

Replace `crates/stim-parser/src/extended/mod.rs` with:

```rust
//! Extended Stim dialect — interprets PPVM tag-based extensions into a
//! typed AST.

pub mod ast;
mod interpret;
pub mod parser;

pub use ast::{Axis, ExtendedInstruction, ExtendedProgram};
pub use parser::{ExtendedParseError, parse_extended};
```

In `crates/stim-parser/src/lib.rs`, add `pub mod extended;` to the module list. The file currently reads:

```rust
pub mod ast;
mod grammar;
mod line_map;
mod parser;
mod table;

use line_map::LineMap;

pub mod prelude {
    pub use crate::ast::{
        AnnotationKind, GateName, MeasureName, NoiseName, ParseError, Program, RawInstruction, Tag,
        TagParam,
    };
    pub use crate::parser::parse;
}
```

Update to:

```rust
pub mod ast;
pub mod extended;
mod grammar;
mod line_map;
mod parser;
mod table;

use line_map::LineMap;

pub mod prelude {
    pub use crate::ast::{
        AnnotationKind, GateName, MeasureName, NoiseName, ParseError, Program, RawInstruction, Tag,
        TagParam,
    };
    pub use crate::parser::parse;
}
```

(Prelude exports for the extended types come in Task 7 — keep this commit small.)

- [ ] **Step 7: Run tests to verify they pass**

Run: `cargo test -p stim-parser --test extended`
Expected: PASS — 9 tests pass.

Then run: `cargo check --workspace`
Expected: clean.

- [ ] **Step 8: Format and commit**

```bash
cargo fmt -p stim-parser
git add crates/stim-parser/src/extended/ crates/stim-parser/src/lib.rs crates/stim-parser/tests/extended.rs
git commit -m "$(cat <<'EOF'
Add stim-parser extended-dialect scaffold

Adds the extended/ module with ExtendedProgram, ExtendedInstruction
(Flavor B with promoted variants), Axis, ExtendedParseError, and a
parse_extended entry point that runs an interpret post-pass over a
vanilla Program. This commit only handles vanilla pass-through; the
strict-host recognition for S[T], I[R_*], I[U3], and I_ERROR[*] arrives
in subsequent commits.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Recognize `S[T]` and `S_DAG[T]`

Add the simplest extension shapes: empty-param `[T]` tag on `S` or `S_DAG` promotes to a `T` / `TDag` instruction. Strict — any other tag on these hosts errors.

**Files:**
- Modify: `crates/stim-parser/src/extended/interpret.rs`
- Modify: `crates/stim-parser/tests/extended.rs`

- [ ] **Step 1: Write the failing tests**

Append to `crates/stim-parser/tests/extended.rs`:

```rust
fn parse_err(src: &str) -> ExtendedParseError {
    parse_extended(src).expect_err("must reject")
}

#[test]
fn s_t_promotes_to_t() {
    let p = parse_ok("S[T] 0 1\n");
    match &p.instructions[0] {
        ExtendedInstruction::T { targets, line } => {
            assert_eq!(targets, &vec![0, 1]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn s_dag_t_promotes_to_t_dag() {
    let p = parse_ok("S_DAG[T] 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::TDag { targets, line } => {
            assert_eq!(targets, &vec![0]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn s_with_no_tag_is_vanilla_gate() {
    let p = parse_ok("S 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Gate { name, tags, .. } => {
            assert_eq!(*name, GateName::S);
            assert!(tags.is_empty());
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn s_with_unknown_tag_errors() {
    let err = parse_err("S[X] 0\n");
    match err {
        ExtendedParseError::InvalidTag { tag, instruction, line, .. } => {
            assert_eq!(tag, "X");
            assert_eq!(instruction, "S");
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn s_dag_with_unknown_tag_errors() {
    let err = parse_err("S_DAG[X] 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn s_t_with_params_errors() {
    // [T] must have no params.
    let err = parse_err("S[T(0.5)] 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p stim-parser --test extended s_`
Expected: FAIL — promotion tests fail because S[T] currently passes through as a vanilla `Gate` with the tag retained.

- [ ] **Step 3: Implement S/S_DAG recognition**

In `crates/stim-parser/src/extended/interpret.rs`, replace the `RawInstruction::Gate { … }` arm in `interpret_one` with a delegation to a new `interpret_gate` helper, and add the helper plus a `require_no_params` helper. Add `use` lines at the top for `GateName` and `Tag`.

Top-of-file imports:

```rust
use crate::ast::{GateName, Program, RawInstruction, Tag};
use crate::extended::ast::{ExtendedInstruction, ExtendedProgram};
use crate::extended::parser::ExtendedParseError;
```

Replace the `Gate` arm in `interpret_one`:

```rust
        RawInstruction::Gate { name, tags, args, targets, line } => {
            interpret_gate(*name, tags, args, targets, *line)
        }
```

Add `interpret_gate` and `require_no_params` after `interpret_one`:

```rust
fn interpret_gate(
    name: GateName,
    tags: &[Tag],
    args: &[f64],
    targets: &[usize],
    line: usize,
) -> Result<ExtendedInstruction, ExtendedParseError> {
    use GateName::*;
    match (name, tags) {
        // S[T] / S_DAG[T] — empty-param tag promotes to T / TDag.
        (S, [t]) if t.name == "T" => {
            require_no_params(t, "S", line)?;
            Ok(ExtendedInstruction::T { targets: targets.to_vec(), line })
        }
        (SDag, [t]) if t.name == "T" => {
            require_no_params(t, "S_DAG", line)?;
            Ok(ExtendedInstruction::TDag { targets: targets.to_vec(), line })
        }

        // Strict-host rejection for unknown tags on S / S_DAG.
        (S, [t]) => Err(ExtendedParseError::InvalidTag {
            tag: t.name.clone(),
            instruction: "S".into(),
            line,
            message: "expected [T]".into(),
        }),
        (SDag, [t]) => Err(ExtendedParseError::InvalidTag {
            tag: t.name.clone(),
            instruction: "S_DAG".into(),
            line,
            message: "expected [T]".into(),
        }),
        (S, _) | (SDag, _) if !tags.is_empty() => Err(ExtendedParseError::InvalidTag {
            tag: tags[0].name.clone(),
            instruction: if matches!(name, S) { "S" } else { "S_DAG" }.into(),
            line,
            message: "expected exactly one tag".into(),
        }),

        // Lenient pass-through for everything else (including S/S_DAG with no tag).
        _ => Ok(ExtendedInstruction::Gate {
            name,
            tags: tags.to_vec(),
            args: args.to_vec(),
            targets: targets.to_vec(),
            line,
        }),
    }
}

fn require_no_params(
    tag: &Tag,
    instruction: &str,
    line: usize,
) -> Result<(), ExtendedParseError> {
    if !tag.params.is_empty() {
        return Err(ExtendedParseError::InvalidTag {
            tag: tag.name.clone(),
            instruction: instruction.to_string(),
            line,
            message: "tag must have no parameters".into(),
        });
    }
    Ok(())
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p stim-parser --test extended`
Expected: PASS — all tests so far (Task 1 + Task 2) pass.

- [ ] **Step 5: Format and commit**

```bash
cargo fmt -p stim-parser
git add crates/stim-parser/src/extended/interpret.rs crates/stim-parser/tests/extended.rs
git commit -m "$(cat <<'EOF'
Recognize S[T] and S_DAG[T] in extended dialect

Promotes [T]-tagged S and S_DAG to ExtendedInstruction::T and ::TDag.
Unknown tags on S / S_DAG produce InvalidTag at parse time.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Recognize `I[R_X]`, `I[R_Y]`, `I[R_Z]`

Identity with a single-axis rotation tag promotes to `Rotation`. Each tag requires a named `theta` parameter.

**Files:**
- Modify: `crates/stim-parser/src/extended/interpret.rs`
- Modify: `crates/stim-parser/tests/extended.rs`

- [ ] **Step 1: Write the failing tests**

Append to `crates/stim-parser/tests/extended.rs`:

```rust
fn approx_eq(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} vs {b}");
}

#[test]
fn i_r_x_promotes_to_rotation_x() {
    let p = parse_ok("I[R_X(theta=0.5*pi)] 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Rotation { axis, theta, targets, line } => {
            assert!(matches!(axis, Axis::X));
            approx_eq(*theta, 0.5 * std::f64::consts::PI);
            assert_eq!(targets, &vec![0]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_r_y_promotes_to_rotation_y() {
    let p = parse_ok("I[R_Y(theta=0.25)] 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Rotation { axis, theta, .. } => {
            assert!(matches!(axis, Axis::Y));
            approx_eq(*theta, 0.25);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_r_z_promotes_to_rotation_z() {
    let p = parse_ok("I[R_Z(theta=0.1)] 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Rotation { axis, theta, .. } => {
            assert!(matches!(axis, Axis::Z));
            approx_eq(*theta, 0.1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_with_no_tag_is_vanilla_identity() {
    let p = parse_ok("I 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Gate { name, tags, .. } => {
            assert_eq!(*name, GateName::Identity);
            assert!(tags.is_empty());
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_r_x_missing_theta_errors() {
    let err = parse_err("I[R_X] 0\n");
    match err {
        ExtendedParseError::InvalidTag { tag, instruction, line, .. } => {
            assert_eq!(tag, "R_X");
            assert_eq!(instruction, "I");
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p stim-parser --test extended i_r_`
Expected: FAIL — `I[R_X]` etc. currently pass through as vanilla `Gate { name: Identity, tags: [...] }`.

- [ ] **Step 3: Implement Identity-with-rotation-tag recognition**

In `crates/stim-parser/src/extended/interpret.rs`, add `Axis` to the imports:

```rust
use crate::extended::ast::{Axis, ExtendedInstruction, ExtendedProgram};
```

In `interpret_gate`, add an `(Identity, [t])` arm above the lenient pass-through. Insert after the `(S, _) | (SDag, _) if !tags.is_empty() => …` arm:

```rust
        // I[<extension-tag>] dispatch. Strict: any unrecognized single
        // tag on Identity is rejected. Bare I and multi-tag I are
        // handled below (lenient / strict respectively).
        (Identity, [t]) => interpret_identity_tag(t, targets, line),
        (Identity, _) if !tags.is_empty() => Err(ExtendedParseError::InvalidTag {
            tag: tags[0].name.clone(),
            instruction: "I".into(),
            line,
            message: "expected exactly one tag".into(),
        }),
```

Add `interpret_identity_tag` and `lookup_named` helpers after `require_no_params`:

```rust
fn interpret_identity_tag(
    tag: &Tag,
    targets: &[usize],
    line: usize,
) -> Result<ExtendedInstruction, ExtendedParseError> {
    let axis = match tag.name.as_str() {
        "R_X" => Some(Axis::X),
        "R_Y" => Some(Axis::Y),
        "R_Z" => Some(Axis::Z),
        _ => None,
    };

    if let Some(axis) = axis {
        let theta = lookup_named(tag, "theta", "I", line)?;
        return Ok(ExtendedInstruction::Rotation {
            axis,
            theta,
            targets: targets.to_vec(),
            line,
        });
    }

    // U3 lands here in the next task. Until then, anything that isn't a
    // rotation is rejected.
    Err(ExtendedParseError::InvalidTag {
        tag: tag.name.clone(),
        instruction: "I".into(),
        line,
        message: "unrecognized tag (expected R_X / R_Y / R_Z)".into(),
    })
}

fn lookup_named(
    tag: &Tag,
    key: &str,
    instruction: &str,
    line: usize,
) -> Result<f64, ExtendedParseError> {
    use crate::ast::TagParam;
    tag.params
        .iter()
        .find_map(|p| match p {
            TagParam::Named { key: k, value } if k == key => Some(*value),
            _ => None,
        })
        .ok_or(ExtendedParseError::InvalidTag {
            tag: tag.name.clone(),
            instruction: instruction.to_string(),
            line,
            message: format!("missing required named parameter '{key}'"),
        })
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p stim-parser --test extended`
Expected: PASS.

- [ ] **Step 5: Format and commit**

```bash
cargo fmt -p stim-parser
git add crates/stim-parser/src/extended/interpret.rs crates/stim-parser/tests/extended.rs
git commit -m "$(cat <<'EOF'
Recognize I[R_X/R_Y/R_Z] in extended dialect

Promotes I with an R_X / R_Y / R_Z tag and a named 'theta' parameter
to ExtendedInstruction::Rotation. Bare I passes through as vanilla
Identity. I with an unrecognized single tag, multiple tags, or
R_* missing 'theta' produces InvalidTag.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Recognize `I[U3(theta, phi, lambda)]`

Identity with a U3 tag carrying three named params promotes to `U3`.

**Files:**
- Modify: `crates/stim-parser/src/extended/interpret.rs`
- Modify: `crates/stim-parser/tests/extended.rs`

- [ ] **Step 1: Write the failing tests**

Append to `crates/stim-parser/tests/extended.rs`:

```rust
#[test]
fn i_u3_promotes_to_u3() {
    let p = parse_ok("I[U3(theta=0.1, phi=0.2, lambda=0.3)] 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::U3 { theta, phi, lambda, targets, line } => {
            approx_eq(*theta, 0.1);
            approx_eq(*phi, 0.2);
            approx_eq(*lambda, 0.3);
            assert_eq!(targets, &vec![0]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_u3_missing_phi_errors() {
    let err = parse_err("I[U3(theta=0.1, lambda=0.2)] 0\n");
    match err {
        ExtendedParseError::InvalidTag { tag, instruction, .. } => {
            assert_eq!(tag, "U3");
            assert_eq!(instruction, "I");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_unrecognized_tag_errors() {
    let err = parse_err("I[FOO] 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p stim-parser --test extended i_u3`
Expected: FAIL — `I[U3(...)]` currently rejected by Task 3's "unrecognized tag" arm.

- [ ] **Step 3: Implement U3 recognition**

In `crates/stim-parser/src/extended/interpret.rs`, modify `interpret_identity_tag` to handle `U3` before the unrecognized-tag fallback:

```rust
fn interpret_identity_tag(
    tag: &Tag,
    targets: &[usize],
    line: usize,
) -> Result<ExtendedInstruction, ExtendedParseError> {
    let axis = match tag.name.as_str() {
        "R_X" => Some(Axis::X),
        "R_Y" => Some(Axis::Y),
        "R_Z" => Some(Axis::Z),
        _ => None,
    };

    if let Some(axis) = axis {
        let theta = lookup_named(tag, "theta", "I", line)?;
        return Ok(ExtendedInstruction::Rotation {
            axis,
            theta,
            targets: targets.to_vec(),
            line,
        });
    }

    if tag.name == "U3" {
        let theta = lookup_named(tag, "theta", "I", line)?;
        let phi = lookup_named(tag, "phi", "I", line)?;
        let lambda = lookup_named(tag, "lambda", "I", line)?;
        return Ok(ExtendedInstruction::U3 {
            theta,
            phi,
            lambda,
            targets: targets.to_vec(),
            line,
        });
    }

    Err(ExtendedParseError::InvalidTag {
        tag: tag.name.clone(),
        instruction: "I".into(),
        line,
        message: "unrecognized tag (expected R_X / R_Y / R_Z / U3)".into(),
    })
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p stim-parser --test extended`
Expected: PASS.

- [ ] **Step 5: Format and commit**

```bash
cargo fmt -p stim-parser
git add crates/stim-parser/src/extended/interpret.rs crates/stim-parser/tests/extended.rs
git commit -m "$(cat <<'EOF'
Recognize I[U3(theta, phi, lambda)] in extended dialect

Promotes I with a U3 tag carrying three named parameters to
ExtendedInstruction::U3. Missing any of theta / phi / lambda produces
InvalidTag.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Recognize `I_ERROR[loss]`

`I_ERROR` with a `[loss]` tag and exactly one instruction-arg promotes to `Loss`.

**Files:**
- Modify: `crates/stim-parser/src/extended/interpret.rs`
- Modify: `crates/stim-parser/tests/extended.rs`

- [ ] **Step 1: Write the failing tests**

Append to `crates/stim-parser/tests/extended.rs`:

```rust
#[test]
fn i_error_loss_promotes_to_loss() {
    let p = parse_ok("I_ERROR[loss](0.01) 0\n");
    match &p.instructions[0] {
        ExtendedInstruction::Loss { p, targets, line } => {
            approx_eq(*p, 0.01);
            assert_eq!(targets, &vec![0]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_error_with_no_tag_errors() {
    let err = parse_err("I_ERROR(0.1) 0\n");
    match err {
        ExtendedParseError::InvalidTag { instruction, .. } => {
            assert_eq!(instruction, "I_ERROR");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_error_loss_wrong_arg_count_errors() {
    let err = parse_err("I_ERROR[loss](0.1, 0.2) 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}

#[test]
fn i_error_unknown_tag_errors() {
    let err = parse_err("I_ERROR[bogus](0.1) 0\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p stim-parser --test extended i_error`
Expected: FAIL — `I_ERROR` currently passes through as vanilla `Noise`.

- [ ] **Step 3: Implement I_ERROR / loss recognition**

In `crates/stim-parser/src/extended/interpret.rs`, add `NoiseName` to imports:

```rust
use crate::ast::{GateName, NoiseName, Program, RawInstruction, Tag};
```

In `interpret_one`, replace the `Noise` arm with delegation:

```rust
        RawInstruction::Noise { name, tags, args, targets, line } => {
            interpret_noise(*name, tags, args, targets, *line)
        }
```

Add `interpret_noise` after `interpret_identity_tag`:

```rust
fn interpret_noise(
    name: NoiseName,
    tags: &[Tag],
    args: &[f64],
    targets: &[usize],
    line: usize,
) -> Result<ExtendedInstruction, ExtendedParseError> {
    use NoiseName::*;
    match (name, tags) {
        // I_ERROR[loss] — 1 arg required.
        (IError, [t]) if t.name == "loss" => {
            require_no_params(t, "I_ERROR", line)?;
            if args.len() != 1 {
                return Err(ExtendedParseError::InvalidTag {
                    tag: "loss".into(),
                    instruction: "I_ERROR".into(),
                    line,
                    message: format!("[loss] expects 1 arg, got {}", args.len()),
                });
            }
            Ok(ExtendedInstruction::Loss {
                p: args[0],
                targets: targets.to_vec(),
                line,
            })
        }

        // I_ERROR with no tag is rejected — there is no defined meaning
        // for an untagged I_ERROR in PPVM.
        (IError, []) => Err(ExtendedParseError::InvalidTag {
            tag: String::new(),
            instruction: "I_ERROR".into(),
            line,
            message: "I_ERROR requires a [loss] or [correlated_loss] tag".into(),
        }),

        // I_ERROR with a single unrecognized tag.
        (IError, [t]) => Err(ExtendedParseError::InvalidTag {
            tag: t.name.clone(),
            instruction: "I_ERROR".into(),
            line,
            message: "expected [loss] or [correlated_loss]".into(),
        }),

        // I_ERROR with multiple tags.
        (IError, _) => Err(ExtendedParseError::InvalidTag {
            tag: tags[0].name.clone(),
            instruction: "I_ERROR".into(),
            line,
            message: "expected exactly one tag".into(),
        }),

        // Other noise channels: lenient pass-through.
        _ => Ok(ExtendedInstruction::Noise {
            name,
            tags: tags.to_vec(),
            args: args.to_vec(),
            targets: targets.to_vec(),
            line,
        }),
    }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p stim-parser --test extended`
Expected: PASS — note `i_error_correlated_loss_*` tests don't exist yet (Task 6).

- [ ] **Step 5: Format and commit**

```bash
cargo fmt -p stim-parser
git add crates/stim-parser/src/extended/interpret.rs crates/stim-parser/tests/extended.rs
git commit -m "$(cat <<'EOF'
Recognize I_ERROR[loss] in extended dialect

Promotes I_ERROR[loss](p) to ExtendedInstruction::Loss. Tagless
I_ERROR, I_ERROR with multiple tags, I_ERROR[unknown], and
I_ERROR[loss] with the wrong arg count all produce InvalidTag.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Recognize `I_ERROR[correlated_loss]`

`I_ERROR[correlated_loss]` with 1 or 3 instruction-args promotes to `CorrelatedLoss`. The 1-arg shorthand expands to `[p, 0.0, 0.0]`.

**Files:**
- Modify: `crates/stim-parser/src/extended/interpret.rs`
- Modify: `crates/stim-parser/tests/extended.rs`

- [ ] **Step 1: Write the failing tests**

Append to `crates/stim-parser/tests/extended.rs`:

```rust
#[test]
fn i_error_correlated_loss_one_arg_expands() {
    let p = parse_ok("I_ERROR[correlated_loss](0.5) 0 1\n");
    match &p.instructions[0] {
        ExtendedInstruction::CorrelatedLoss { ps, targets, line } => {
            approx_eq(ps[0], 0.5);
            approx_eq(ps[1], 0.0);
            approx_eq(ps[2], 0.0);
            assert_eq!(targets, &vec![0, 1]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_error_correlated_loss_three_args() {
    let p = parse_ok("I_ERROR[correlated_loss](0.1, 0.2, 0.3) 0 1\n");
    match &p.instructions[0] {
        ExtendedInstruction::CorrelatedLoss { ps, .. } => {
            approx_eq(ps[0], 0.1);
            approx_eq(ps[1], 0.2);
            approx_eq(ps[2], 0.3);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn i_error_correlated_loss_two_args_errors() {
    let err = parse_err("I_ERROR[correlated_loss](0.1, 0.2) 0 1\n");
    assert!(matches!(err, ExtendedParseError::InvalidTag { .. }));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p stim-parser --test extended i_error_correlated`
Expected: FAIL — `correlated_loss` currently hits the "expected [loss] or [correlated_loss]" arm in `interpret_noise` (because Task 5 only added a `loss` arm).

- [ ] **Step 3: Implement correlated_loss recognition**

In `crates/stim-parser/src/extended/interpret.rs`, add a `correlated_loss` arm to `interpret_noise` immediately after the `loss` arm (before the `(IError, [])` arm):

```rust
        // I_ERROR[correlated_loss] — 1 or 3 args; 1 arg expands to [p, 0, 0].
        (IError, [t]) if t.name == "correlated_loss" => {
            require_no_params(t, "I_ERROR", line)?;
            let ps = match args.len() {
                1 => [args[0], 0.0, 0.0],
                3 => [args[0], args[1], args[2]],
                n => {
                    return Err(ExtendedParseError::InvalidTag {
                        tag: "correlated_loss".into(),
                        instruction: "I_ERROR".into(),
                        line,
                        message: format!("[correlated_loss] expects 1 or 3 args, got {n}"),
                    });
                }
            };
            Ok(ExtendedInstruction::CorrelatedLoss {
                ps,
                targets: targets.to_vec(),
                line,
            })
        }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p stim-parser --test extended`
Expected: PASS — all extended tests now pass.

- [ ] **Step 5: Format and commit**

```bash
cargo fmt -p stim-parser
git add crates/stim-parser/src/extended/interpret.rs crates/stim-parser/tests/extended.rs
git commit -m "$(cat <<'EOF'
Recognize I_ERROR[correlated_loss] in extended dialect

Promotes I_ERROR[correlated_loss] with 1 or 3 instruction-args to
ExtendedInstruction::CorrelatedLoss. The 1-arg shorthand expands to
[p, 0, 0] (matching the existing normalize.rs behavior). Any other
arg count produces InvalidTag.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Wire prelude exports

Add the extended-dialect items to `stim_parser::prelude` so `ppvm-stim` (which already does `pub use stim_parser::prelude::*`) picks them up automatically.

**Files:**
- Modify: `crates/stim-parser/src/lib.rs:9-15`

- [ ] **Step 1: Write the failing test**

Append to `crates/stim-parser/tests/extended.rs`:

```rust
#[test]
fn prelude_exposes_parse_extended_and_types() {
    use stim_parser::prelude::*;
    let p: ExtendedProgram = parse_extended("H 0\n").unwrap();
    assert_eq!(p.instructions.len(), 1);
    // ExtendedParseError, ExtendedInstruction, Axis must also be in scope.
    fn _is_axis(_: Axis) {}
    fn _is_inst(_: &ExtendedInstruction) {}
    fn _is_err(_: ExtendedParseError) {}
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p stim-parser --test extended prelude_exposes`
Expected: FAIL with "cannot find type `ExtendedProgram` in scope" or similar.

- [ ] **Step 3: Update the prelude**

In `crates/stim-parser/src/lib.rs`, replace the `prelude` module:

```rust
pub mod prelude {
    pub use crate::ast::{
        AnnotationKind, GateName, MeasureName, NoiseName, ParseError, Program, RawInstruction, Tag,
        TagParam,
    };
    pub use crate::extended::{
        Axis, ExtendedInstruction, ExtendedParseError, ExtendedProgram, parse_extended,
    };
    pub use crate::parser::parse;
}
```

- [ ] **Step 4: Run test and full workspace check**

Run: `cargo test -p stim-parser --test extended prelude_exposes`
Expected: PASS.

Run: `cargo check --workspace`
Expected: clean — `ppvm-stim`'s `pub use stim_parser::prelude::*` now re-exports the extended items, but no consumer uses them yet.

- [ ] **Step 5: Format and commit**

```bash
cargo fmt -p stim-parser
git add crates/stim-parser/src/lib.rs crates/stim-parser/tests/extended.rs
git commit -m "$(cat <<'EOF'
Export extended-dialect items from stim_parser::prelude

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Migrate `ppvm-stim` to `ExtendedProgram` (atomic refactor)

Switch `to_tableau` from `&Program` to `&ExtendedProgram`, update `run_string` / `run_file` / `Error` in `lib.rs`, update the module-level doc-test, and update every test/bench file that calls `parse` directly. All in one atomic commit because the changes are interlocked: the moment `to_tableau`'s signature changes, every caller of `parse + to_tableau` must switch to `parse_extended` for the crate to compile.

This is the biggest task by line count but is mostly mechanical: every match arm in `normalize.rs` either stays the same (for vanilla pass-through) or gets simpler (for typed extensions); every test/bench `parse(...)` call becomes `parse_extended(...)`.

**Files:**
- Modify: `crates/ppvm-stim/src/normalize.rs` (full rewrite)
- Modify: `crates/ppvm-stim/src/lib.rs` (module doc, `Error` enum, `run_string` body)
- Modify: `crates/ppvm-stim/tests/normalize.rs` (helpers + delete moved tests)
- Modify: `crates/ppvm-stim/tests/executor.rs` (imports + `parse` → `parse_extended`)
- Modify: `crates/ppvm-stim/tests/run.rs` (imports + `Error::Parse` matching)
- Modify: `crates/ppvm-stim/tests/stim_corpus.rs` (imports + `parse` → `parse_extended`)
- Modify: `crates/ppvm-stim/benches/tableau-msd-stim.rs` (imports + `parse` → `parse_extended`)

- [ ] **Step 1: Update `tests/normalize.rs` helpers and delete moved tests**

Replace `crates/ppvm-stim/tests/normalize.rs` lines 1-13 (the imports and helpers) with:

```rust
use ppvm_stim::{
    GateKind, Instruction, MeasureKind, NoiseKind, NormalizeError, TableauProgram, normalize,
    parse_extended,
};

fn norm(src: &str) -> TableauProgram {
    let prog = parse_extended(src).expect("parse_extended");
    normalize::to_tableau(&prog).expect("normalize")
}

fn norm_err(src: &str) -> NormalizeError {
    let prog = parse_extended(src).expect("parse_extended");
    normalize::to_tableau(&prog).expect_err("must reject")
}

fn approx_eq(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} vs {b}");
}
```

Then delete the following tests from `crates/ppvm-stim/tests/normalize.rs` (now covered by `crates/stim-parser/tests/extended.rs`):

- `s_t_tag_maps_to_gate_t` (lines 31-38)
- `s_dag_t_tag_maps_to_gate_t_dag` (lines 40-47)
- `i_rx_tag_maps_to_rx` (lines 49-60)
- `i_u3_tag_maps_to_u3` (lines 62-75)
- `i_error_loss_tag_maps_to_loss` (lines 77-85)
- `i_error_correlated_loss_one_arg_expands_to_three` (lines 87-98)
- `i_error_correlated_loss_three_args_passthrough` (lines 100-110)
- `malformed_rx_tag_missing_theta_rejected` (lines 257-269)
- `untagged_i_error_rejected_as_invalid_tag` (lines 271-275)

Keep all other tests. The `approx_eq` helper is no longer used after the deletions but is harmless if dead-code warnings are not -D-warnings; if they are, also delete the `approx_eq` function (it was only used by the removed extension tests).

Verify dead-code state: after the deletions, search for `approx_eq` in the file. If it has no callers, delete the function.

- [ ] **Step 2: Run tests to verify the affected ones fail**

Run: `cargo test -p ppvm-stim --test normalize`
Expected: FAIL — `parse_extended` is in scope (via `pub use stim_parser::prelude::*`) but `normalize::to_tableau` still expects `&Program`, not `&ExtendedProgram`. Compilation error in the test file.

- [ ] **Step 3: Rewrite `crates/ppvm-stim/src/normalize.rs`**

Replace the entire file with:

```rust
use stim_parser::ast::AnnotationKind;
use stim_parser::extended::{Axis, ExtendedInstruction, ExtendedProgram};

use crate::tableau_program::{GateKind, Instruction, MeasureKind, NoiseKind, TableauProgram};

#[derive(Debug, thiserror::Error, PartialEq)]
#[non_exhaustive]
pub enum NormalizeError {
    #[error("unsupported instruction '{name}' at line {line} (phase 1)")]
    Unsupported { name: String, line: usize },

    #[error("'MPAD' at line {line} target #{index} = {value}, must be 0 or 1")]
    InvalidMPadTarget {
        line: usize,
        index: usize,
        value: usize,
    },
}

pub fn to_tableau(program: &ExtendedProgram) -> Result<TableauProgram, NormalizeError> {
    let mut out = Vec::with_capacity(program.instructions.len());
    let mut count = 0usize;
    normalize_slice(&program.instructions, &mut out, &mut count, 1)?;
    Ok(TableauProgram {
        instructions: out,
        expected_measurement_count: count,
    })
}

fn normalize_slice(
    src: &[ExtendedInstruction],
    out: &mut Vec<Instruction>,
    measure_count: &mut usize,
    enclosing_repeat_factor: u64,
) -> Result<(), NormalizeError> {
    use stim_parser::ast::{GateName, MeasureName, NoiseName};

    for ext in src {
        match ext {
            ExtendedInstruction::Gate { name, targets, line, .. } => {
                let kind = vanilla_gate_to_kind(*name, *line)?;
                out.push(Instruction::Gate {
                    kind,
                    targets: targets.clone(),
                    line: *line,
                });
            }
            ExtendedInstruction::Noise { name, args, targets, line, .. } => {
                let kind = vanilla_noise_to_kind(*name, *line)?;
                out.push(Instruction::Noise {
                    kind,
                    targets: targets.clone(),
                    args: args.clone(),
                    line: *line,
                });
            }
            ExtendedInstruction::Measure { name, args, targets, line, .. } => {
                let kind = match name {
                    MeasureName::M | MeasureName::MZ => MeasureKind::M,
                    MeasureName::MR => MeasureKind::MR,
                    other => {
                        return Err(NormalizeError::Unsupported {
                            name: other.canonical_name().to_string(),
                            line: *line,
                        });
                    }
                };
                let noise = args.first().copied().unwrap_or(0.0);
                *measure_count = measure_count.saturating_add(
                    targets
                        .len()
                        .saturating_mul(enclosing_repeat_factor as usize),
                );
                out.push(Instruction::Measure {
                    kind,
                    targets: targets.clone(),
                    noise,
                    line: *line,
                });
            }
            ExtendedInstruction::Annotation { line, .. } => {
                out.push(Instruction::Annotation { line: *line });
            }
            ExtendedInstruction::MPad { prob, bits, line, .. } => {
                let mut converted = Vec::with_capacity(bits.len());
                for (index, value) in bits.iter().copied().enumerate() {
                    let bit = match value {
                        0 => false,
                        1 => true,
                        _ => {
                            return Err(NormalizeError::InvalidMPadTarget {
                                line: *line,
                                index,
                                value,
                            });
                        }
                    };
                    converted.push(bit);
                }
                *measure_count = measure_count.saturating_add(
                    converted
                        .len()
                        .saturating_mul(enclosing_repeat_factor as usize),
                );
                out.push(Instruction::MPad {
                    bits: converted,
                    noise: prob.unwrap_or(0.0),
                    line: *line,
                });
            }
            ExtendedInstruction::Repeat { count, body, line } => {
                let mut inner = Vec::with_capacity(body.len());
                normalize_slice(
                    body,
                    &mut inner,
                    measure_count,
                    enclosing_repeat_factor.saturating_mul(*count),
                )?;
                out.push(Instruction::Repeat {
                    count: *count,
                    body: inner,
                    line: *line,
                });
            }

            // Promoted extension variants → direct mapping.
            ExtendedInstruction::T { targets, line } => {
                out.push(Instruction::Gate {
                    kind: GateKind::T,
                    targets: targets.clone(),
                    line: *line,
                });
            }
            ExtendedInstruction::TDag { targets, line } => {
                out.push(Instruction::Gate {
                    kind: GateKind::TDag,
                    targets: targets.clone(),
                    line: *line,
                });
            }
            ExtendedInstruction::Rotation { axis, theta, targets, line } => {
                let kind = match axis {
                    Axis::X => GateKind::RX { theta: *theta },
                    Axis::Y => GateKind::RY { theta: *theta },
                    Axis::Z => GateKind::RZ { theta: *theta },
                };
                out.push(Instruction::Gate {
                    kind,
                    targets: targets.clone(),
                    line: *line,
                });
            }
            ExtendedInstruction::U3 { theta, phi, lambda, targets, line } => {
                out.push(Instruction::Gate {
                    kind: GateKind::U3 {
                        theta: *theta,
                        phi: *phi,
                        lambda: *lambda,
                    },
                    targets: targets.clone(),
                    line: *line,
                });
            }
            ExtendedInstruction::Loss { p, targets, line } => {
                out.push(Instruction::Noise {
                    kind: NoiseKind::Loss,
                    targets: targets.clone(),
                    args: vec![*p],
                    line: *line,
                });
            }
            ExtendedInstruction::CorrelatedLoss { ps, targets, line } => {
                out.push(Instruction::Noise {
                    kind: NoiseKind::CorrelatedLoss,
                    targets: targets.clone(),
                    args: ps.to_vec(),
                    line: *line,
                });
            }
        }
    }
    Ok(())
}

fn vanilla_gate_to_kind(
    name: stim_parser::ast::GateName,
    line: usize,
) -> Result<GateKind, NormalizeError> {
    use stim_parser::ast::GateName::*;
    Ok(match name {
        Reset | ResetZ => GateKind::Reset,
        X => GateKind::X,
        Y => GateKind::Y,
        Z => GateKind::Z,
        H | HXZ => GateKind::H,
        S | SqrtZ => GateKind::S,
        SDag | SqrtZDag => GateKind::SDag,
        SqrtX => GateKind::SqrtX,
        SqrtXDag => GateKind::SqrtXDag,
        SqrtY => GateKind::SqrtY,
        SqrtYDag => GateKind::SqrtYDag,
        Identity => GateKind::I,
        CX | ZCX | CNot => GateKind::CX,
        CY | ZCY => GateKind::CY,
        CZ | ZCZ => GateKind::CZ,
        Swap | ISwap | ISwapDag | SqrtXX | SqrtYY | SqrtZZ | CXSwap | SwapCX | XCX | XCY | XCZ
        | YCX | YCY | YCZ | CXYZ | CZYX | HXY | HYZ => {
            return Err(NormalizeError::Unsupported {
                name: name.canonical_name().to_string(),
                line,
            });
        }
    })
}

fn vanilla_noise_to_kind(
    name: stim_parser::ast::NoiseName,
    line: usize,
) -> Result<NoiseKind, NormalizeError> {
    use stim_parser::ast::NoiseName::*;
    Ok(match name {
        Depolarize1 => NoiseKind::Depolarize1,
        Depolarize2 => NoiseKind::Depolarize2,
        PauliChannel1 => NoiseKind::PauliChannel1,
        PauliChannel2 => NoiseKind::PauliChannel2,
        XError => NoiseKind::XError,
        YError => NoiseKind::YError,
        ZError => NoiseKind::ZError,
        // IError reaching here means the parse_extended pass passed through
        // an unrecognized I_ERROR — but we made I_ERROR strict, so this is
        // unreachable in practice. Treat as Unsupported defensively.
        IError | HeraldedErase | HeraldedPauliChannel1 | CorrelatedError | ElseCorrelatedError => {
            return Err(NormalizeError::Unsupported {
                name: name.canonical_name().to_string(),
                line,
            });
        }
    })
}
```

Note: this version drops `find_tag`, `require_no_params`, `identity_to_kind`, `gate_to_kind`'s tag-handling, and `noise_to_kind`'s tag-handling — they all moved to `stim-parser::extended::interpret`. `NormalizeError::InvalidTag` is gone too.

The unused `_args` binding in the `Gate` arm is intentional: vanilla Stim gates have no args, so we drop them. (Today's normalizer also ignores them.)

- [ ] **Step 4: Update `crates/ppvm-stim/src/lib.rs`**

Replace the `Error` enum (currently lines 46-65, beginning with `#[derive(Debug, thiserror::Error)]`) with:

```rust
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    #[error(transparent)]
    Parse(#[from] ExtendedParseError),
    #[error(transparent)]
    Normalize(#[from] NormalizeError),
    #[error(transparent)]
    Exec(#[from] ExecError),
    #[error("failed to read stim file {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
}
```

Replace the body of `run_string` (only the `let prog = parse(src)?;` line changes; `let tprog = normalize::to_tableau(&prog)?;` is unchanged at the source level — `&prog` is now `&ExtendedProgram`):

```rust
    let prog = parse_extended(src)?;
    let tprog = normalize::to_tableau(&prog)?;
    let results = execute(&tprog, tab)?;
    Ok(results)
```

Update the module-level doc-test (lines 1-30 of the file). Replace lines 5-8 of the doc:

```rust
//! 1. [`parse_extended`] — `&str` → [`ExtendedProgram`] (re-exported from [`stim_parser`]).
//! 2. [`normalize::to_tableau`] — [`ExtendedProgram`] → [`TableauProgram`] (dialect resolved,
//!    unsupported instructions rejected).
//! 3. [`execute`] / [`sample`] — apply a [`TableauProgram`] to a [`GeneralizedTableau`].
```

Replace line 10:

```rust
//! Multi-shot usage should call [`parse_extended`] and [`normalize::to_tableau`] once and
```

Replace the doc-test code block (lines 17-27):

```rust
//! ```ignore
//! use ppvm_stim::{parse_extended, normalize, sample};
//! use ppvm_tableau::prelude::*;
//!
//! let prog = parse_extended(circuit_src)?;
//! let tprog = normalize::to_tableau(&prog)?;
//! let shots = sample(&tprog, 10_000, || {
//!     GeneralizedTableau::<_, usize, _>::new(n_qubits, 1e-10)
//! })?;
//! # Ok::<(), ppvm_stim::Error>(())
//! ```
```

- [ ] **Step 5: Update test imports — `executor.rs`, `stim_corpus.rs`**

In `crates/ppvm-stim/tests/executor.rs:1-3`, change:

```rust
use ppvm_runtime::config::indexmap::ByteFxHashF64;
use ppvm_stim::{execute, normalize, parse};
use ppvm_tableau::prelude::*;
```

to:

```rust
use ppvm_runtime::config::indexmap::ByteFxHashF64;
use ppvm_stim::{execute, normalize, parse_extended};
use ppvm_tableau::prelude::*;
```

Then in the file body, replace every `parse(` call with `parse_extended(`. Use `cargo check -p ppvm-stim --tests` afterward to confirm no stragglers.

In `crates/ppvm-stim/tests/stim_corpus.rs:1-5`, identical change: replace `parse` with `parse_extended` on line 4 and in every body call site.

- [ ] **Step 6: Update `crates/ppvm-stim/tests/run.rs` imports and `Error::Parse` matching**

Replace lines 1-3:

```rust
use ppvm_runtime::config::indexmap::ByteFxHashF64;
use ppvm_stim::{Error, ExtendedParseError, NormalizeError, ParseError, run_file, run_string};
use ppvm_tableau::prelude::*;
```

The existing test `run_string_propagates_parse_error` (lines 29-37) matches on `Error::Parse(ParseError::UnknownInstruction { .. })`. Update to wrap through `ExtendedParseError`:

```rust
#[test]
fn run_string_propagates_parse_error() {
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    let err = run_string("FROBNICATE 0", &mut tab).unwrap_err();
    assert!(matches!(
        err,
        Error::Parse(ExtendedParseError::Parse(ParseError::UnknownInstruction { .. }))
    ));
}
```

Append a new test at the end of the file that exercises an extension-specific parse error:

```rust
#[test]
fn run_string_propagates_extended_parse_error() {
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-10);
    let err = run_string("I[FOO] 0\n", &mut tab).unwrap_err();
    assert!(matches!(
        err,
        Error::Parse(ExtendedParseError::InvalidTag { .. })
    ));
}
```

- [ ] **Step 7: Update `crates/ppvm-stim/benches/tableau-msd-stim.rs`**

Change line 5 from:

```rust
use ppvm_stim::{execute, normalize, parse};
```

to:

```rust
use ppvm_stim::{execute, normalize, parse_extended};
```

Replace every `parse(...)` call in the file body with `parse_extended(...)`. The `TableauProgram` type at the function-signature level is unchanged.

- [ ] **Step 8: Run the ppvm-stim test suite**

Run: `cargo test -p ppvm-stim`
Expected: PASS — every test in `normalize.rs`, `executor.rs`, `run.rs`, `stim_corpus.rs` passes.

Run: `cargo build -p ppvm-stim --benches`
Expected: clean — the `tableau-msd-stim` bench compiles.

Note: `cargo check --workspace` will FAIL at this point because `ppvm-python-native` still calls the old API (`parse`, `ppvm_stim::ParseError`). That is fixed in Task 9. Verify `ppvm-stim` itself with `cargo check -p ppvm-stim` if you want to spot-check.

- [ ] **Step 9: Format and commit**

```bash
cargo fmt -p ppvm-stim
git add crates/ppvm-stim/src/normalize.rs crates/ppvm-stim/src/lib.rs crates/ppvm-stim/tests/normalize.rs crates/ppvm-stim/tests/executor.rs crates/ppvm-stim/tests/run.rs crates/ppvm-stim/tests/stim_corpus.rs crates/ppvm-stim/benches/tableau-msd-stim.rs
git commit -m "$(cat <<'EOF'
Migrate ppvm-stim to consume ExtendedProgram

Atomic refactor that switches to_tableau's input from &Program to
&ExtendedProgram and threads parse_extended / ExtendedParseError
through every consumer in the crate (run_string / run_file / Error,
executor / run / stim_corpus tests, the tableau-msd-stim bench).

Replaces tag-recognition arms in normalize.rs with direct mappings on
the typed extension variants (T, TDag, Rotation, U3, Loss,
CorrelatedLoss). Deletes find_tag, require_no_params, identity_to_kind,
and NormalizeError::InvalidTag — recognition lives in
stim-parser::extended now.

Removes the extension-recognition tests from tests/normalize.rs;
equivalent coverage lives in crates/stim-parser/tests/extended.rs.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Update `ppvm-python-native::stim_program`

Switch the Python-facing wrapper to call `parse_extended` and update its error-conversion helper. The `PyStimProgram.inner` is still a `TableauProgram`, so the Python API surface is unchanged.

**Files:**
- Modify: `crates/ppvm-python-native/src/stim_program.rs`

- [ ] **Step 1: Read the current file to confirm shape**

Read `crates/ppvm-python-native/src/stim_program.rs` (46 lines). Confirm the parse path in `parse()` (lines 16-20) calls `parse(src)` then `normalize::to_tableau(&prog)`, and that `stim_to_pyerr` (lines 39-41) takes `ppvm_stim::ParseError`.

- [ ] **Step 2: Write a parser test (compile-only)**

The Python-side test surface is in `ppvm-python/tests`, but a compile-time check that the new wrapper builds is enough. We rely on `cargo build --workspace` for verification.

- [ ] **Step 3: Replace the file**

Replace `crates/ppvm-python-native/src/stim_program.rs` with:

```rust
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;

use ppvm_stim::{TableauProgram, normalize, parse_extended};

/// Python-facing wrapper around `ppvm_stim::TableauProgram`.
#[pyclass(name = "StimProgram", module = "ppvm_python_native")]
pub struct PyStimProgram {
    pub(crate) inner: TableauProgram,
}

#[pymethods]
impl PyStimProgram {
    /// Parse + normalize a Stim circuit string.
    #[staticmethod]
    pub fn parse(src: &str) -> PyResult<Self> {
        let prog = parse_extended(src).map_err(stim_to_pyerr)?;
        let tprog = normalize::to_tableau(&prog).map_err(stim_to_pyerr_norm)?;
        Ok(Self { inner: tprog })
    }

    /// Read a `.stim` file and parse it.
    #[staticmethod]
    pub fn from_file(path: &str) -> PyResult<Self> {
        let src = std::fs::read_to_string(path)
            .map_err(|e| PyIOError::new_err(format!("failed to read {path}: {e}")))?;
        Self::parse(&src)
    }

    fn __repr__(&self) -> String {
        format!(
            "<StimProgram instructions={} measurements={}>",
            self.inner.instructions.len(),
            self.inner.expected_measurement_count
        )
    }
}

fn stim_to_pyerr(e: ppvm_stim::ExtendedParseError) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}

fn stim_to_pyerr_norm(e: ppvm_stim::NormalizeError) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}
```

Two changes: `parse` → `parse_extended` (line 4 import + line 17 call), and `ppvm_stim::ParseError` → `ppvm_stim::ExtendedParseError` in `stim_to_pyerr` (line 39).

- [ ] **Step 4: Run workspace check**

Run: `cargo build --workspace`
Expected: clean — `ppvm-python-native` builds against the new `ppvm-stim` API.

Run: `cargo test --workspace`
Expected: PASS — every Rust test in the workspace passes.

- [ ] **Step 5: Format and commit**

```bash
cargo fmt -p ppvm-python-native
git add crates/ppvm-python-native/src/stim_program.rs
git commit -m "$(cat <<'EOF'
Switch ppvm-python-native to parse_extended

PyStimProgram.parse now calls parse_extended and threads
ExtendedParseError through stim_to_pyerr. The Python-facing API is
unchanged (PyStimProgram.inner remains a TableauProgram).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: Workspace fmt + final test sweep

Final cleanup pass: run `cargo fmt --all`, `cargo test --workspace`, and `cargo doc --workspace --no-deps` to confirm intra-doc links still resolve.

**Files:**
- (No file changes if previous tasks formatted as they went; this task is a verification pass.)

- [ ] **Step 1: Run workspace fmt**

Run: `cargo fmt --all`
Expected: no diff — every previous task ran `cargo fmt -p <crate>` before committing.

If `cargo fmt --all` produces a diff, inspect the affected files. If they are part of an earlier task's scope, that's a signal the prior fmt step was missed; stage and commit:

```bash
git add -u
git commit -m "$(cat <<'EOF'
cargo fmt --all sweep

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

If the diff is empty, this commit step is skipped.

- [ ] **Step 2: Run full workspace test**

Run: `cargo test --workspace`
Expected: PASS — every test across `stim-parser`, `ppvm-stim`, `ppvm-tableau`, `ppvm-runtime`, `ppvm-sym`, `ppvm-python-native` passes.

- [ ] **Step 3: Run rustdoc with warnings as errors**

Run: `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps`
Expected: clean — every intra-doc link resolves. (`ppvm-stim/src/lib.rs` references `[parse_extended]`, `[ExtendedProgram]`, etc.; this verifies they all resolve through the prelude re-exports.)

- [ ] **Step 4: Done**

The feature is complete. Summary of the final state:
- `stim-parser::extended` owns recognition of all PPVM tag-based extensions.
- `ppvm-stim::normalize` is a near-1:1 translation from `ExtendedInstruction` to `TableauProgram::Instruction`.
- `ppvm-stim::run_string` / `run_file` and `ppvm-python-native::PyStimProgram::parse` call `parse_extended`.
- `NormalizeError::InvalidTag` is gone; tag errors fire at parse time as `ExtendedParseError::InvalidTag`.

---

## Spec coverage check

| Spec section | Task(s) |
|---|---|
| `extended/ast.rs` types | Task 1, Step 3 |
| `extended/parser.rs` (`ExtendedParseError`, `parse_extended`) | Task 1, Step 4 |
| `extended/interpret.rs` vanilla pass-through + Repeat recursion | Task 1, Step 5 |
| `S[T]` / `S_DAG[T]` recognition | Task 2 |
| `I[R_X]` / `I[R_Y]` / `I[R_Z]` recognition | Task 3 |
| `I[U3]` recognition | Task 4 |
| `I_ERROR[loss]` recognition | Task 5 |
| `I_ERROR[correlated_loss]` recognition (1-arg + 3-arg) | Task 6 |
| `mod.rs` re-exports + `lib.rs` prelude | Task 1, Step 6 + Task 7 |
| Tests in `stim-parser/tests/extended.rs` | Tasks 1–7 (incremental) |
| Move extension-recognition tests out of `ppvm-stim/tests/normalize.rs` | Task 8, Step 1 |
| Migrate `ppvm-stim/src/normalize.rs` to `&ExtendedProgram` | Task 8, Step 3 |
| Delete dead helpers (`find_tag`, `require_no_params`, `identity_to_kind`) | Task 8, Step 3 (full file rewrite drops them) |
| Drop `NormalizeError::InvalidTag` | Task 8, Step 3 |
| Update `ppvm-stim::run_string` / `run_file` / `Error` | Task 8, Step 4 |
| Update doc-test in `ppvm-stim/src/lib.rs` | Task 8, Step 4 |
| Update `ppvm-stim` test/bench `parse` → `parse_extended` | Task 8, Steps 5–7 |
| Update `ppvm-python-native/src/stim_program.rs` | Task 9 |
| `cargo fmt --all` + workspace test sweep | Task 10 |

All sections covered.
