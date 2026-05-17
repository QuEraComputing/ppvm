---
title: Stim Extended Dialect — Tag-Based Extensions Design
date: 2026-05-05
branch: david/ppvm-stim-2
---

# Stim Extended Dialect — Tag-Based Extensions Design

## Summary

PPVM's tag-based extensions to the Stim text format (`S[T]`, `I[R_X(theta=…)]`, `I[R_Y(...)]`, `I[R_Z(...)]`, `I[U3(theta=…, phi=…, lambda=…)]`, `I_ERROR[loss]`, `I_ERROR[correlated_loss]`) are currently recognized inside `ppvm-stim/src/normalize.rs`. Per [issue #44](https://github.com/QuEraComputing/ppvm/issues/44), they should be "handled as grammar extensions in the chumsky parser" so that any consumer of stim-parser (FLAIR pipeline, Pauli Propagation runtime, future bytecode lowering) gets typed extension nodes without re-implementing the recognition.

This spec moves recognition into a new `stim-parser::extended` module that produces a typed `ExtendedProgram` from the existing `Program`, and migrates `ppvm-stim` to consume it. Extensions become first-class `ExtendedInstruction` variants (Flavor B promotion). Recognition is implemented as a post-pass over `Program` rather than by extending the chumsky grammar directly — the post-pass approach has equivalent observable behavior with substantially less combinator complexity (see "Decisions" #2).

## Context

### Current state

`stim-parser` produces a vanilla `Program` with `RawInstruction` variants where tags are preserved as opaque `Tag { name, params }` blobs. `ppvm-stim/src/normalize.rs` then dispatches on host-instruction name and pattern-matches on tags to translate extensions into `TableauProgram::Instruction` variants. Concretely:

- `gate_to_kind` ([crates/ppvm-stim/src/normalize.rs:150-195](../../crates/ppvm-stim/src/normalize.rs#L150-L195)) handles `S[T]`/`S_DAG[T]` and delegates to `identity_to_kind` for `I[R_X/R_Y/R_Z/U3]`.
- `identity_to_kind` ([crates/ppvm-stim/src/normalize.rs:197-242](../../crates/ppvm-stim/src/normalize.rs#L197-L242)) does the named-param lookup for rotation and U3 tags.
- `noise_to_kind` ([crates/ppvm-stim/src/normalize.rs:244-318](../../crates/ppvm-stim/src/normalize.rs#L244-L318)) handles `I_ERROR[loss]`/`I_ERROR[correlated_loss]`.
- Helpers `find_tag`, `require_no_params` ([normalize.rs:334-348](../../crates/ppvm-stim/src/normalize.rs#L334-L348)) support all of the above.
- `NormalizeError::InvalidTag` is the error variant for malformed tag shapes.

The parser is unaware of these extensions; an `R_X` tag missing its `theta` named param fails as a `NormalizeError` after parsing succeeds, with a line number but not a tag-span.

### Why move it

1. **Multi-consumer reuse.** Issue #44 names FLAIR, Pauli Propagation runtime, and bytecode lowering as consumers. Today only `ppvm-stim` recognizes the extensions; any other consumer would re-derive the same matching.
2. **Errors closer to source.** `R_X` without `theta` belongs to the parser layer's "this isn't a well-formed extension" stance, not to the normalizer's "this can't run on a tableau" stance.
3. **Issue #44 explicitly calls for this** in its closing line: "Our custom extensions like `I[...]` and `I_ERROR[...]` … will be handled as grammar extensions in the chumsky parser."

### What "extended dialect" covers

| Host instruction | Tag(s) | Extension |
|---|---|---|
| `S` | `[T]` (no params) | T gate |
| `S_DAG` | `[T]` (no params) | T† gate |
| `I` | `[R_X(theta=<f64>)]` | RX rotation |
| `I` | `[R_Y(theta=<f64>)]` | RY rotation |
| `I` | `[R_Z(theta=<f64>)]` | RZ rotation |
| `I` | `[U3(theta=<f64>, phi=<f64>, lambda=<f64>)]` | U3 rotation |
| `I_ERROR` | `[loss]`, 1 instruction-arg | Single-qubit loss |
| `I_ERROR` | `[correlated_loss]`, 1 or 3 instruction-args | Correlated loss |

These four host instructions (`S`, `S_DAG`, `I`, `I_ERROR`) are the *strict* hosts — any tag on them that isn't in the recognized set is an error. All other instructions retain their tags verbatim and are passed through unchanged.

## Decisions

| # | Decision |
|---|---|
| 1 | Extensions are promoted to first-class `ExtendedInstruction` variants (Flavor B), not retained as typed tags. ppvm-stim's normalize matches on the typed variants directly, eliminating the per-tag recognition logic. |
| 2 | Recognition is a post-pass over `Program` rather than chumsky grammar combinators. Same observable behavior; simpler implementation; no backtracking complexity; named-param tags (U3 in particular) require post-validation either way, so option 1's purity is partially fictional. |
| 3 | Strict for known hosts (`S`, `S_DAG`, `I`, `I_ERROR`); lenient elsewhere. Catches typos like `I[R_x]` or `I_ERROR[Loss]` at parse time, where the user expects diagnostic precision, while leaving room for unrelated tags on other instructions to pass through untouched. |
| 4 | Vanilla pass-through `ExtendedInstruction` variants retain `tags: Vec<Tag>` on hosts where tags are lenient. MPad and Measure carry `tags` opaquely (no extensions defined for them today). |
| 5 | `parse_extended` and `ExtendedProgram` go in `stim_parser::prelude::*` alongside `parse` and `Program`. Power-user separation isn't worth the friction; ppvm-stim and other consumers grab everything from one prelude. |
| 6 | The migration of `ppvm-stim/src/normalize.rs` to consume `ExtendedProgram` happens in this same plan. The boundary between the new interpret pass and the simplified normalize is what makes the design coherent; bundling means dead code is deleted in the same change that justifies removing it. |
| 7 | `NormalizeError::InvalidTag` is deleted. The variant becomes unreachable once recognition moves upstream. |
| 8 | `Tag` and `TagParam` AST types do **not** gain spans in this round. Error messages use the host instruction's `line` number (already on every `RawInstruction`). Span-precision improvements are out of scope. |

## Architecture

### Module layout

```
crates/stim-parser/src/
  ast.rs               (unchanged)
  grammar.rs           (unchanged)
  line_map.rs          (unchanged)
  parser.rs            (unchanged — vanilla parse)
  table.rs             (unchanged)
  lib.rs               (adds `pub mod extended` and prelude exports)
  extended/
    mod.rs             — declarations + re-exports
    ast.rs             — ExtendedProgram, ExtendedInstruction, Axis
    parser.rs          — parse_extended() entry point + ExtendedParseError
    interpret.rs       — Program → ExtendedProgram post-pass
```

### Pipeline

```
src ──parse──▶ Program ──interpret──▶ ExtendedProgram
                                         │
                                         ▼
                            ppvm-stim::normalize::to_tableau
                                         │
                                         ▼
                                 TableauProgram
                                         │
                                         ▼
                                 ppvm-stim::execute
```

`parse_extended` orchestrates `parse` + `interpret`. ppvm-stim's `normalize::to_tableau` switches input from `&Program` to `&ExtendedProgram` and becomes a near-1:1 translation pass.

### `ExtendedInstruction` AST

```rust
// stim-parser/src/extended/ast.rs

pub struct ExtendedProgram {
    pub instructions: Vec<ExtendedInstruction>,
}

pub enum ExtendedInstruction {
    // Vanilla pass-through. Tags retained for unrecognized hosts.
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

pub enum Axis { X, Y, Z }
```

`#[non_exhaustive]` on both `ExtendedInstruction` and `Axis` to allow future extensions without breaking changes.

### Strictness rules

For each strict host, the recognized tag shapes are exhaustive:

- **`S`**: `[T]` (no params, no instruction-args) → `T`. Any other tag → `InvalidTag`. No tag → vanilla `Gate { name: S, … }`.
- **`S_DAG`**: `[T]` (no params, no instruction-args) → `TDag`. Any other tag → `InvalidTag`. No tag → vanilla `Gate { name: SDag, … }`.
- **`I`**: with exactly one tag, recognized tags are `R_X(theta=…)`, `R_Y(theta=…)`, `R_Z(theta=…)`, `U3(theta=…, phi=…, lambda=…)`. Any other tag, or more than one tag, or unrecognized tag name → `InvalidTag`. No tag → vanilla `Gate { name: Identity, … }`.
- **`I_ERROR`**: with exactly one tag, recognized tags are `loss` (no params, requires 1 instruction-arg) and `correlated_loss` (no params, requires 1 or 3 instruction-args). Any other tag, more than one tag, or wrong arg count → `InvalidTag`. No tag → `InvalidTag` (an `I_ERROR` without a recognized tag is malformed).

Lenient hosts (everything else): pass through with `tags` retained on the vanilla variant.

### Errors

```rust
// stim-parser/src/extended/parser.rs

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
```

The `InvalidTag` variant has the same shape as today's `NormalizeError::InvalidTag` to ease migration. Error messages remain prose explanations; no machine-readable error codes.

### `parse_extended` entry point

```rust
// stim-parser/src/extended/parser.rs

pub fn parse_extended(src: &str) -> Result<ExtendedProgram, ExtendedParseError> {
    let prog = crate::parse(src)?;
    Ok(crate::extended::interpret::interpret(&prog)?)
}
```

`interpret` is `pub(crate)` and takes `&Program`, returns `Result<ExtendedProgram, ExtendedParseError>`. Recursion into `Repeat` bodies is straightforward; the post-pass does not preserve any state across instructions.

### Public API

`stim_parser::prelude` exports:

```rust
pub mod prelude {
    pub use crate::ast::{
        AnnotationKind, GateName, MeasureName, NoiseName, ParseError, Program,
        RawInstruction, Tag, TagParam,
    };
    pub use crate::parser::parse;
    pub use crate::extended::{
        Axis, ExtendedInstruction, ExtendedParseError, ExtendedProgram, parse_extended,
    };
}
```

`stim_parser::extended::*` is also accessible directly for consumers that want to be explicit.

## ppvm-stim migration

`ppvm-stim/src/normalize.rs` switches input from `&Program` to `&ExtendedProgram`:

- `normalize_slice` matches on `&[ExtendedInstruction]` instead of `&[RawInstruction]`. Each typed variant maps directly to a `TableauProgram::Instruction`:
  - `ExtendedInstruction::T { targets, line }` → `Instruction::Gate { kind: GateKind::T, … }`
  - `ExtendedInstruction::TDag { … }` → `GateKind::TDag`
  - `ExtendedInstruction::Rotation { axis: Axis::X, theta, … }` → `GateKind::RX { theta }`, similarly for Y/Z
  - `ExtendedInstruction::U3 { theta, phi, lambda, … }` → `GateKind::U3 { theta, phi, lambda }`
  - `ExtendedInstruction::Loss { p, … }` → `Instruction::Noise { kind: NoiseKind::Loss, args: vec![p], … }`
  - `ExtendedInstruction::CorrelatedLoss { ps, … }` → `Instruction::Noise { kind: NoiseKind::CorrelatedLoss, args: ps.to_vec(), … }`
- `gate_to_kind` shrinks: only the standard Stim Clifford gates remain. `S[T]` and `S_DAG[T]` arms are gone.
- `identity_to_kind` is deleted entirely.
- `noise_to_kind` shrinks: `I_ERROR` arms are gone. The remaining noise channels are direct mappings.
- Helpers `find_tag`, `require_no_params` are deleted (no callers).
- `NormalizeError::InvalidTag` variant is deleted. Remaining variants: `Unsupported`, `InvalidMPadTarget`.
- `to_tableau` signature: `pub fn to_tableau(program: &ExtendedProgram) -> Result<TableauProgram, NormalizeError>`.

`ppvm-stim/src/lib.rs`:

- `run_string` calls `parse_extended` instead of `parse`. Its `Error` enum's `Parse` variant holds `ExtendedParseError` instead of `ParseError`. The change is local; consumers using `run_string` get a wider error type.
- The doc-test example switches `parse` + `normalize::to_tableau` to `parse_extended` + `normalize::to_tableau` (the second call now takes `&ExtendedProgram`).

`ppvm-stim/src/lib.rs` re-exports: existing `pub use stim_parser::prelude::*;` already brings in the new types after the prelude update — no manual re-export changes needed.

`ppvm-python-native/src/stim_program.rs` uses `ppvm_stim::{TableauProgram, normalize, parse}` ([crates/ppvm-python-native/src/stim_program.rs:4](../../crates/ppvm-python-native/src/stim_program.rs#L4)). Switch to `parse_extended` and update the error-conversion helper. The Python-facing API does not change.

## Testing

### New tests in `stim-parser/tests/extended.rs`

For each recognized extension shape, one happy-path test and one strict-host error test:

| Test | Shape under test |
|---|---|
| `recognizes_s_t` | `S[T] 0 1` → `ExtendedInstruction::T` |
| `recognizes_s_dag_t` | `S_DAG[T] 0` → `ExtendedInstruction::TDag` |
| `recognizes_i_r_x` | `I[R_X(theta=0.5)] 0` → `ExtendedInstruction::Rotation { axis: X, theta: 0.5 }` |
| `recognizes_i_r_y` / `recognizes_i_r_z` | analogous |
| `recognizes_i_u3` | `I[U3(theta=0.1, phi=0.2, lambda=0.3)] 0` → `ExtendedInstruction::U3` |
| `recognizes_i_error_loss` | `I_ERROR[loss](0.01) 0` → `ExtendedInstruction::Loss` |
| `recognizes_i_error_correlated_loss_one_arg` | `I_ERROR[correlated_loss](0.01) 0 1` → `ExtendedInstruction::CorrelatedLoss { ps: [0.01, 0.0, 0.0] }` |
| `recognizes_i_error_correlated_loss_three_args` | `(0.01, 0.02, 0.03)` → `ps: [0.01, 0.02, 0.03]` |

Strict-host error tests:

| Test | Input | Expected |
|---|---|---|
| `s_with_unknown_tag_errors` | `S[X] 0` | `InvalidTag` |
| `i_with_no_tag_passes_through` | `I 0` | `Gate { name: Identity, tags: [], … }` (lenient when no tag) |
| `i_with_unknown_tag_errors` | `I[FOO] 0` | `InvalidTag` |
| `i_r_x_missing_theta_errors` | `I[R_X] 0` | `InvalidTag` |
| `i_u3_missing_phi_errors` | `I[U3(theta=0.1, lambda=0.2)] 0` | `InvalidTag` |
| `i_error_with_no_tag_errors` | `I_ERROR(0.1) 0` | `InvalidTag` |
| `i_error_loss_wrong_arg_count_errors` | `I_ERROR[loss](0.1, 0.2) 0` | `InvalidTag` |
| `i_error_correlated_loss_two_args_errors` | `I_ERROR[correlated_loss](0.1, 0.2) 0 1` | `InvalidTag` |

Lenient pass-through tests:

| Test | Input | Expected |
|---|---|---|
| `unknown_tag_on_h_passes_through` | `H[unrelated] 0` | `Gate { name: H, tags: [Tag { name: "unrelated", … }], … }` |
| `tag_on_m_passes_through` | `M[foo] 0` | `Measure { name: M, tags: [Tag { name: "foo", … }], … }` |
| `repeat_recurses_into_body` | `REPEAT 2 { I[R_X(theta=pi)] 0 }` | `Repeat { count: 2, body: [Rotation { … }] }` |

### Tests moving from `ppvm-stim/tests/normalize.rs` to `stim-parser/tests/extended.rs`

Tests in [crates/ppvm-stim/tests/normalize.rs](../../crates/ppvm-stim/tests/normalize.rs) that exercise extension-tag recognition (e.g., `s_t_promotes_to_t`, `i_r_x_normalizes`, `i_error_loss_normalizes`, etc.) move with their input/expected-output adjusted to the new types. Their *intent* remains the same: validate that the recognition produces the correct typed variant.

Tests that exercise normalize-specific concerns (MPAD bit validation, `Unsupported` rejection of `Swap`/`ISwap`/`HeraldedErase`, etc.) stay in `ppvm-stim/tests/normalize.rs`.

### Pre-existing executor and stim_corpus tests

[crates/ppvm-stim/tests/executor.rs](../../crates/ppvm-stim/tests/executor.rs) and [crates/ppvm-stim/tests/stim_corpus.rs](../../crates/ppvm-stim/tests/stim_corpus.rs) drive the full pipeline via `parse` + `normalize::to_tableau` + `execute`. They switch to `parse_extended`. Test bodies otherwise unchanged; output bits should be bit-identical.

## Out of scope

- Adding spans to `Tag` / `TagParam`. Error messages continue to use the host instruction's `line` number.
- Promoting any further Stim instructions to extension status. The eight extensions in the table above are the complete set.
- Lenient mode for known hosts (a flag-driven "accept unknown tags on `I`/`I_ERROR`/`S`/`S_DAG` and pass them through"). Strict is the only mode.
- A dialect-flag on `parse` itself. The two-layer `parse` + `parse_extended` design is final.
- Recognition for new tag shapes on `M`/`MZ`/`MR`, `MPad`, or annotations. These remain opaque-tag-bearing pass-through variants.
- Any change to `ppvm-vihaco`. That crate's `Circuit` continues to operate against `GeneralizedTableau` directly and is not affected.

## Implementation order

1. Add `extended/ast.rs` with `ExtendedProgram`, `ExtendedInstruction`, `Axis`. No interpretation yet — types only.
2. Add `extended/parser.rs` with `ExtendedParseError` (`Parse` and `InvalidTag` variants).
3. Implement `extended/interpret.rs` with vanilla pass-through (all non-extension instructions, including `Repeat` recursion). `parse_extended` lives in `extended/parser.rs` and calls `interpret`. Tests cover pass-through and lenient hosts.
4. Add recognition for `S[T]` / `S_DAG[T]`. Tests + strict error cases.
5. Add recognition for `I[R_X]` / `I[R_Y]` / `I[R_Z]`. Tests.
6. Add recognition for `I[U3]`. Tests.
7. Add recognition for `I_ERROR[loss]` / `I_ERROR[correlated_loss]`. Tests.
8. Wire `extended/mod.rs` re-exports; update `lib.rs` and `prelude`.
9. Move the extension-flavored tests from `ppvm-stim/tests/normalize.rs` to `stim-parser/tests/extended.rs`. Re-shape inputs/expected outputs to the new types.
10. Switch `ppvm-stim/src/normalize.rs` to consume `&ExtendedProgram`. Delete the dead helpers (`find_tag`, `require_no_params`, `identity_to_kind`) and the now-unreachable extension arms in `gate_to_kind` and `noise_to_kind`. Drop `NormalizeError::InvalidTag`.
11. Update `ppvm-stim/src/lib.rs` (`run_string`/`run_file`, `Error::Parse` variant) and the doc-test example. Update `ppvm-python-native/src/stim_program.rs` to call `parse_extended`.
12. `cargo fmt --all` workspace-wide. Run `cargo test --workspace` to confirm full test suite passes. Commit.

## Open questions

None. Decisions on MPad/Measure tag opacity, strictness scope, and migration bundling were resolved during design.
