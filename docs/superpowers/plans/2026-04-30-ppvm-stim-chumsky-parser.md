# ppvm-stim Chumsky Parser Rewrite Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the hand-written line-based Stim parser in `crates/ppvm-stim/src/parser/mod.rs` with a chumsky 0.12 grammar, while preserving the public API (`parse(&str) -> Result<Program, ParseError>`) and the existing `Program` AST.

**Architecture:** Two-stage internal pipeline. (1) chumsky grammar in `parser/grammar.rs` produces an internal `Vec<RawSyntaxNode>` (raw name strings, raw target lexemes, evaluated args, parsed tags, byte spans) and forwards chumsky's `Rich<'static, char>` for syntax errors. (2) A post-pass `validate_program` in `parser/mod.rs` walks the raw tree and emits typed `ParseError::{UnknownInstruction, ArgCount, TargetCount}` validation errors plus the final `Program` AST. `LineMap` (already present, currently dead-coded) is re-activated for line/col conversion.

**Tech Stack:** Rust 2024, chumsky 0.12 (already in `Cargo.toml`), thiserror 2.

**Spec:** `docs/superpowers/specs/2026-04-30-ppvm-stim-chumsky-parser-design.md`.

**Branch:** `david/ppvm-stim-2` (existing). All work continues on this branch.

**Standing rules for every commit:**
- Run `cargo fmt --all` *before* every `git commit`.
- Test corpus (`tests/data/`, `tests/regen-stim/`, `tests/stim_corpus.rs`) is **not** modified.

---

## File-by-file overview

| File | Action | Responsibility after rewrite |
|---|---|---|
| `crates/ppvm-stim/src/parser/mod.rs` | Heavily modify | `pub fn parse`, `pub struct LineMap`, internal raw types (`RawSyntaxNode`, `RawTag`, `RawTagParam`, `RawTarget`), `validate_program` post-pass. ~250 LOC. |
| `crates/ppvm-stim/src/parser/grammar.rs` | Create | chumsky 0.12 combinators, top-to-bottom: pad → numbers → ident → tags → args → targets → instruction → repeat → program. Single `pub(crate) fn program_parser()` export. ~300 LOC. |
| `crates/ppvm-stim/src/parser/ast.rs` | Modify | Add `SyntaxError` struct (`Rich<'static, char>` + `Arc<LineMap>` + `Display`). Change `ParseError::Syntax { line, col, message }` to `ParseError::Syntax(SyntaxError)`. Drop `PartialEq` on `ParseError`. Other AST items unchanged. |
| `crates/ppvm-stim/src/parser/table.rs` | No change | Source of truth for instruction names + arity rules; consumed by validator. |
| `crates/ppvm-stim/src/lib.rs` | No change | Re-exports already cover `ParseError`, `parse`, etc. |
| `crates/ppvm-stim/tests/parser_errors.rs` | Pattern updates | `ParseError::Syntax { .. }` → `ParseError::Syntax(_)`. Add new line/col-pinning test. |
| All other `tests/*.rs` | No change | Run as-is. |

---

## Reference: chumsky 0.12 API summary

These are the chumsky API shapes used throughout the plan. If any spelling diverges in the actual chumsky 0.12 release, adjust the call sites and keep the structure.

```rust
use chumsky::prelude::*;
use chumsky::error::Rich;
use chumsky::extra;
use chumsky::span::SimpleSpan;

// Parser type alias for our use case.
type Extra<'src> = extra::Err<Rich<'src, char>>;

// Constructing a custom error from a span and message:
let span: SimpleSpan<usize> = (byte_start..byte_end).into();
let rich = Rich::<char>::custom(span, "some message");

// Converting a borrowed Rich to 'static:
let owned: Rich<'static, char> = borrowed_rich.into_owned();

// Accessing a Rich's span:
let span: &SimpleSpan<usize> = rich.span();

// Recursive parser:
recursive(|body_ref| {
    // ... use body_ref where the parser refers to itself ...
})

// Span access in map:
parser.map_with(|out, e| (out, e.span()))
```

If `map_with` does not give a `MapExtra` with `.span()` in your chumsky 0.12 release (older chumsky used `map_with_span(|out, span| ...)`), use whichever spelling compiles. The goal is: **for each parsed instruction head, capture the byte span so the validator can derive a line number.**

---

### Task 1: Add `LineMap::line_col` and make `LineMap` `pub`

Activate the dead-coded `LineMap`, give it a `line_col` method, and expose it for use by `SyntaxError`.

**Files:**
- Modify: `crates/ppvm-stim/src/parser/mod.rs:97-122` (the existing `LineMap` block)

- [ ] **Step 1: Add a small unit test for `LineMap::line_col`**

Add this `#[cfg(test)] mod tests` block at the bottom of `crates/ppvm-stim/src/parser/mod.rs`. (Do not collide with existing private tests — append to file.)

```rust
#[cfg(test)]
mod line_map_tests {
    use super::LineMap;

    #[test]
    fn line_col_at_line_start() {
        let m = LineMap::new("abc\ndef\nghi");
        assert_eq!(m.line_col(0), (1, 1));
        assert_eq!(m.line_col(4), (2, 1));
        assert_eq!(m.line_col(8), (3, 1));
    }

    #[test]
    fn line_col_mid_line() {
        let m = LineMap::new("abc\ndef\nghi");
        assert_eq!(m.line_col(2), (1, 3));
        assert_eq!(m.line_col(6), (2, 3));
    }

    #[test]
    fn line_col_at_eof() {
        let m = LineMap::new("abc\ndef");
        assert_eq!(m.line_col(7), (2, 4));
    }
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test -p ppvm-stim line_map_tests -- --nocapture`
Expected: compile error — `line_col` does not exist.

- [ ] **Step 3: Implement `line_col` and make `LineMap` public**

Replace the entire existing `LineMap` block (lines 97–122 of `mod.rs`) with:

```rust
/// Maps byte offsets in source to 1-indexed line/column positions.
pub struct LineMap {
    /// `starts[i]` = byte offset of the start of line (i+1).
    starts: Vec<usize>,
}

impl std::fmt::Debug for LineMap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LineMap")
            .field("lines", &self.starts.len())
            .finish()
    }
}

impl LineMap {
    /// Build a `LineMap` for `src`.
    pub fn new(src: &str) -> Self {
        let mut starts = vec![0];
        for (i, b) in src.bytes().enumerate() {
            if b == b'\n' {
                starts.push(i + 1);
            }
        }
        Self { starts }
    }

    /// 1-indexed line number for a byte offset.
    pub fn line_of(&self, byte_offset: usize) -> usize {
        match self.starts.binary_search(&byte_offset) {
            Ok(i) => i + 1,
            Err(i) => i, // i is the insertion index; start of line `i` is at starts[i-1].
        }
    }

    /// 1-indexed `(line, col)` for a byte offset.
    pub fn line_col(&self, byte_offset: usize) -> (usize, usize) {
        let line = self.line_of(byte_offset);
        let line_start = self.starts[line - 1];
        let col = byte_offset - line_start + 1;
        (line, col)
    }
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cargo test -p ppvm-stim line_map_tests`
Expected: all three tests PASS.

- [ ] **Step 5: Run the rest of the crate's tests to confirm no regressions**

Run: `cargo test -p ppvm-stim`
Expected: all PASS (the existing parser still constructs `ParseError::Syntax { line, col, message }`; nothing else has changed).

- [ ] **Step 6: Format and commit**

```bash
cd /Users/david/git/claude-repos/ppvm
cargo fmt --all
git add crates/ppvm-stim/src/parser/mod.rs
git commit -m "$(cat <<'EOF'
Expose LineMap and add line_col helper

Activates the previously dead-coded LineMap (which was marked
allow(dead_code) for a planned chumsky integration) by making it
public and adding a 1-indexed (line, col) lookup that the upcoming
SyntaxError variant will use for error display.
EOF
)"
```

---

### Task 2: Define `SyntaxError` struct and change `ParseError::Syntax` to a tuple variant

The error-shape API change is isolated here. The hand-written parser keeps working — its `ParseError::Syntax { line, col, message }` construction sites are converted to use `SyntaxError::synth(line, col, message, &line_map)` which builds a `Rich::custom` from a synthetic span. After this task, the public API shape matches the spec; the new chumsky parser is built on this shape in later tasks.

**Files:**
- Modify: `crates/ppvm-stim/src/parser/ast.rs:187-215` (the `ParseError` enum)
- Add to: `crates/ppvm-stim/src/parser/ast.rs` (top-of-file imports + new `SyntaxError` struct + `Display` impl)
- Modify: `crates/ppvm-stim/src/parser/mod.rs` (every `ParseError::Syntax { ... }` construction site, plus thread `Arc<LineMap>` into the parser fns)
- Modify: `crates/ppvm-stim/tests/parser_errors.rs` (pattern updates)

- [ ] **Step 1: Update the failing test in `parser_errors.rs` and add a new line/col pin**

Replace the contents of `crates/ppvm-stim/tests/parser_errors.rs` with:

```rust
use ppvm_stim::{ParseError, parse};

#[test]
fn unknown_instruction_returns_error() {
    let err = parse("FROBNICATE 0").unwrap_err();
    match err {
        ParseError::UnknownInstruction { name, line } => {
            assert_eq!(name, "FROBNICATE");
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn arg_count_mismatch() {
    let err = parse("DEPOLARIZE1(0.1, 0.2) 0").unwrap_err();
    match err {
        ParseError::ArgCount {
            name,
            expected,
            found,
            line,
        } => {
            assert_eq!(name, "DEPOLARIZE1");
            assert_eq!(expected, 1);
            assert_eq!(found, 2);
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn target_pair_violation() {
    let err = parse("CX 0 1 2").unwrap_err();
    match err {
        ParseError::TargetCount {
            name,
            divisor,
            found,
            line,
        } => {
            assert_eq!(name, "CX");
            assert_eq!(divisor, 2);
            assert_eq!(found, 3);
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn at_least_one_target_required_for_h() {
    let err = parse("H").unwrap_err();
    assert!(matches!(
        err,
        ParseError::TargetCount { .. } | ParseError::Syntax(_)
    ));
}

#[test]
fn invalid_target_yields_syntax_error() {
    let err = parse("H abc").unwrap_err();
    assert!(matches!(err, ParseError::Syntax(_)));
}

#[test]
fn unclosed_bracket_yields_syntax_error() {
    let err = parse("S[T 0").unwrap_err();
    assert!(matches!(err, ParseError::Syntax(_)));
}

#[test]
fn line_numbers_in_errors_are_correct() {
    let err = parse("X 0\nY 0\nFROBNICATE 0").unwrap_err();
    match err {
        ParseError::UnknownInstruction { line, .. } => assert_eq!(line, 3),
        other => panic!("{other:?}"),
    }
}

#[test]
fn syntax_error_includes_line_and_col() {
    // Pin the new Display behavior: line and col both appear in the message.
    let err = parse("H 0\nH abc").unwrap_err();
    let s = err.to_string();
    assert!(s.contains("line 2"), "message was: {s}");
    assert!(s.contains("col"), "message was: {s}");
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test -p ppvm-stim --test parser_errors`
Expected: compile error — pattern `ParseError::Syntax(_)` does not match the current struct-style variant.

- [ ] **Step 3: Add imports + `SyntaxError` to `ast.rs`**

At the very top of `crates/ppvm-stim/src/parser/ast.rs`, add (or extend) imports:

```rust
//! Pure-Stim AST. Tags are preserved verbatim; dialect resolution
//! happens in `crate::normalize`.

use std::fmt;
use std::sync::Arc;

use chumsky::error::Rich;

use super::LineMap;
```

Then, immediately above the existing `ParseError` enum, add the `SyntaxError` struct + `Display` impl:

```rust
/// Carries a chumsky 0.12 `Rich<char>` error plus a shared `LineMap`
/// so that `Display` formats `line:col` consistently with the typed
/// validation variants of [`ParseError`].
#[derive(Debug, Clone)]
pub struct SyntaxError {
    pub(crate) rich: Rich<'static, char>,
    pub(crate) line_map: Arc<LineMap>,
}

impl SyntaxError {
    /// Synthesise a `SyntaxError` from a (line, col) position. Used by
    /// the validator when it needs to emit a syntax error for an issue
    /// the grammar could not catch (e.g. an annotation-tolerated target
    /// that fails `usize` parsing for a non-annotation instruction).
    pub(crate) fn synth(
        line: usize,
        col: usize,
        message: impl Into<String>,
        line_map: Arc<LineMap>,
    ) -> Self {
        let line_idx = line.saturating_sub(1);
        let line_start = line_map.starts_at(line_idx).unwrap_or(0);
        let byte = line_start + col.saturating_sub(1);
        let span = chumsky::span::SimpleSpan::from(byte..byte);
        let rich = Rich::<char>::custom(span, message.into());
        SyntaxError { rich, line_map }
    }

    /// Construct from a chumsky `Rich` (with any borrow lifetime) and
    /// a shared `LineMap`. `into_owned` widens the lifetime to `'static`.
    pub(crate) fn from_rich<'src>(
        rich: Rich<'src, char>,
        line_map: Arc<LineMap>,
    ) -> Self {
        SyntaxError {
            rich: rich.into_owned(),
            line_map,
        }
    }
}

impl fmt::Display for SyntaxError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let span = self.rich.span();
        let (line, col) = self.line_map.line_col(span.start);
        write!(f, "syntax error at line {line}, col {col}: {}", self.rich)
    }
}
```

`SyntaxError::synth` references `LineMap::starts_at`, which we add in the next step.

- [ ] **Step 4: Add `LineMap::starts_at`**

In `crates/ppvm-stim/src/parser/mod.rs`, add a method to `impl LineMap`:

```rust
impl LineMap {
    // ... existing methods ...

    /// Byte offset of the start of line `(line_idx + 1)`. `None` for
    /// out-of-range indices.
    pub fn starts_at(&self, line_idx: usize) -> Option<usize> {
        self.starts.get(line_idx).copied()
    }
}
```

- [ ] **Step 5: Replace `ParseError::Syntax` variant**

In `crates/ppvm-stim/src/parser/ast.rs`, replace the existing `ParseError` enum (lines 187–215 of the original file) with:

```rust
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ParseError {
    #[error("{0}")]
    Syntax(SyntaxError),

    #[error("unknown instruction '{name}' at line {line}")]
    UnknownInstruction { name: String, line: usize },

    #[error("'{name}' at line {line} expected {expected} args, got {found}")]
    ArgCount {
        name: String,
        expected: usize,
        found: usize,
        line: usize,
    },

    #[error("'{name}' at line {line} expected target count divisible by {divisor}, got {found}")]
    TargetCount {
        name: String,
        divisor: usize,
        found: usize,
        line: usize,
    },
}
```

Note: `PartialEq` is dropped (chumsky's `Rich` is not `PartialEq`). Existing tests use `matches!` rather than `==`; verify in Step 8.

- [ ] **Step 6: Thread `Arc<LineMap>` through the hand-written parser and convert every `Syntax` construction site**

In `crates/ppvm-stim/src/parser/mod.rs`:

(a) Update `pub fn parse`:

```rust
pub fn parse(src: &str) -> Result<Program, ParseError> {
    use std::sync::Arc;
    let line_map = Arc::new(LineMap::new(src));
    let tokens = tokenize_lines(src);
    let (instructions, rest) = parse_block(&tokens, &line_map, false)?;
    if !rest.is_empty() {
        return Err(syntax(
            tokens[tokens.len() - rest.len()].line,
            1,
            "unexpected '}' without matching REPEAT".to_string(),
            &line_map,
        ));
    }
    Ok(Program { instructions })
}
```

(b) Add a helper near the top of `mod.rs` (just below the imports block):

```rust
use std::sync::Arc;

use crate::parser::ast::SyntaxError;

fn syntax(
    line: usize,
    col: usize,
    message: impl Into<String>,
    line_map: &Arc<LineMap>,
) -> ParseError {
    ParseError::Syntax(SyntaxError::synth(
        line,
        col,
        message,
        Arc::clone(line_map),
    ))
}
```

(c) Replace every existing `ParseError::Syntax { line, col, message }` site with `syntax(line, col, message, &line_map)`. The sites are (search for `ParseError::Syntax`):

- `parse_pi_expr` — needs `&Arc<LineMap>` parameter; thread through callers.
- `parse_tags` — needs `&Arc<LineMap>` parameter; thread through callers.
- `parse_block` — needs `&Arc<LineMap>` parameter; already has `line_map` in scope (re-use).
- `parse_line` — needs `&Arc<LineMap>` parameter (today it takes `_line_map: &LineMap`; change to `line_map: &Arc<LineMap>` and use it).
- `parse_head` — needs `&Arc<LineMap>` parameter.

Concretely, change every signature like `fn parse_X(... ) -> Result<Y, ParseError>` to take an additional `line_map: &Arc<LineMap>` parameter, and update each `ParseError::Syntax { line, col, message }` call to `syntax(line, col, message, line_map)`. Update callers to pass it through.

This is a mechanical refactor. After it compiles, behavior is unchanged but the error variant uses `SyntaxError` for syntax errors.

(d) Remove the now-unused `_line_map: &LineMap` parameter rename — `parse_line` now takes `line_map: &Arc<LineMap>` instead.

- [ ] **Step 7: Run all tests in `ppvm-stim`**

Run: `cargo test -p ppvm-stim`
Expected: all PASS, including the new `syntax_error_includes_line_and_col` test (because hand-written `parse_line` reports `line: line_no` for `H abc` at line 2 with `col: 1`, so the formatted message contains "line 2" and "col").

If `syntax_error_includes_line_and_col` fails because the message says "col 1" — that's fine, it still contains "col". If any other test fails because of `==` on `ParseError`, fix by switching to `matches!`.

- [ ] **Step 8: Run the wider workspace tests**

Run: `cargo test --workspace`
Expected: all PASS. This catches any external consumer of `ParseError` that relied on `PartialEq`.

If anything fails, fix it (typically: replace `==` checks with `matches!`).

- [ ] **Step 9: Format and commit**

```bash
cd /Users/david/git/claude-repos/ppvm
cargo fmt --all
git add crates/ppvm-stim/src/parser/ast.rs \
        crates/ppvm-stim/src/parser/mod.rs \
        crates/ppvm-stim/tests/parser_errors.rs
git commit -m "$(cat <<'EOF'
Switch ParseError::Syntax to tuple variant carrying SyntaxError

Replaces the struct variant `ParseError::Syntax { line, col, message }`
with a tuple variant `ParseError::Syntax(SyntaxError)` that wraps a
chumsky `Rich<'static, char>` and a shared `Arc<LineMap>`. Display
formats `line:col` from the LineMap so the new shape matches the
existing typed-variant error messages.

The hand-written parser is still in use; its existing Syntax sites
synthesise a Rich via `SyntaxError::synth`. The chumsky grammar
introduced in subsequent commits will populate Rich directly.

PartialEq is dropped from ParseError (Rich is not PartialEq); existing
tests use `matches!`, so this is non-breaking.
EOF
)"
```

---

### Task 3: Define internal raw types `RawSyntaxNode`, `RawTag`, `RawTagParam`, `RawTarget`

Internal `pub(crate)` types that the chumsky grammar will produce and the validator will consume. Defined in `mod.rs` so both `grammar.rs` and the validator can `use super::*`.

**Files:**
- Modify: `crates/ppvm-stim/src/parser/mod.rs` (add types near the top, below `LineMap`)

- [ ] **Step 1: Add the raw types**

Add this block in `crates/ppvm-stim/src/parser/mod.rs`, immediately after the `LineMap` impl block:

```rust
use chumsky::span::SimpleSpan;

/// Raw syntactic tree produced by the chumsky grammar before
/// table-driven validation. `pub(crate)` because it is plumbing
/// between `grammar.rs` and the validator post-pass; not part of the
/// public API.
#[derive(Debug, Clone)]
pub(crate) enum RawSyntaxNode {
    Instruction {
        name: String,
        tags: Vec<RawTag>,
        args: Vec<f64>,
        targets: Vec<RawTarget>,
        span: SimpleSpan<usize>,
    },
    Repeat {
        count: u64,
        body: Vec<RawSyntaxNode>,
        span: SimpleSpan<usize>,
    },
}

#[derive(Debug, Clone)]
pub(crate) struct RawTag {
    pub name: String,
    pub params: Vec<RawTagParam>,
}

#[derive(Debug, Clone)]
pub(crate) enum RawTagParam {
    Positional(f64),
    Named { key: String, value: f64 },
}

#[derive(Debug, Clone)]
pub(crate) struct RawTarget {
    pub text: String,
    pub span: SimpleSpan<usize>,
}
```

- [ ] **Step 2: Verify it compiles**

Run: `cargo build -p ppvm-stim`
Expected: builds clean (the types are unused — `pub(crate)` may yield "never used" warnings for now; that is OK and resolves in later tasks).

If a "never used" warning is emitted, do not silence it with `#[allow(dead_code)]` — the next tasks consume these types.

- [ ] **Step 3: Run all tests**

Run: `cargo test -p ppvm-stim`
Expected: all PASS (no behavioral change).

- [ ] **Step 4: Format and commit**

```bash
cd /Users/david/git/claude-repos/ppvm
cargo fmt --all
git add crates/ppvm-stim/src/parser/mod.rs
git commit -m "$(cat <<'EOF'
Define raw syntactic-tree types for the chumsky grammar

Adds internal RawSyntaxNode / RawTag / RawTagParam / RawTarget types
that the chumsky grammar will produce and the validator post-pass
will consume. Plumbing only; no behavior change.
EOF
)"
```

---

### Task 4: Stand up `parser/grammar.rs` with foundation combinators (pad, ident, signed_float, pi_expr)

Create `grammar.rs` with the smallest set of combinators and inline tests that exercise the chumsky 0.12 API surface we depend on. This is the smoke-test moment — if the chumsky API differs from this plan, surface it now.

**Files:**
- Create: `crates/ppvm-stim/src/parser/grammar.rs`
- Modify: `crates/ppvm-stim/src/parser/mod.rs` (add `mod grammar;` declaration)

- [ ] **Step 1: Create `grammar.rs` with foundations + inline tests**

Create `crates/ppvm-stim/src/parser/grammar.rs` with the following content:

```rust
//! Chumsky 0.12 grammar for Stim source.
//!
//! Reads top-to-bottom: whitespace/comments → numbers → pi-expressions →
//! identifiers → tags → args → targets → instruction line → REPEAT block →
//! program. Pure syntax; no table lookups.

use chumsky::error::Rich;
use chumsky::extra;
use chumsky::prelude::*;

type Extra<'src> = extra::Err<Rich<'src, char>>;

/// `# ...` comment, stopping before `\n` if present.
fn line_comment<'src>() -> impl Parser<'src, &'src str, (), Extra<'src>> + Clone {
    just('#')
        .ignore_then(any().filter(|c: &char| *c != '\n').repeated())
        .ignored()
}

/// Pad: zero or more whitespace characters or `# ...` comments.
/// Includes newlines. Used between instructions.
pub(super) fn pad<'src>() -> impl Parser<'src, &'src str, (), Extra<'src>> + Clone {
    let ws = any().filter(|c: &char| c.is_whitespace()).ignored();
    choice((line_comment(), ws)).repeated().ignored()
}

/// Inline pad: spaces/tabs/CR only. Excludes comments and `\n`.
pub(super) fn inline_pad<'src>() -> impl Parser<'src, &'src str, (), Extra<'src>> + Clone {
    let ws = any()
        .filter(|c: &char| matches!(*c, ' ' | '\t' | '\r'))
        .ignored();
    ws.repeated().ignored()
}

/// At least one inline whitespace character. Used before each target.
pub(super) fn inline_ws1<'src>() -> impl Parser<'src, &'src str, (), Extra<'src>> + Clone {
    any()
        .filter(|c: &char| matches!(*c, ' ' | '\t' | '\r'))
        .repeated()
        .at_least(1)
        .ignored()
}

/// Optional trailing spaces/tabs plus an optional line comment.
pub(super) fn trailing_pad<'src>() -> impl Parser<'src, &'src str, (), Extra<'src>> + Clone {
    inline_pad().then(line_comment().or_not()).ignored()
}

/// Identifier: `[A-Za-z_][A-Za-z0-9_]*`. Returns owned `String`.
pub(super) fn ident<'src>() -> impl Parser<'src, &'src str, String, Extra<'src>> + Clone {
    any()
        .filter(|c: &char| c.is_ascii_alphabetic() || *c == '_')
        .then(
            any()
                .filter(|c: &char| c.is_ascii_alphanumeric() || *c == '_')
                .repeated(),
        )
        .to_slice()
        .map(|s: &str| s.to_string())
}

/// Signed float: `[+-]? digits ('.' digits)? ([eE] [+-]? digits)?`.
pub(super) fn signed_float<'src>() -> impl Parser<'src, &'src str, f64, Extra<'src>> + Clone {
    let digits = any()
        .filter(|c: &char| c.is_ascii_digit())
        .repeated()
        .at_least(1);
    let int_part = digits.clone();
    let frac_part = just('.').then(digits.clone());
    let exp_part = one_of("eE")
        .then(one_of("+-").or_not())
        .then(digits.clone());
    one_of("+-")
        .or_not()
        .then(int_part)
        .then(frac_part.or_not())
        .then(exp_part.or_not())
        .to_slice()
        .map(|s: &str| s.parse::<f64>().expect("validated by combinator shape"))
}

/// Pi-expression: `pi`, `<num>*pi`, or plain number. Evaluates to f64.
pub(super) fn pi_expr<'src>() -> impl Parser<'src, &'src str, f64, Extra<'src>> + Clone {
    let pi_kw = just("pi").to(std::f64::consts::PI);
    let num_then_pi = signed_float()
        .then(just("*pi").or_not())
        .map(|(n, suffix)| {
            if suffix.is_some() {
                n * std::f64::consts::PI
            } else {
                n
            }
        });
    choice((pi_kw, num_then_pi))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run<'src, T>(p: impl Parser<'src, &'src str, T, Extra<'src>>, src: &'src str) -> T {
        p.parse(src).into_result().expect("parse failed")
    }

    #[test]
    fn ident_matches_alpha_then_alphanumeric() {
        assert_eq!(run(ident(), "H"), "H");
        assert_eq!(run(ident(), "DEPOLARIZE1"), "DEPOLARIZE1");
        assert_eq!(run(ident(), "_x"), "_x");
        assert_eq!(run(ident(), "R_X"), "R_X");
    }

    #[test]
    fn signed_float_parses_common_shapes() {
        assert_eq!(run(signed_float(), "0"), 0.0);
        assert_eq!(run(signed_float(), "0.5"), 0.5);
        assert_eq!(run(signed_float(), "-0.5"), -0.5);
        assert_eq!(run(signed_float(), "+1.0e-3"), 1.0e-3);
        assert_eq!(run(signed_float(), "42"), 42.0);
    }

    #[test]
    fn pi_expr_parses_pi_keyword_coeff_and_plain_number() {
        assert_eq!(run(pi_expr(), "pi"), std::f64::consts::PI);
        assert_eq!(run(pi_expr(), "0.5*pi"), 0.5 * std::f64::consts::PI);
        assert_eq!(run(pi_expr(), "-2*pi"), -2.0 * std::f64::consts::PI);
        assert_eq!(run(pi_expr(), "0.5"), 0.5);
    }

    #[test]
    fn pad_eats_comments_and_newlines() {
        // pad's output is (); we just verify it doesn't error.
        let p = pad().then_ignore(end());
        assert!(p.parse("").into_result().is_ok());
        assert!(p.parse("   \n\t# comment\n").into_result().is_ok());
    }

    #[test]
    fn inline_pad_does_not_consume_comment_or_newline() {
        // inline_pad should leave both comment starts and newlines to the
        // line-ending parser.
        let p = inline_pad().then_ignore(just('#'));
        assert!(p.parse(" \t#").into_result().is_ok());
        let p = inline_pad().then_ignore(just('\n'));
        assert!(p.parse(" \t\n").into_result().is_ok());
    }

    #[test]
    fn trailing_pad_consumes_trailing_comment() {
        let p = trailing_pad()
            .then_ignore(just('\n').or_not())
            .then_ignore(end());
        assert!(p.parse(" # comment\n").into_result().is_ok());
        assert!(p.parse("   ").into_result().is_ok());
    }
}
```

- [ ] **Step 2: Wire `mod grammar;` in `mod.rs`**

In `crates/ppvm-stim/src/parser/mod.rs`, add at the top of the file (after the existing `pub mod ast;` and `pub mod table;`):

```rust
mod grammar;
```

- [ ] **Step 3: Run grammar tests**

Run: `cargo test -p ppvm-stim parser::grammar::tests`
Expected: all 6 tests PASS.

If any chumsky API spelling is wrong (e.g. `to_slice` is named differently in your chumsky 0.12), fix the call sites here. The shape is correct; only method names may need adjustment.

- [ ] **Step 4: Run all crate tests to confirm no regressions**

Run: `cargo test -p ppvm-stim`
Expected: all PASS.

- [ ] **Step 5: Format and commit**

```bash
cd /Users/david/git/claude-repos/ppvm
cargo fmt --all
git add crates/ppvm-stim/src/parser/grammar.rs crates/ppvm-stim/src/parser/mod.rs
git commit -m "$(cat <<'EOF'
Add chumsky grammar foundation combinators

Introduces parser/grammar.rs with pad, inline_pad, ident, signed_float,
and pi_expr combinators plus inline unit tests. Validates the chumsky
0.12 API shape we depend on before building the full grammar.
EOF
)"
```

---

### Task 5: Add tag and arg combinators (`tag_param`, `tag`, `tags_block`, `args_block`)

**Files:**
- Modify: `crates/ppvm-stim/src/parser/grammar.rs` (append below `pi_expr`)

- [ ] **Step 1: Add the failing tests**

Append to the existing `#[cfg(test)] mod tests` block in `grammar.rs`:

```rust
    #[test]
    fn tag_with_no_params() {
        let t = run(tag(), "T");
        assert_eq!(t.name, "T");
        assert!(t.params.is_empty());
    }

    #[test]
    fn tag_with_positional_params() {
        let t = run(tag(), "R(0.5, 1.0)");
        assert_eq!(t.name, "R");
        assert_eq!(t.params.len(), 2);
        assert!(matches!(&t.params[0], super::super::RawTagParam::Positional(v) if (*v - 0.5).abs() < 1e-12));
    }

    #[test]
    fn tag_with_named_param() {
        let t = run(tag(), "R_X(theta=0.5*pi)");
        assert_eq!(t.name, "R_X");
        assert_eq!(t.params.len(), 1);
        match &t.params[0] {
            super::super::RawTagParam::Named { key, value } => {
                assert_eq!(key, "theta");
                assert!((value - 0.5 * std::f64::consts::PI).abs() < 1e-12);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn tags_block_parses_multiple_tags() {
        let ts = run(tags_block(), "[T, R(0.5)]");
        assert_eq!(ts.len(), 2);
        assert_eq!(ts[0].name, "T");
        assert_eq!(ts[1].name, "R");
    }

    #[test]
    fn args_block_parses_csv_floats() {
        let a = run(args_block(), "(0.1, 0.2, 0.3)");
        assert_eq!(a, vec![0.1, 0.2, 0.3]);
    }

    #[test]
    fn args_block_with_pi_exprs() {
        let a = run(args_block(), "(pi, 0.5*pi)");
        assert!((a[0] - std::f64::consts::PI).abs() < 1e-12);
        assert!((a[1] - 0.5 * std::f64::consts::PI).abs() < 1e-12);
    }
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p ppvm-stim parser::grammar::tests`
Expected: compile error — `tag`, `tags_block`, `args_block`, `RawTagParam` not in scope.

- [ ] **Step 3: Add the combinators**

Append to `grammar.rs` (above the `#[cfg(test)] mod tests` block):

```rust
use super::{RawTag, RawTagParam};

/// `<ident>=<pi_expr>` (Named) or `<pi_expr>` (Positional).
pub(super) fn tag_param<'src>() -> impl Parser<'src, &'src str, RawTagParam, Extra<'src>> + Clone {
    let named = ident()
        .then_ignore(inline_pad())
        .then_ignore(just('='))
        .then_ignore(inline_pad())
        .then(pi_expr())
        .map(|(key, value)| RawTagParam::Named { key, value });
    let positional = pi_expr().map(RawTagParam::Positional);
    choice((named, positional))
}

/// Tag: `<ident>` or `<ident>(<tag_param>, ...)`.
pub(super) fn tag<'src>() -> impl Parser<'src, &'src str, RawTag, Extra<'src>> + Clone {
    let params = tag_param()
        .separated_by(inline_pad().then(just(',')).then(inline_pad()))
        .allow_trailing()
        .collect::<Vec<_>>()
        .delimited_by(
            just('(').then(inline_pad()),
            inline_pad().then(just(')')),
        );
    ident()
        .then(params.or_not())
        .map(|(name, params)| RawTag {
            name,
            params: params.unwrap_or_default(),
        })
}

/// `[tag, tag, ...]`.
pub(super) fn tags_block<'src>() -> impl Parser<'src, &'src str, Vec<RawTag>, Extra<'src>> + Clone {
    tag()
        .separated_by(inline_pad().then(just(',')).then(inline_pad()))
        .allow_trailing()
        .collect::<Vec<_>>()
        .delimited_by(
            just('[').then(inline_pad()),
            inline_pad().then(just(']')),
        )
}

/// `(pi_expr, pi_expr, ...)`.
pub(super) fn args_block<'src>() -> impl Parser<'src, &'src str, Vec<f64>, Extra<'src>> + Clone {
    pi_expr()
        .separated_by(inline_pad().then(just(',')).then(inline_pad()))
        .allow_trailing()
        .collect::<Vec<_>>()
        .delimited_by(
            just('(').then(inline_pad()),
            inline_pad().then(just(')')),
        )
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test -p ppvm-stim parser::grammar::tests`
Expected: all PASS.

- [ ] **Step 5: Run all crate tests**

Run: `cargo test -p ppvm-stim`
Expected: all PASS.

- [ ] **Step 6: Format and commit**

```bash
cd /Users/david/git/claude-repos/ppvm
cargo fmt --all
git add crates/ppvm-stim/src/parser/grammar.rs
git commit -m "$(cat <<'EOF'
Add tag and arg combinators to chumsky grammar

Implements tag_param (named/positional), tag, tags_block, and
args_block combinators with unit tests covering the Stim shapes:
plain idents, parenthesised positional params, named params, and
CSV-of-pi-expr arg lists.
EOF
)"
```

---

### Task 6: Add target and instruction-head combinators (`target_lexeme`, `instruction_head`)

**Files:**
- Modify: `crates/ppvm-stim/src/parser/grammar.rs` (append)

- [ ] **Step 1: Add tests**

Append to the inline tests:

```rust
    #[test]
    fn target_lexeme_reads_a_non_whitespace_run() {
        let t = run(target_lexeme(), "0");
        assert_eq!(t.text, "0");
        let t = run(target_lexeme(), "rec[-1]");
        assert_eq!(t.text, "rec[-1]");
    }

    #[test]
    fn target_lexeme_stops_at_brace() {
        // target_lexeme should not consume `}`.
        let p = target_lexeme().then_ignore(just('}'));
        let t = p.parse("0}").into_result().expect("parse failed");
        assert_eq!(t.text, "0");
    }

    #[test]
    fn instruction_head_with_tags_and_args() {
        let (name, tags, args, _span) = run(instruction_head(), "S[T](0.5)");
        assert_eq!(name, "S");
        assert_eq!(tags.len(), 1);
        assert_eq!(tags[0].name, "T");
        assert_eq!(args, vec![0.5]);
    }

    #[test]
    fn instruction_head_no_tags_no_args() {
        let (name, tags, args, _span) = run(instruction_head(), "H");
        assert_eq!(name, "H");
        assert!(tags.is_empty());
        assert!(args.is_empty());
    }
```

- [ ] **Step 2: Verify tests fail**

Run: `cargo test -p ppvm-stim parser::grammar::tests`
Expected: compile error — `target_lexeme`, `instruction_head` undefined.

- [ ] **Step 3: Add the combinators**

Append to `grammar.rs` (still above the test module):

```rust
use chumsky::span::SimpleSpan;
use super::RawTarget;

/// One non-whitespace, non-`#`, non-`{`/`}` lexeme. Captures the span
/// so the validator can derive a line number for invalid-target errors.
pub(super) fn target_lexeme<'src>() -> impl Parser<'src, &'src str, RawTarget, Extra<'src>> + Clone
{
    any()
        .filter(|c: &char| !c.is_whitespace() && *c != '#' && *c != '{' && *c != '}')
        .repeated()
        .at_least(1)
        .to_slice()
        .map_with(|s: &str, e| RawTarget {
            text: s.to_string(),
            span: e.span(),
        })
}

/// `<ident> [<tags>]? (<args>)?`. Returns name, tags, args, and the
/// span of the identifier (used for line-number reporting).
pub(super) fn instruction_head<'src>(
) -> impl Parser<'src, &'src str, (String, Vec<RawTag>, Vec<f64>, SimpleSpan<usize>), Extra<'src>>
       + Clone {
    ident()
        .map_with(|name, e| (name, e.span()))
        .then(tags_block().or_not())
        .then(args_block().or_not())
        .map(|(((name, span), tags), args)| {
            (name, tags.unwrap_or_default(), args.unwrap_or_default(), span)
        })
}
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p ppvm-stim parser::grammar::tests`
Expected: all PASS.

If `map_with` does not work in your chumsky 0.12 release, the older spelling is `map_with_span(|out, span| ...)` returning `(out, span)`. Use whichever compiles. The data shape — `(name, span)` then chained — must remain.

- [ ] **Step 5: Run all crate tests**

Run: `cargo test -p ppvm-stim`
Expected: all PASS.

- [ ] **Step 6: Format and commit**

```bash
cd /Users/david/git/claude-repos/ppvm
cargo fmt --all
git add crates/ppvm-stim/src/parser/grammar.rs
git commit -m "$(cat <<'EOF'
Add target_lexeme and instruction_head combinators

Reads target lexemes as raw text + span (validator handles usize
parsing and annotation tolerance), and reads instruction heads of
shape `<ident> [tags]? (args)?` with the ident span captured for
line-number reporting.
EOF
)"
```

---

### Task 7: Add `instruction_line`, `repeat_block`, and `program` combinators

**Files:**
- Modify: `crates/ppvm-stim/src/parser/grammar.rs` (append)

- [ ] **Step 1: Add tests**

Append to the inline tests:

```rust
    use super::super::RawSyntaxNode;

    #[test]
    fn instruction_line_with_targets() {
        let n = run(instruction_line(), "CX 0 1 2 3");
        match n {
            RawSyntaxNode::Instruction { name, targets, .. } => {
                assert_eq!(name, "CX");
                let texts: Vec<_> = targets.iter().map(|t| t.text.clone()).collect();
                assert_eq!(texts, vec!["0", "1", "2", "3"]);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn instruction_line_no_targets() {
        let n = run(instruction_line(), "TICK");
        match n {
            RawSyntaxNode::Instruction { name, targets, .. } => {
                assert_eq!(name, "TICK");
                assert!(targets.is_empty());
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn program_parses_multiline() {
        let p = program_parser();
        let nodes = p.parse("X 0\nY 1\n").into_result().expect("parse failed");
        assert_eq!(nodes.len(), 2);
    }

    #[test]
    fn program_handles_blank_lines_and_comments() {
        let p = program_parser();
        let nodes = p
            .parse("\n# header\nX 0\n# mid\nY 1\n")
            .into_result()
            .expect("parse failed");
        assert_eq!(nodes.len(), 2);
    }

    #[test]
    fn program_handles_trailing_comments() {
        let p = program_parser();
        let nodes = p
            .parse("X 0 # flip\nY 1 # measure\n")
            .into_result()
            .expect("parse failed");
        assert_eq!(nodes.len(), 2);
    }

    #[test]
    fn program_parses_repeat_block() {
        let p = program_parser();
        let nodes = p
            .parse("REPEAT 3 {\n    X 0\n    M 0\n}\n")
            .into_result()
            .expect("parse failed");
        assert_eq!(nodes.len(), 1);
        match &nodes[0] {
            RawSyntaxNode::Repeat { count, body, .. } => {
                assert_eq!(*count, 3);
                assert_eq!(body.len(), 2);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn program_parses_nested_repeat() {
        let p = program_parser();
        let nodes = p
            .parse("REPEAT 2 {\n  REPEAT 3 {\n    H 0\n  }\n}\n")
            .into_result()
            .expect("parse failed");
        match &nodes[0] {
            RawSyntaxNode::Repeat { body, .. } => match &body[0] {
                RawSyntaxNode::Repeat { count, .. } => assert_eq!(*count, 3),
                other => panic!("{other:?}"),
            },
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn program_parses_one_line_repeat() {
        let p = program_parser();
        let nodes = p
            .parse("REPEAT 5 { H 0 }")
            .into_result()
            .expect("parse failed");
        match &nodes[0] {
            RawSyntaxNode::Repeat { count, body, .. } => {
                assert_eq!(*count, 5);
                assert_eq!(body.len(), 1);
            }
            other => panic!("{other:?}"),
        }
    }
```

- [ ] **Step 2: Verify tests fail**

Run: `cargo test -p ppvm-stim parser::grammar::tests`
Expected: compile error — `instruction_line`, `program_parser`, `RawSyntaxNode` (in module path) undefined.

- [ ] **Step 3: Add the combinators**

Append to `grammar.rs` (above the test module):

```rust
use super::RawSyntaxNode;

/// End of an instruction: optional trailing spaces/comment followed by
/// newline, EOF, or a `}` that belongs to the enclosing REPEAT parser.
pub(super) fn line_end<'src>() -> impl Parser<'src, &'src str, (), Extra<'src>> + Clone {
    trailing_pad()
        .then(choice((
            just('\n').ignored(),
            just('}').rewind().ignored(),
            end(),
        )))
        .ignored()
}

/// Instruction line: head + space-separated raw targets, terminated by
/// a newline / EOF / `}`. Newline is consumed; `}` is only peeked so
/// the enclosing REPEAT parser can consume it.
pub(super) fn instruction_line<'src>(
) -> impl Parser<'src, &'src str, RawSyntaxNode, Extra<'src>> + Clone {
    instruction_head()
        .then(
            inline_ws1()
                .ignore_then(target_lexeme())
                .repeated()
                .collect::<Vec<RawTarget>>(),
        )
        .map(|((name, tags, args, span), targets)| RawSyntaxNode::Instruction {
            name,
            tags,
            args,
            targets,
            span,
        })
        .then_ignore(line_end())
}

/// `REPEAT <count> { <body> }`. `body` is the parser for a list of
/// nodes (instructions and nested REPEATs).
fn repeat_block<'src>(
    body: impl Parser<'src, &'src str, Vec<RawSyntaxNode>, Extra<'src>> + Clone + 'src,
) -> impl Parser<'src, &'src str, RawSyntaxNode, Extra<'src>> + Clone + 'src {
    let digits = any()
        .filter(|c: &char| c.is_ascii_digit())
        .repeated()
        .at_least(1)
        .to_slice()
        .map(|s: &str| s.parse::<u64>().expect("digit-only"));
    just("REPEAT")
        .map_with(|_, e| e.span())
        .then_ignore(inline_pad())
        .then(digits)
        .then_ignore(inline_pad())
        .then_ignore(just('{'))
        .then_ignore(pad())
        .then(body)
        .then_ignore(pad())
        .then_ignore(just('}'))
        .map(|((span, count), body)| RawSyntaxNode::Repeat { count, body, span })
}

/// Top-level program parser. Recursively defines the body shared by
/// the program and REPEAT blocks.
pub(crate) fn program_parser<'src>(
) -> impl Parser<'src, &'src str, Vec<RawSyntaxNode>, Extra<'src>> {
    recursive(|body| {
        let item = choice((repeat_block(body.clone()), instruction_line()));
        pad()
            .ignore_then(item)
            .repeated()
            .collect::<Vec<RawSyntaxNode>>()
            .then_ignore(pad())
    })
    .then_ignore(end())
}
```

- [ ] **Step 4: Run grammar tests**

Run: `cargo test -p ppvm-stim parser::grammar::tests`
Expected: all PASS.

If `recursive` requires a different closure shape in your chumsky 0.12, adapt. The intended structure is: `recursive(|body| <body parser>)` where `body` is the parser being defined and is referenced inside `repeat_block(body.clone())`.

- [ ] **Step 5: Run all crate tests**

Run: `cargo test -p ppvm-stim`
Expected: all PASS (existing parser still in use; grammar tests are isolated).

- [ ] **Step 6: Format and commit**

```bash
cd /Users/david/git/claude-repos/ppvm
cargo fmt --all
git add crates/ppvm-stim/src/parser/grammar.rs
git commit -m "$(cat <<'EOF'
Add instruction_line, repeat_block, and program combinators

Completes the chumsky grammar with the top-level program parser.
Newlines and trailing comments terminate instructions via line_end,
matching Stim's line-oriented format. REPEAT blocks recurse via
chumsky's `recursive` combinator.
EOF
)"
```

---

### Task 8: Implement `validate_program` in `mod.rs`

The validator walks `Vec<RawSyntaxNode>` and produces `Vec<RawInstruction>`, doing table lookup, arg-count check, target-arity check, and target `usize` parsing with annotation tolerance. This logic mirrors the existing `parse_line` code; the new function operates on a pre-parsed tree instead of raw lines.

**Files:**
- Modify: `crates/ppvm-stim/src/parser/mod.rs` (add `validate_program` near the end of the file, **without** wiring `parse()` to it yet)

- [ ] **Step 1: Add a unit test for `validate_program`**

Append to the existing `#[cfg(test)] mod tests` block (or `line_map_tests` block — it doesn't matter; rename the inline mod to `mod_tests` if there's a name collision) in `mod.rs`:

```rust
#[cfg(test)]
mod validate_tests {
    use super::*;
    use crate::parser::ast::{GateName, MeasureName, RawInstruction};
    use chumsky::span::SimpleSpan;
    use std::sync::Arc;

    fn lm() -> Arc<LineMap> {
        Arc::new(LineMap::new("H 0\nM 0"))
    }

    fn instr(name: &str, args: Vec<f64>, targets: Vec<&str>, span: (usize, usize)) -> RawSyntaxNode {
        RawSyntaxNode::Instruction {
            name: name.to_string(),
            tags: vec![],
            args,
            targets: targets
                .into_iter()
                .map(|t| RawTarget {
                    text: t.to_string(),
                    span: SimpleSpan::from(span.0..span.1),
                })
                .collect(),
            span: SimpleSpan::from(span.0..span.1),
        }
    }

    #[test]
    fn validates_simple_gate() {
        let nodes = vec![instr("H", vec![], vec!["0"], (0, 1))];
        let result = validate_program(nodes, &lm()).unwrap();
        assert!(matches!(
            &result[0],
            RawInstruction::Gate { name: GateName::H, line: 1, .. }
        ));
    }

    #[test]
    fn validates_measure() {
        let nodes = vec![instr("M", vec![], vec!["0"], (4, 5))];
        let result = validate_program(nodes, &lm()).unwrap();
        assert!(matches!(
            &result[0],
            RawInstruction::Measure { name: MeasureName::M, line: 2, .. }
        ));
    }

    #[test]
    fn unknown_instruction_errors() {
        let nodes = vec![instr("FROBNICATE", vec![], vec!["0"], (0, 10))];
        let err = validate_program(nodes, &lm()).unwrap_err();
        match err {
            ParseError::UnknownInstruction { name, line } => {
                assert_eq!(name, "FROBNICATE");
                assert_eq!(line, 1);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn arg_count_errors() {
        // DEPOLARIZE1 expects exactly 1 arg.
        let nodes = vec![instr("DEPOLARIZE1", vec![0.1, 0.2], vec!["0"], (0, 11))];
        let err = validate_program(nodes, &lm()).unwrap_err();
        assert!(matches!(err, ParseError::ArgCount { .. }));
    }

    #[test]
    fn target_pair_errors() {
        let nodes = vec![instr("CX", vec![], vec!["0", "1", "2"], (0, 2))];
        let err = validate_program(nodes, &lm()).unwrap_err();
        assert!(matches!(err, ParseError::TargetCount { divisor: 2, found: 3, .. }));
    }

    #[test]
    fn invalid_target_for_gate_is_syntax_error() {
        let nodes = vec![instr("H", vec![], vec!["abc"], (0, 1))];
        let err = validate_program(nodes, &lm()).unwrap_err();
        assert!(matches!(err, ParseError::Syntax(_)));
    }

    #[test]
    fn annotation_tolerates_non_numeric_targets() {
        // DETECTOR is an annotation: rec[-1] should be silently dropped.
        let nodes = vec![instr("DETECTOR", vec![], vec!["rec[-1]"], (0, 8))];
        let result = validate_program(nodes, &lm()).unwrap();
        assert!(matches!(&result[0], RawInstruction::Annotation { .. }));
    }
}
```

- [ ] **Step 2: Verify tests fail**

Run: `cargo test -p ppvm-stim validate_tests`
Expected: compile error — `validate_program` undefined.

- [ ] **Step 3: Implement `validate_program`**

Append to `crates/ppvm-stim/src/parser/mod.rs` (place above the `#[cfg(test)] mod` blocks):

```rust
use crate::parser::ast::{RawInstruction, Tag, TagParam};
use crate::parser::table::{TableEntry, lookup};

/// Walk the raw syntactic tree and emit validated instructions.
fn validate_program(
    nodes: Vec<RawSyntaxNode>,
    line_map: &Arc<LineMap>,
) -> Result<Vec<RawInstruction>, ParseError> {
    let mut out = Vec::with_capacity(nodes.len());
    for node in nodes {
        out.push(validate_node(node, line_map)?);
    }
    Ok(out)
}

fn validate_node(
    node: RawSyntaxNode,
    line_map: &Arc<LineMap>,
) -> Result<RawInstruction, ParseError> {
    match node {
        RawSyntaxNode::Instruction {
            name,
            tags,
            args,
            targets,
            span,
        } => {
            let line = line_map.line_of(span.start);
            let entry = lookup(&name).ok_or(ParseError::UnknownInstruction {
                name: name.clone(),
                line,
            })?;
            let (arg_rule, target_rule, canonical) = arity_of(entry);

            // Convert RawTag → Tag (no validation; tag semantics live in `normalize`).
            let tags: Vec<Tag> = tags
                .into_iter()
                .map(|t| Tag {
                    name: t.name,
                    params: t
                        .params
                        .into_iter()
                        .map(|p| match p {
                            RawTagParam::Named { key, value } => TagParam::Named { key, value },
                            RawTagParam::Positional(v) => TagParam::Positional(v),
                        })
                        .collect(),
                })
                .collect();

            // Arg-count validation.
            let skip_arg_validation = matches!(entry, TableEntry::Annotation { .. })
                || matches!(
                    entry,
                    TableEntry::Noise {
                        args: ArgCount::None,
                        ..
                    }
                );
            if !skip_arg_validation {
                match arg_rule {
                    ArgCount::None if !args.is_empty() => {
                        return Err(ParseError::ArgCount {
                            name: canonical.to_string(),
                            expected: 0,
                            found: args.len(),
                            line,
                        });
                    }
                    ArgCount::Exact(n) if args.len() != n => {
                        return Err(ParseError::ArgCount {
                            name: canonical.to_string(),
                            expected: n,
                            found: args.len(),
                            line,
                        });
                    }
                    ArgCount::Optional(n) if !args.is_empty() && args.len() != n => {
                        return Err(ParseError::ArgCount {
                            name: canonical.to_string(),
                            expected: n,
                            found: args.len(),
                            line,
                        });
                    }
                    _ => {}
                }
            }

            // Target validation.
            let tolerate_non_numeric = matches!(entry, TableEntry::Annotation { .. });
            let mut numeric_targets: Vec<usize> = Vec::with_capacity(targets.len());
            for t in &targets {
                match t.text.parse::<usize>() {
                    Ok(n) => numeric_targets.push(n),
                    Err(_) if tolerate_non_numeric => {}
                    Err(_) => {
                        let (l, c) = line_map.line_col(t.span.start);
                        return Err(syntax(
                            l,
                            c,
                            format!("invalid target {:?}", t.text),
                            line_map,
                        ));
                    }
                }
            }

            // Target-arity validation.
            match target_rule {
                TargetArity::Any => {}
                TargetArity::AtLeastOne if numeric_targets.is_empty() => {
                    return Err(ParseError::TargetCount {
                        name: canonical.to_string(),
                        divisor: 1,
                        found: 0,
                        line,
                    });
                }
                TargetArity::Pairs
                    if !numeric_targets.len().is_multiple_of(2) || numeric_targets.is_empty() =>
                {
                    return Err(ParseError::TargetCount {
                        name: canonical.to_string(),
                        divisor: 2,
                        found: numeric_targets.len(),
                        line,
                    });
                }
                TargetArity::Quadruples
                    if !numeric_targets.len().is_multiple_of(4) || numeric_targets.is_empty() =>
                {
                    return Err(ParseError::TargetCount {
                        name: canonical.to_string(),
                        divisor: 4,
                        found: numeric_targets.len(),
                        line,
                    });
                }
                _ => {}
            }

            Ok(build_instruction_v2(entry, tags, args, numeric_targets, line))
        }
        RawSyntaxNode::Repeat { count, body, span } => {
            let line = line_map.line_of(span.start);
            let body_validated = validate_program(body, line_map)?;
            Ok(RawInstruction::Repeat {
                count,
                body: body_validated,
                line,
            })
        }
    }
}

fn build_instruction_v2(
    entry: TableEntry,
    tags: Vec<Tag>,
    args: Vec<f64>,
    targets: Vec<usize>,
    line: usize,
) -> RawInstruction {
    match entry {
        TableEntry::Gate { name, .. } => RawInstruction::Gate {
            name,
            tags,
            args,
            targets,
            line,
        },
        TableEntry::Noise { name, .. } => RawInstruction::Noise {
            name,
            tags,
            args,
            targets,
            line,
        },
        TableEntry::Measure { name, .. } => RawInstruction::Measure {
            name,
            tags,
            args,
            targets,
            line,
        },
        TableEntry::Annotation { kind, .. } => RawInstruction::Annotation {
            kind,
            args,
            targets,
            line,
        },
    }
}
```

`build_instruction_v2` is named `_v2` only to avoid colliding with the existing `build_instruction` (which still exists for the hand-written parser). Task 9 deletes the `_v2` suffix and the old `build_instruction` together.

`arity_of` already exists in `mod.rs` — re-use it.

Note the imports: `GateName`, `MeasureName`, `NoiseName`, `AnnotationKind` are imported but may already be in scope. Adjust the `use` block to avoid duplicate-import warnings.

- [ ] **Step 4: Run validator tests**

Run: `cargo test -p ppvm-stim validate_tests`
Expected: all 7 tests PASS.

- [ ] **Step 5: Run all crate tests**

Run: `cargo test -p ppvm-stim`
Expected: all PASS. The hand-written parser is still in use for `parse()`; `validate_program` is exercised only by the new unit tests.

- [ ] **Step 6: Format and commit**

```bash
cd /Users/david/git/claude-repos/ppvm
cargo fmt --all
git add crates/ppvm-stim/src/parser/mod.rs
git commit -m "$(cat <<'EOF'
Add validate_program post-pass

Walks the chumsky grammar's RawSyntaxNode tree and emits the public
RawInstruction tree, doing table lookups, arg-count checks, target
arity checks, and target usize parsing with annotation tolerance for
DETECTOR / OBSERVABLE_INCLUDE / etc. Mirrors the validation logic of
the existing parse_line; will replace it in the next commit.

Not yet wired into pub fn parse — exercised only by unit tests.
EOF
)"
```

---

### Task 9: Cutover — wire `parse()` to chumsky + validator and delete the hand-written parser

Replace the body of `parse()` to call `grammar::program_parser()` and `validate_program()`. Delete the now-unused hand-written code (`tokenize_lines`, `parse_block`, `parse_line`, `parse_head`, `split_head_and_targets`, `split_commas_shallow`, `parse_pi_expr`, `parse_tags`, `build_instruction`, `LineToken`, the `arity_of` if duplicated, etc.).

**Files:**
- Modify: `crates/ppvm-stim/src/parser/mod.rs` (rewrite `parse()`, delete legacy code)

- [ ] **Step 1: Rewrite `pub fn parse`**

Replace the existing `pub fn parse` body in `crates/ppvm-stim/src/parser/mod.rs` with:

```rust
/// Parse Stim source into a [`Program`].
pub fn parse(src: &str) -> Result<Program, ParseError> {
    let line_map = Arc::new(LineMap::new(src));
    let parse_result = grammar::program_parser().parse(src);
    let nodes = match parse_result.into_result() {
        Ok(nodes) => nodes,
        Err(errors) => {
            // Forward the first chumsky Rich error.
            let first = errors.into_iter().next().expect("non-empty error list");
            return Err(ParseError::Syntax(SyntaxError::from_rich(
                first,
                Arc::clone(&line_map),
            )));
        }
    };
    let instructions = validate_program(nodes, &line_map)?;
    Ok(Program { instructions })
}
```

- [ ] **Step 2: Delete hand-written parser code**

Delete every one of the following items from `crates/ppvm-stim/src/parser/mod.rs`. Remove their bodies *and* their `use` imports if those imports were only used by deleted code.

Delete:

- `parse_pi_expr`
- `split_commas_shallow`
- `parse_tags`
- `tokenize_lines`
- `LineToken` struct
- `parse_block`
- `parse_line`
- `parse_head`
- `split_head_and_targets`
- `build_instruction` (the old one — keep `build_instruction_v2`)

Keep:

- `LineMap` (with `pub` and the new `line_col` / `starts_at` methods)
- `RawSyntaxNode`, `RawTag`, `RawTagParam`, `RawTarget`
- `arity_of`
- `validate_program`, `validate_node`, `build_instruction_v2`
- `syntax` helper
- `pub fn parse`
- `mod grammar;` declaration

Rename `build_instruction_v2` back to `build_instruction` (the old one is gone now).

The legacy `syntax` helper signature stays as `fn syntax(line, col, message, &line_map) -> ParseError`; this is still used by `validate_node` for invalid-target errors.

- [ ] **Step 3: Run all `ppvm-stim` tests**

Run: `cargo test -p ppvm-stim`
Expected: every test PASSES — the existing parser_*.rs tests are now exercising the chumsky parser.

If any test fails, the failure pinpoints a behavioral divergence. Common ones to expect and how to fix:

- **`parser_syntax::comments_and_blank_lines_skipped`** — checks `line == 4` for the first instruction. The grammar must skip leading blank/comment lines via `pad` before the first item. Verify `pad()` is called via `pad().ignore_then(item)` inside the `recursive` body.
- **`parser_syntax::leading_indentation_tolerated`** — checks `"    H 0\n\tH 1"` parses two instructions. `pad` includes leading whitespace, so this should work.
- **`parser_syntax::parse_repeat_one_line`** — `REPEAT 5 { H 0 }`. `inline_pad` doesn't consume `\n` but does consume spaces; the grammar handles single-line REPEAT because `pad()` (with newlines) wraps the body content.
- **`parser_errors::at_least_one_target_required_for_h`** — `parse("H")`. Grammar parses `H` with zero targets. Validator checks `TargetArity::AtLeastOne` and returns `TargetCount { divisor: 1, found: 0, .. }`. Test allows either `TargetCount` or `Syntax`; both are acceptable.

- [ ] **Step 4: Run the full workspace test suite**

Run: `cargo test --workspace`
Expected: all PASS.

- [ ] **Step 5: Format and commit**

```bash
cd /Users/david/git/claude-repos/ppvm
cargo fmt --all
git add crates/ppvm-stim/src/parser/mod.rs
git commit -m "$(cat <<'EOF'
Wire chumsky grammar into pub fn parse and delete legacy parser

Replaces the hand-written line-based parser with the chumsky 0.12
grammar in parser/grammar.rs plus the validate_program post-pass.
The public API (parse, ParseError, Program) is unchanged; the
ParseError::Syntax variant now carries chumsky's Rich error verbatim
for genuine syntax failures.

Deletes tokenize_lines, parse_block, parse_line, parse_head,
split_head_and_targets, split_commas_shallow, parse_pi_expr,
parse_tags, and the old build_instruction.
EOF
)"
```

---

### Task 10: Clippy + fmt sweep

**Files:** none (verification only)

- [ ] **Step 1: Run clippy on the workspace**

Run: `cargo clippy --workspace --all-targets -- -D warnings`
Expected: no warnings.

If clippy flags issues introduced by the rewrite, fix them. Typical patterns:
- Unnecessary clone of `Arc<LineMap>` — replace with `Arc::clone(&x)` or restructure to borrow.
- Dead-code warnings on `RawTagParam::Positional` if a test inadvertently never constructs it — confirm the unit tests cover it.
- `needless_collect` on `.repeated().collect()` chains — chumsky's `collect` is idiomatic; if clippy complains, consider `#[allow(clippy::collect_into_iter)]` locally or refactor.

- [ ] **Step 2: Run fmt check**

Run: `cargo fmt --all -- --check`
Expected: clean (no diff).

- [ ] **Step 3: If clippy fixes were applied, format and commit**

```bash
cd /Users/david/git/claude-repos/ppvm
cargo fmt --all
git add -u
git commit -m "Address clippy warnings from chumsky parser rewrite"
```

(Skip this step if there were no clippy fixes.)

---

### Task 11: Python verification

Rebuild the native extension and run pytest.

**Files:** none (verification only — no commits)

- [ ] **Step 1: Reinstall the Python package against the new Rust code**

Run: `cd ppvm-python && uv sync --reinstall`
Expected: builds clean. The native extension is rebuilt against the chumsky-based parser.

- [ ] **Step 2: Run pytest**

Run: `cd ppvm-python && uv run pytest`
Expected: all PASS.

If a test fails, investigate. Likely causes are:
- Stim source string that the new grammar handles differently. Check the failing input against the grammar.
- A Python wrapper relying on the old `ParseError` Display format. Inspect with `pytest -s` and look at the message.

Fix any genuine bugs in `parser/grammar.rs` or `parser/mod.rs`. Add a Rust unit test in `tests/parser_*.rs` that pins the bug for the future. Re-run from Task 11 Step 1.

If fixes were needed:

```bash
cd /Users/david/git/claude-repos/ppvm
cargo fmt --all
git add -u
git commit -m "Fix chumsky parser regression discovered by Python tests"
```

(No commit if pytest was clean on first try.)

---

### Task 12: Extended-tests bug-finding sweep

Run the broader `.stim` corpus from `extended_tests/stim/run_all_programs.py`. This is a bug-finding net — any output artifacts must **not** be committed.

**Files:** none (verification only — no commits unless bugs are found and fixed)

- [ ] **Step 1: Run the extended tests**

Run: `cd ppvm-python && uv run python extended_tests/stim/run_all_programs.py`
Expected: completes without unexpected failures relative to the pre-rewrite baseline. (If the script does not have a clean baseline, capture stderr/stdout and review for new failure modes.)

- [ ] **Step 2: Check git status**

Run: `git status`
Expected: clean — no new files or modifications. If anything was generated, review and discard with `git clean -fd <path>` *only* on confirmed artifacts (do not blanket-clean).

- [ ] **Step 3: Investigate and fix any genuine bugs**

For each new failure (relative to the pre-rewrite baseline):
- Reproduce in isolation by feeding the offending `.stim` snippet through `parse()` in a minimal Rust test.
- Add the test to `crates/ppvm-stim/tests/parser_*.rs` (using the most fitting file).
- Fix the grammar or validator.
- Re-run Tasks 11 and 12.

If fixes are needed, commit each fix with a tight message:

```bash
cd /Users/david/git/claude-repos/ppvm
cargo fmt --all
git add -u
git commit -m "Fix <one-line bug summary> found by extended-tests sweep"
```

---

### Task 13: Final acceptance check

**Files:** none (verification only)

- [ ] **Step 1: Run the full acceptance matrix**

```bash
cd /Users/david/git/claude-repos/ppvm
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cd ppvm-python && uv run pytest
```

Expected: every command exits 0.

- [ ] **Step 2: Confirm legacy parser code is gone**

Run: `grep -E "tokenize_lines|parse_block|parse_line|parse_head|split_head_and_targets|split_commas_shallow|parse_pi_expr|parse_tags" crates/ppvm-stim/src/parser/mod.rs`
Expected: no matches.

- [ ] **Step 3: Confirm test corpus untouched**

Run: `git log --oneline crates/ppvm-stim/tests/data/ crates/ppvm-stim/tests/regen-stim/ crates/ppvm-stim/tests/stim_corpus.rs origin/david/ppvm-stim-2..HEAD`
Expected: no commits in this PR touch these paths.

- [ ] **Step 4: Confirm branch is the same one we started on**

Run: `git rev-parse --abbrev-ref HEAD`
Expected: `david/ppvm-stim-2`.

---

## Self-Review

**Spec coverage:** Decisions D1-D12 are covered: typed validation variants stay in `validate_program`; `SyntaxError` carries `Rich<'static, char>` plus `Arc<LineMap>`; grammar lives in `parser/grammar.rs`; validation stays as a post-pass; `LineMap` is public; target parsing is raw-lexeme-first; corpus paths are explicitly protected; final checks confirm branch and corpus history.

**Placeholder scan:** No red-flag placeholder markers remain. Verification-only tasks describe exact commands and expected outcomes.

**Type consistency:** `RawSyntaxNode`, `RawTag`, `RawTagParam`, `RawTarget`, `SyntaxError`, `LineMap::line_col`, `LineMap::starts_at`, and `validate_program` names are consistent across tasks. `SimpleSpan` examples use `SimpleSpan::from(start..end)`, matching chumsky 0.12.

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-04-30-ppvm-stim-chumsky-parser.md`. Two execution options:

**1. Subagent-Driven (recommended)** - Dispatch a fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** - Execute tasks in this session using executing-plans, batch execution with checkpoints.
