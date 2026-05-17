---
title: ppvm-stim chumsky Parser Rewrite Design
date: 2026-04-30
branch: david/ppvm-stim-2
---

# ppvm-stim chumsky Parser Rewrite Design

## Summary

Replace the hand-written line-based Stim parser in `crates/ppvm-stim/src/parser/mod.rs` with a chumsky 0.12 grammar. The original ppvm-stim spec (`2026-04-21-ppvm-stim-design.md`) called for chumsky from the start, but the implementation plan deviated to a hand-written parser to ship sooner; this spec captures the agreed-upon swap.

The rewrite is a drop-in replacement for the public API: `parse(&str) -> Result<Program, ParseError>` keeps the same signature and the same `Program` AST. The change surfaces are:

- `ParseError::Syntax { line, col, message }` (struct variant) becomes `ParseError::Syntax(SyntaxError)` (tuple variant) where `SyntaxError` carries chumsky's `Rich<'static, char>` plus an `Arc<LineMap>` for `Display` line/col formatting.
- The other typed validation variants (`UnknownInstruction`, `ArgCount`, `TargetCount`) are unchanged.
- Existing tests that match `ParseError::Syntax { … }` need their pattern updated to `ParseError::Syntax(_)`. All other tests pass without modification.

## Context

### Current state

`crates/ppvm-stim/src/parser/mod.rs` is a 515-line hand-written parser:

- Per-line tokenizer (`tokenize_lines`) that emits one `LineToken` per non-comment, non-blank line, plus standalone `{`/`}` tokens for REPEAT block boundaries.
- Recursive `parse_block` for top-level + REPEAT bodies.
- Per-line `parse_line` that splits head from targets, parses tags via `parse_head` + `parse_tags`, parses pi-expressions via `parse_pi_expr`, and runs table-driven validation (instruction lookup, arg count, target arity).
- A `LineMap` exists but is dead-coded (line numbers are tracked through tokenization).

`chumsky = "0.12"` is already in `Cargo.toml` dependencies.

### Why rewrite

- The original ppvm-stim design spec called for chumsky 0.12.
- chumsky-based grammars are easier to extend with structured spans, multi-token error context, and richer error reporting once we lift the parser into its own crate (a goal already noted in the original spec).
- Today's parser hardcodes `col: 1` in most syntax errors; chumsky gives real spans for free.

## Decisions

| # | Decision |
|---|---|
| D1 | Keep typed validation variants (`UnknownInstruction`, `ArgCount`, `TargetCount`). Do not collapse them into chumsky errors. |
| D2 | Replace `ParseError::Syntax { line, col, message }` with `ParseError::Syntax(SyntaxError)` carrying `Rich<'static, char>` + `Arc<LineMap>`. |
| D3 | `Display` of `SyntaxError` formats as `syntax error at line {l}, col {c}: {rich}` using the `LineMap` to convert byte spans. |
| D4 | Validation is a separate post-pass over the grammar's raw tree. Do **not** drive it from chumsky `try_map_with` / custom `chumsky::error::Error` impl. |
| D5 | Grammar lives in `parser/grammar.rs` (combinators only, top-to-bottom). Public `parse`, `LineMap`, intermediate raw types, and validator live in `parser/mod.rs`. |
| D6 | Single-shot rewrite. No side-by-side parser / parity harness. Existing tests are the safety net. |
| D7 | Newlines are significant in the grammar — each instruction is a logical line, matching today's semantics and the Stim format. Comments consume to EOL. |
| D8 | `LineMap` lives in `parser/mod.rs` as a `pub` type. `ast.rs` imports it via `super::LineMap` for use in `SyntaxError`. |
| D9 | Targets are parsed as raw lexemes by the grammar; the validator does the `usize` parse and the annotation tolerance for non-numeric tokens (e.g. `rec[-1]`). |
| D10 | `cargo fmt --all` runs before every commit. |
| D11 | Test corpus (`tests/data/`, `tests/regen-stim/`, `stim_corpus.rs`) is untouched. |
| D12 | Work on the existing `david/ppvm-stim-2` branch. |

## Architecture

The pipeline shape stays the same:

```
&str ── parse() ──▶ Program ── normalize ──▶ TableauProgram ── execute ──▶ shots
```

Internal structure of `parse(src)`:

```
                ┌──────────────┐
                │ build LineMap│
                └──────┬───────┘
                       │
                       ▼
          ┌────────────────────────┐
          │  chumsky grammar       │  grammar.rs
          │  src ─▶ Vec<RawNode>   │  pure syntax,
          │  errors: Vec<Rich>     │  Rich<char> errors
          └──────┬─────────────────┘
                 │ ok                ▼ err
                 ▼              SyntaxError(rich, line_map)
       ┌─────────────────┐
       │ validate_program│  mod.rs
       │ table lookup,   │  emits typed
       │ arg/target check│  validation errors
       │ Vec<RawNode>    │
       │   ─▶ Program    │
       └────────┬────────┘
                │
                ▼
             Program
```

### File layout (`crates/ppvm-stim/src/parser/`)

```
parser/
├── mod.rs       # pub fn parse, LineMap, RawSyntaxNode types, validator
├── grammar.rs   # chumsky combinators only, top-to-bottom (NEW)
├── ast.rs       # AST + ParseError (modified) + SyntaxError (NEW) + Display impls
└── table.rs     # unchanged
```

### Module responsibilities

**`parser/mod.rs`** — public entry point + glue.

- `LineMap` type. Builds the byte-offset → line index. Adds `line_col(byte: usize) -> (usize, usize)`.
- `RawSyntaxNode` — internal `pub(crate)` enum representing the grammar's output before validation. One variant per top-level shape: `Instruction { name, tags, args, targets, span }`, `Repeat { count, body, span }`. `RawTag`, `RawTagParam`, `RawTarget` (which is just `String` + `SimpleSpan`) live here too.
- `pub fn parse(src: &str) -> Result<Program, ParseError>` — builds `Arc<LineMap>`, runs the chumsky parser via `grammar::program_parser()`, on grammar failure wraps the first `Rich` into `ParseError::Syntax`, on success calls `validate_program(raw, &line_map) -> Result<Program, ParseError>`.
- `validate_program` — walks `Vec<RawSyntaxNode>` and emits typed `RawInstruction`s. Drives table lookup, arg-count check, target-arity check, target-`usize` parse with annotation tolerance. Recurses into REPEAT bodies.

**`parser/grammar.rs`** — chumsky combinators only.

- One `pub(crate) fn program_parser<'src>() -> impl Parser<'src, &'src str, Vec<RawSyntaxNode>, extra::Err<Rich<'src, char>>>`.
- Private helper combinators above it, top-to-bottom: `pad`, `signed_float`, `pi_expr`, `ident`, `tag_param`, `tag`, `tags_block`, `args_block`, `target_lexeme`, `instruction_head`, `instruction_line`, `repeat_block`, `program`.
- No knowledge of `parser::table`, no validation, no error reshaping. Pure syntax.

**`parser/ast.rs`** — public AST + error types.

- `Program`, `RawInstruction`, `Tag`, `TagParam`, `GateName`, `NoiseName`, `MeasureName`, `AnnotationKind`, `ArgCount`, `TargetArity` — unchanged.
- `ParseError`: the `Syntax` variant changes shape; others unchanged.
  ```rust
  #[derive(Debug, thiserror::Error)]
  #[non_exhaustive]
  pub enum ParseError {
      #[error("{0}")]
      Syntax(SyntaxError),

      #[error("unknown instruction '{name}' at line {line}")]
      UnknownInstruction { name: String, line: usize },

      #[error("'{name}' at line {line} expected {expected} args, got {found}")]
      ArgCount { name: String, expected: usize, found: usize, line: usize },

      #[error("'{name}' at line {line} expected target count divisible by {divisor}, got {found}")]
      TargetCount { name: String, divisor: usize, found: usize, line: usize },
  }
  ```
- `SyntaxError`:
  ```rust
  #[derive(Debug, Clone)]
  pub struct SyntaxError {
      pub(crate) rich: chumsky::error::Rich<'static, char>,
      pub(crate) line_map: Arc<LineMap>,
  }
  ```
- `Display for SyntaxError` — formats as `syntax error at line {l}, col {c}: {rich}` using `line_map.line_col(rich.span().start)`.
- `PartialEq` on `ParseError` is dropped (or implemented manually for the `Syntax` variant by comparing on stringified `rich` + spans). Existing tests rely on `matches!(...)` rather than `==`, so dropping the derive is safe; verify during implementation.

## Grammar shape

`grammar.rs` reads top-to-bottom:

1. **`pad`** — zero or more of: `' '`, `'\t'`, `'\r'`, `'\n'`, or `# … \n` (comment to EOL).
2. **`inline_pad`** — like `pad` but excludes `'\n'`. Used between target lexemes on a single line.
3. **`signed_float`** — `[+-]? digits ('.' digits)? ([eE] [+-]? digits)?`. Returns `f64`.
4. **`pi_expr`** — `("pi") | (signed_float "*pi") | signed_float`. Returns `f64`. Maps `pi` to `std::f64::consts::PI`, and `c*pi` to `c * std::f64::consts::PI`.
5. **`ident`** — `[A-Za-z_][A-Za-z0-9_]*`. Used for instruction names, tag names, tag-param keys. Case-distinction is the table's job.
6. **`tag_param`** — `ident "=" pi_expr` (Named) or `pi_expr` (Positional).
7. **`tag`** — `ident ("(" tag_param ("," tag_param)* ")")?`.
8. **`tags_block`** — `"[" tag ("," tag)* "]"`.
9. **`args_block`** — `"(" pi_expr ("," pi_expr)* ")"`.
10. **`target_lexeme`** — one or more non-whitespace, non-`#` chars. Returns `(String, SimpleSpan)`. The validator parses to `usize` later.
11. **`instruction_head`** — `ident tags_block? args_block?`.
12. **`instruction_line`** — `instruction_head (inline_pad target_lexeme)* &(newline | eof | "}")`. Targets are read until end-of-line; trailing `# comment` is consumed by `pad` between this line and the next.
13. **`repeat_block`** — `recursive(|body| "REPEAT" inline_pad digits inline_pad "{" pad (instruction_line | repeat_block)* pad "}")`.
14. **`program`** — `pad (instruction_line | repeat_block)* (pad)? then_ignore(end())`.

### Newline handling (Decision D7)

Newlines are significant only as instruction terminators. `instruction_line` reads targets greedily until it hits a newline, EOF, or `}`. This matches Stim's de facto line-oriented format and avoids the `H X 1` ambiguity (where `X` is both a valid gate and a valid lexeme). Comments and blank lines between instructions are consumed by `pad`.

## Validation pass

`fn validate_program(nodes: Vec<RawSyntaxNode>, line_map: &LineMap) -> Result<Vec<RawInstruction>, ParseError>`:

For each `RawSyntaxNode::Instruction { name, tags, args, targets, span }`:

1. **Lookup** — `table::lookup(&name)` → `TableEntry` or `ParseError::UnknownInstruction { name, line: line_map.line_of(span.start) }`.
2. **Tags** — convert `Vec<RawTag>` → `Vec<Tag>` (no validation; tag semantics live in `normalize`). `pi_expr` was already evaluated by the grammar.
3. **Args** — already `Vec<f64>` from the grammar. Apply the existing arg-count validation (with skip rules for annotations and `Noise { args: ArgCount::None }`).
4. **Targets** — for each `RawTarget(s, t_span)`:
   - `s.parse::<usize>()` → push to `targets: Vec<usize>`.
   - Parse fails + entry is `TableEntry::Annotation` → drop silently (matches today).
   - Parse fails + entry is non-annotation → `ParseError::Syntax(SyntaxError::synth(t_span, format!("invalid target {s:?}"), line_map.clone()))`.
5. **Target arity** — `Any | AtLeastOne | Pairs | Quadruples` checks → `ParseError::TargetCount`.
6. **Build** — dispatch on `TableEntry` into `RawInstruction::{Gate, Noise, Measure, Annotation}` exactly as today's `build_instruction`.

For `RawSyntaxNode::Repeat { count, body, span }` — recurse on `body`, build `RawInstruction::Repeat { count, body, line: line_map.line_of(span.start) }`.

`SyntaxError::synth(span, msg, line_map)` constructs a `Rich<'static, char>::custom(span, msg)` and wraps it. This is the only place the validator produces a `Syntax` variant — for "invalid target" tokens that the grammar accepted as raw lexemes but failed `usize` parsing.

## Error display

```rust
impl fmt::Display for SyntaxError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let span = self.rich.span();
        let (line, col) = self.line_map.line_col(span.start);
        write!(f, "syntax error at line {line}, col {col}: {}", self.rich)
    }
}
```

`Rich<'static, char>` is used because the input `&str` does not outlive `parse()`'s return value. `Rich::into_owned()` is called once when wrapping. The implementer should confirm this is the right call (chumsky 0.12's exact API); if `into_owned` is not the method name, the equivalent owned conversion suffices.

`LineMap::line_col(byte_offset)` returns `(line: usize, col: usize)` — both 1-indexed. Implementation uses the existing `starts: Vec<usize>` (one entry per line, byte offset of line start). `line` = `line_of(byte_offset)`; `col` = `byte_offset - starts[line - 1] + 1`.

## Testing

### Existing tests

| File | Status |
|---|---|
| `tests/parser_gates.rs` | Pass unchanged. |
| `tests/parser_noise.rs` | Pass unchanged. |
| `tests/parser_measure.rs` | Pass unchanged. |
| `tests/parser_tags.rs` | Pass unchanged. |
| `tests/parser_syntax.rs` | Pass unchanged. |
| `tests/parser_errors.rs` | Pattern updates only: `ParseError::Syntax { … }` → `ParseError::Syntax(_)`. Other variants unchanged. |
| `tests/normalize.rs` | Pass unchanged. |
| `tests/executor.rs` | Pass unchanged. |
| `tests/run.rs` | Pass unchanged. |
| `tests/stim_corpus.rs` | **Untouched per Decision D11.** |

### New tests

One new test pinning the new error-display behavior:

```rust
// tests/parser_errors.rs
#[test]
fn syntax_error_includes_line_and_col() {
    let err = parse("H 0\nFROBNICATE )").unwrap_err();
    let s = err.to_string();
    assert!(s.contains("line 2"), "{s}");
    assert!(s.contains("col"), "{s}");
}
```

(Plus the existing "unknown-instruction reports correct line" test in `parser_errors.rs` continues to cover line-number accuracy for the typed variants.)

### Verification matrix

The implementer runs all of these before declaring done. Outputs from the extended-tests script are **never committed** (per user instruction).

| Check | Command | When |
|---|---|---|
| Format | `cargo fmt --all` | Before every commit (D10). |
| Format check | `cargo fmt --all -- --check` | Final sweep. |
| Clippy | `cargo clippy --workspace --all-targets -- -D warnings` | Final sweep. |
| Rust unit + integration | `cargo test --workspace` | Per task. |
| Python rebuild | `cd ppvm-python && uv sync --reinstall` | After Rust changes are stable. |
| Python pytest | `cd ppvm-python && uv run pytest` | After rebuild. |
| Extended Python tests | `cd ppvm-python && uv run python extended_tests/stim/run_all_programs.py` | Bug-finding sweep before declaring done. Failures investigated, not silenced. Output artifacts not committed. |

If extended-test failures surface real bugs, fix them and consider distilling the failing input into a Rust unit test in `parser_*.rs` (with the offending `.stim` snippet inline, not as a corpus file).

## Out of scope

- Test corpus changes (`tests/data/`, `tests/regen-stim/`, `stim_corpus.rs`) — explicitly untouched.
- Richer error UX (multi-span, `ariadne`-style pretty printing) — possible future follow-up.
- Lifting the parser into its own crate — orthogonal future work; this rewrite makes that easier but doesn't do it.
- Changing the AST shape, the table, the normalizer, the executor, or the Python bindings.
- Driving validation from inside chumsky via custom `chumsky::error::Error` impl — explicitly rejected (D4).

## Risks

- **chumsky 0.12 API surface drift.** chumsky's API has shifted between minor versions; the implementer should write a tiny smoke parser first (e.g., `signed_float` alone, or just `ident`) to confirm `Rich`, `extra::Err`, `SimpleSpan`, `recursive`, and `then_ignore(end())` work as described before wiring the full grammar.
- **`Rich::into_owned` exact API.** If chumsky 0.12 spells this differently (e.g., `to_owned`, or requires a manual `Rich::map_token` step), the implementer adapts the `SyntaxError` construction site without changing the design.
- **Newline-significant grammar in chumsky.** Whitespace skipping that excludes `\n` in some contexts and includes it in others is a known footgun. The plan should include a focused test for `H 0\nX 0` parsing as two instructions, not one.
- **Annotation tolerance regressions.** `DETECTOR rec[-1]` and `OBSERVABLE_INCLUDE(0) rec[-3] rec[-1]` need to keep parsing (with the non-numeric tokens silently dropped). This is covered by the existing tests + the extended-test corpus, but it's the most fragile semantic to preserve.

## Acceptance

- `cargo test --workspace` green.
- `cargo clippy --workspace --all-targets -- -D warnings` green.
- `cargo fmt --all -- --check` green.
- `cd ppvm-python && uv sync --reinstall && uv run pytest` green.
- `cd ppvm-python && uv run python extended_tests/stim/run_all_programs.py` runs cleanly (no new regressions vs. pre-rewrite baseline).
- Old hand-written parser code (`tokenize_lines`, `parse_block`, `parse_line`, `parse_head`, `split_head_and_targets`, `split_commas_shallow`, `parse_pi_expr`, `parse_tags`, `arity_of`, `build_instruction`) is deleted, not commented out.
