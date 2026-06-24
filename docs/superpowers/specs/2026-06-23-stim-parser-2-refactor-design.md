# Design: `stim-parser-2` — typestate lowering pipeline with a diagnostic-sink effect model

- **Date:** 2026-06-23
- **Status:** Approved (design); pending implementation plan
- **Scope:** A ground-up reimplementation of the `stim-parser` crate as a new crate
  `stim-parser-2`, built in parallel with the existing crate and swapped in once
  behavioral feature-parity is verified. No API-compatibility constraints (the
  package is pre-release), so we are free to break the API wherever it improves
  user experience or architectural health.

## 1. Motivation

The current `stim-parser` crate works and is well-tested, but the architectural
review surfaced concentrated duplication and uneven structure:

1. **`RawPassthrough` mirrors `RawInstruction`.** A near-clone of four
   instruction variants plus a ~50-line `into_raw()` and a twin printer
   (`fmt_raw` vs `fmt_raw_passthrough`), all to encode one invariant.
2. **Instruction identity is defined in three drift-prone places** — the name
   enums, the `canonical_name()` match (variant → string), and `TABLE` (string →
   variant + arity). Nothing checks they agree.
3. **`ast.rs` is a grab-bag** mixing AST nodes, name enums, table-rule enums, and
   error types.
4. **Uneven diagnostics** — only `ParseError::Syntax` reports `line:col`; typed
   variants report `line` only. Three fragmented error types
   (`SyntaxError`/`ParseError`/`ExtendedParseError`).
5. **Pretty-printing is internal-only** (behind `Display`) and normalizing.

This redesign keeps the proven pipeline shape (syntax → validate → lower) and the
chumsky syntax engine + its proptest safety net, while restructuring around four
goals the user specified:

- **Object-oriented driver** expressed as an explicit **state machine**.
- **Effects/errors as an effect algebra** handled as **continuations**.
- **Two public AST enums** (vanilla Stim + extended dialect).
- **First-class canonical pretty-printer** enabling round-trips on both layers.

## 2. Locked decisions

| Decision | Choice |
|---|---|
| Delivery | New crate `stim-parser-2`, built in parallel; swap + rename after parity |
| State machine | **Typestate** — each stage a distinct type; transitions consume `self` |
| Effect model | **Diagnostic sink / handler** — the handler is the continuation |
| AST shape | **Two public enums**, sharing per-family payload structs |
| Instruction metadata | **One declarative `const TABLE`** as the single source of truth |
| Pretty-printer | **Canonical / normalizing**, promoted to first-class public API |
| Spans | Each node carries a `Span`; programs own the `LineMap` |
| Syntax engine | **Keep chumsky** + the oversized-stack thread |
| API compatibility | **None required** (pre-release); free to break |

## 3. Architecture

### 3.1 Pipeline (data flow)

```
&str
 │  syntax/        chumsky combinators — PURE SYNTAX, zero semantics
 ▼
RawSyntaxTree      RawSyntaxNode { name:String, tags, args, targets:RawTarget{text,span} } | Repeat
 │  pipeline/validate.rs   table-driven: name→variant, arity checks, target typing, MPP/rec parse
 ▼
Program            two public enum: vanilla Stim AST; tags preserved verbatim
 │  pipeline/lower.rs      tag promotion: S[T]→T, I[R_X(..)]→Rotation, I_ERROR[loss]→Loss …
 ▼
ExtendedProgram    two public enum: typed PPVM dialect
```

The whole drive runs on the oversized parser-stack thread (`run_on_parser_stack`),
preserved from today, with the same `wasm32` inline carve-out.

### 3.2 Crate layout

```
crates/stim-parser-2/
  Cargo.toml
  src/
    lib.rs            public API + prelude; re-exports
    syntax/           Stage 1 — pure syntax (chumsky), no semantics
      mod.rs
      grammar.rs        combinators → RawSyntaxTree
      raw.rs            RawSyntaxNode / RawTarget
    instructions/     single source of truth: spelling ↔ variant ↔ arity ↔ canonical
      mod.rs            name enums (GateName/NoiseName/MeasureName/AnnotationKind),
                        arity enums (ArgCount/TargetArity), const TABLE,
                        lookup() + canonical_name() + completeness tests
    ast/
      mod.rs
      shared.rs         GateOp/NoiseOp/MeasureOp/AnnotationOp/MppOp, Target, PauliFactor, Axis,
                        Tag/TagParam (name enums are imported from instructions/)
      vanilla.rs        Instruction enum + Program
      extended.rs       ExtendedInstruction enum + ExtendedProgram
    pipeline/           Stage driver — the typestate state machine
      mod.rs            Pipeline<State>, transitions: parse / validate / lower / finish
      validate.rs       RawSyntaxTree → Program
      lower.rs          Program → ExtendedProgram (tag promotion)
    diagnostics/        the effect algebra
      mod.rs            Diagnostic, Severity, Flow, DiagnosticSink trait, Diagnostics
      sinks.rs          FailFast, Collect (provided handlers)
      span.rs           Span + LineMap + line:col rendering
    print/              first-class canonical printer
      mod.rs            Printer, PrintOptions, StimPrint trait
  tests/                ported integration tests + proptests + differential parity harness
```

Each module has one clear job and can be understood/tested independently.

## 4. Public API surface

### 4.1 Two tiers

```rust
// Tier 1 — the 99% case. Drives the pipeline internally with a FailFast sink.
pub fn parse(src: &str) -> Result<Program, Diagnostics>;
pub fn parse_extended(src: &str) -> Result<ExtendedProgram, Diagnostics>;

// Tier 2 — explicit typestate machine + caller-supplied sink.
let prog = Pipeline::new(src)        // Pipeline<Source>
    .parse(&mut sink)?               // Pipeline<Parsed>
    .validate(&mut sink)?            // Pipeline<Validated>   → .finish() → vanilla Program
    .lower(&mut sink)?               // Pipeline<Lowered>
    .finish();                       // ExtendedProgram
```

### 4.2 Typestate pipeline

```rust
pub struct Pipeline<S> { state: S }

pub struct Source<'a>   { src: &'a str }
pub struct Parsed       { tree: RawSyntaxTree,     line_map: Arc<LineMap> }
pub struct Validated    { program: Program,        line_map: Arc<LineMap> }
pub struct Lowered      { program: ExtendedProgram, line_map: Arc<LineMap> }

impl<'a> Pipeline<Source<'a>> {
    pub fn new(src: &'a str) -> Self;
    pub fn parse(self, sink: &mut dyn DiagnosticSink) -> Result<Pipeline<Parsed>, Aborted>;
}
impl Pipeline<Parsed>    { pub fn validate(self, sink: &mut dyn DiagnosticSink) -> Result<Pipeline<Validated>, Aborted>; }
impl Pipeline<Validated> { pub fn lower(self, sink: &mut dyn DiagnosticSink) -> Result<Pipeline<Lowered>, Aborted>;
                           pub fn finish(self) -> Program; }
impl Pipeline<Lowered>   { pub fn finish(self) -> ExtendedProgram; }
```

Each transition consumes `self` and returns the next state type, so illegal
orderings do not compile.

## 5. Effect model — the sink is the continuation

```rust
pub enum Severity { Error, Warning }
pub enum Flow { Continue, Abort }                 // the handler's continuation choice

pub struct Diagnostic {
    pub severity: Severity,
    pub span: Span,
    pub message: String,
    pub code: Option<&'static str>,               // stable kind tag, e.g. "unknown-instruction"
}

pub trait DiagnosticSink {
    /// Emit a diagnostic. The returned Flow tells the current stage whether to
    /// keep going (recover) or abort as soon as possible.
    fn emit(&mut self, d: Diagnostic) -> Flow;
}

pub struct FailFast(/* … */);   // returns Abort on first Error  → used by Tier-1 functions
pub struct Collect(/* … */);    // returns Continue, accumulates every diagnostic
```

- The handler's `Flow` return value **is** the continuation decision — this is the
  "effects handled as continuations" model encoded idiomatically in Rust.
- A transition told `Abort` returns `Err(Aborted)` (a tiny marker). The diagnostics
  themselves live in the **caller-owned sink**, so the pipeline stays policy-free.
- chumsky's multi-error output is forwarded into the sink, so even syntax errors can
  be collected (with `Collect`) rather than first-wins.
- `Diagnostics(Vec<Diagnostic>)` is the aggregate returned by the Tier-1 functions on
  failure; its `Display` prints every diagnostic, one per line, each with `line:col`.

## 6. AST model

### 6.1 Shared family payloads (eliminates `RawPassthrough`)

Both enums embed the *same* payload structs for the families that do not diverge:

```rust
// ast/shared.rs — reused by BOTH enums; printed once.
pub struct GateOp       { pub name: GateName,     pub tags: Vec<Tag>, pub args: Vec<f64>, pub targets: Vec<Target>, pub span: Span }
pub struct NoiseOp      { pub name: NoiseName,    pub tags: Vec<Tag>, pub args: Vec<f64>, pub targets: Vec<usize>,  pub span: Span }
pub struct MeasureOp    { pub name: MeasureName,  pub tags: Vec<Tag>, pub args: Vec<f64>, pub targets: Vec<usize>,  pub span: Span }
pub struct AnnotationOp { pub kind: AnnotationKind,                   pub args: Vec<f64>, pub targets: Vec<usize>,  pub span: Span }
pub struct MppOp        { pub tags: Vec<Tag>, pub args: Vec<f64>, pub products: Vec<Vec<PauliFactor>>, pub span: Span }
```

### 6.2 The two enums

```rust
// ast/vanilla.rs
pub enum Instruction {
    Gate(GateOp), Noise(NoiseOp), Measure(MeasureOp), Annotation(AnnotationOp),
    Mpp(MppOp),
    MPad { tags: Vec<Tag>, prob: Option<f64>, bits: Vec<usize>, span: Span },   // raw bits
    Repeat { count: u64, body: Vec<Instruction>, span: Span },
}
pub struct Program { pub instructions: Vec<Instruction>, /* + Arc<LineMap> */ }

// ast/extended.rs
pub enum ExtendedInstruction {
    Gate(GateOp), Noise(NoiseOp), Measure(MeasureOp), Annotation(AnnotationOp),  // SAME structs
    Mpp(MppOp),
    // promoted sugar (genuinely divergent — stays distinct):
    T { targets: Vec<usize>, span: Span },
    TDag { targets: Vec<usize>, span: Span },
    Rotation { axis: Axis, theta: f64, targets: Vec<usize>, span: Span },
    U3 { theta: f64, phi: f64, lambda: f64, targets: Vec<usize>, span: Span },
    Loss { p: f64, targets: Vec<usize>, span: Span },
    CorrelatedLoss { ps: [f64; 3], targets: Vec<(usize, usize)>, span: Span },
    MPad { tags: Vec<Tag>, prob: Option<f64>, bits: Vec<bool>, span: Span },     // validated bits
    Repeat { count: u64, body: Vec<ExtendedInstruction>, span: Span },
}
pub struct ExtendedProgram { pub instructions: Vec<ExtendedInstruction>, /* + Arc<LineMap> */ }
```

The only divergences between the layers are MPad's bit type (`usize` vs `bool`) and
the promoted-sugar variants — everything else is literally the same struct. Result:
no `RawPassthrough`, no `into_raw()`, and the printer formats each family op once.

`ExtendedProgram::measurement_count()` is preserved (pure AST property,
backend-agnostic).

### 6.3 Spans

- Each node carries a `Span { start: usize, end: usize }` (byte range).
- `Program` / `ExtendedProgram` own the `Arc<LineMap>`.
- Helper `node_span.line_col(&line_map) -> (usize, usize)` resolves on demand,
  giving uniform `line:col` in every diagnostic.

## 7. Single-source instruction table

One declarative `const TABLE` is the only place an instruction is defined; it
encodes the three relationships the enum cannot:

```rust
// row = (spelling, family+variant, args rule, target rule, canonical spelling)
("CNOT", Gate(CNot), NoArgs, Pairs, "CNOT"),
("CX",   Gate(CX),   NoArgs, Pairs, "CX"),
("ZCX",  Gate(ZCX),  NoArgs, Pairs, "ZCX"),
// …
```

- `lookup(name: &str) -> Option<TableEntry>` — scan by spelling (for parsing).
- `canonical_name(variant) -> &'static str` — derived from the same rows (for printing).
- Arity columns (`ArgCount` / `TargetArity`) drive validation.

**Drift-proofing tests:** every enum variant has exactly one row; every spelling is
unique; `lookup(canonical_name(v)).variant == v` for every variant. `I_ERROR`
remains `ArgCount::Deferred` (tag-specific arg rules enforced in `lower`).

## 8. Pretty printer (first-class, canonical)

```rust
pub struct PrintOptions { pub indent: Cow<'static, str> }   // default "    "
pub trait StimPrint {
    fn print(&self, out: &mut impl fmt::Write, opts: &PrintOptions, depth: usize) -> fmt::Result;
}

impl Program         { pub fn to_stim(&self) -> String; pub fn to_stim_with(&self, opts: &PrintOptions) -> String; }
impl ExtendedProgram { /* same */ }
impl fmt::Display for Program;          // delegates to to_stim() with defaults
impl fmt::Display for ExtendedProgram;  //
```

`StimPrint` is implemented on the shared family payloads so `GateOp`/`NoiseOp`/… are
formatted exactly once (no `fmt_raw` / `fmt_raw_passthrough` twins). Behavior matches
today's canonical printer: 4-space indent, `FloatLit` always emits a decimal point,
`rec[-k]` and `X0*Y1` targets round-trip, comments dropped. Round-trip is a *semantic*
fixpoint (`parse → print → parse`), not byte-identical.

## 9. Validation responsibility (unchanged contract)

The split with `ppvm-stim` is preserved:

- **`stim-parser-2`** checks dialect-level well-formedness: arity, target shape,
  `rec[-k]` only on gates, MPP factor syntax, MPAD bit values.
- **`ppvm-stim`** checks capability/semantics: supported gates/noise/measures,
  probability ranges, `rec[-k]` only on controlled Paulis, MPP distinct-qubit,
  record range.

## 10. Testing, parity & swap

### 10.1 Build sequence (each step independently testable)

1. `diagnostics/` — Span, LineMap, sink, `FailFast`/`Collect`.
2. `instructions/` — the table + completeness tests.
3. `ast/` — shared payloads + two enums + `measurement_count`.
4. `syntax/` — chumsky grammar → `RawSyntaxTree`.
5. `pipeline/` — typestate transitions parse/validate/lower + Tier-1 `parse`/`parse_extended`.
6. `print/` — `StimPrint` + `Display`.
7. Port `tests/*.rs` and all three proptest suites (roundtrip, ast, parse).
8. Differential parity harness.
9. Swap.

### 10.2 Differential parity harness (the no-regression gate)

A test target that dev-depends on **both** crates. For a corpus (the existing test
programs + proptest-generated programs) it runs old `stim_parser::parse_extended`
and new `stim_parser_2::parse_extended` and asserts two **operationally checkable**
properties (the crates have different AST types, so we compare observable behavior,
not the structs directly):

1. **Same accept/reject decision** — both succeed or both fail on each input.
2. **Byte-identical canonical print output** — old `format!("{}", prog)` (its
   `Display`) equals new `prog.to_stim()`. The new canonical format must therefore
   reproduce the old printer's output exactly on the corpus.

Structural equivalence beyond what printing surfaces (e.g. spans, internal field
shapes) is covered by the **ported unit/proptest suites**, not the differential
harness. Together these prove feature parity before any consumer is touched.

### 10.3 Swap phase (only after parity is green)

1. Repoint `ppvm-stim` at `stim-parser-2`. Mechanical edits:
   - `RawPassthrough::{Gate,Noise,Measure,Annotation}` → the shared `*Op` structs.
   - `*line` field reads → `op.span.line_col(&line_map)` (program owns the `LineMap`).
   - `ParseError`/`ExtendedParseError` → `Diagnostics`.
2. `cargo test --workspace` green.
3. Delete old `stim-parser`; rename `stim-parser-2` → `stim-parser`
   (crate dir, `Cargo.toml` `name`, `ppvm-stim/Cargo.toml` dependency); retire the
   differential harness.

## 11. Non-goals (YAGNI)

- **No lossless/trivia-preserving round-trip.** Canonical normalizing printer only.
  The AST/state-machine shape leaves room to add a trivia layer later, but it is out
  of scope here.
- **No hand-written parser.** Keep chumsky for the syntax stage; a
  dependency-free rewrite (and dropping the stack-thread workaround) is a possible
  follow-up once parity is locked, not part of this work.
- **No streaming/iterator API.** Programs materialize fully, as today.
- **No new accent on runtime dispatch in the AST.** The two enums stay enums (for
  exhaustive matching in `ppvm-stim`); "object-oriented" applies to the pipeline
  driver and the table/printer methods, not to trait-object AST nodes.
```
