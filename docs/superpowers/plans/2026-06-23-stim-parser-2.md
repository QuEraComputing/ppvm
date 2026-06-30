# `stim-parser-2` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reimplement `stim-parser` as a new crate `stim-parser-2` with a typestate lowering pipeline, a diagnostic-sink effect model, two AST enums sharing per-family payloads, a single-source instruction table, and a first-class canonical pretty-printer — then prove parity and swap it in.

**Architecture:** Three-stage pipeline (`syntax → validate → lower`) driven by a typestate `Pipeline<State>` whose transitions consume `self`. Each stage emits `Diagnostic`s to a caller-supplied `DiagnosticSink` whose `Flow` return value is the continuation decision. The chumsky syntax engine is retained; the AST is two enums (`Instruction`, `ExtendedInstruction`) sharing `GateOp`/`NoiseOp`/`MeasureOp`/`AnnotationOp`/`MppOp`.

**Tech Stack:** Rust 2024, chumsky 0.12, thiserror 2, proptest 1.

**Reference implementation:** The existing `crates/stim-parser/` is a working, well-tested reference that stays in the tree until the swap. Many tasks below are *port-and-adapt* from a named source file — open it, copy the logic, apply the stated transformation. This is intentional (DRY across crates) and not a placeholder.

**Design spec:** `docs/superpowers/specs/2026-06-23-stim-parser-2-refactor-design.md`.

---

## File Structure

New crate `crates/stim-parser-2/`:

| File | Responsibility |
|---|---|
| `Cargo.toml` | Crate manifest (name `stim-parser-2`, deps chumsky/thiserror, dev-dep proptest) |
| `src/lib.rs` | Module wiring, `prelude`, Tier-1 `parse`/`parse_extended` |
| `src/diagnostics/mod.rs` | `Severity`, `Flow`, `Diagnostic`, `DiagnosticSink`, `Aborted`, `Diagnostics` |
| `src/diagnostics/span.rs` | `Span`, `LineMap` |
| `src/diagnostics/sinks.rs` | `FailFast`, `Collect` |
| `src/instructions/mod.rs` | Name enums, arity enums, `TableEntry`/`EntryKind`, `TABLE`, `lookup`, `canonical_name`, completeness tests |
| `src/ast/mod.rs` | Re-exports |
| `src/ast/shared.rs` | `GateOp`/`NoiseOp`/`MeasureOp`/`AnnotationOp`/`MppOp`, `Target`, `PauliFactor`, `PauliAxis`, `Axis`, `Tag`, `TagParam` |
| `src/ast/vanilla.rs` | `Instruction`, `Program` |
| `src/ast/extended.rs` | `ExtendedInstruction`, `ExtendedProgram`, `measurement_count` |
| `src/syntax/mod.rs` | Re-exports + `run_on_parser_stack` |
| `src/syntax/raw.rs` | `RawSyntaxNode`, `RawTarget`, `RawSyntaxTree` |
| `src/syntax/grammar.rs` | chumsky combinators → `RawSyntaxTree` |
| `src/pipeline/mod.rs` | `Pipeline<State>`, state structs, transitions |
| `src/pipeline/validate.rs` | `RawSyntaxTree` → `Program` |
| `src/pipeline/lower.rs` | `Program` → `ExtendedProgram` |
| `src/print/mod.rs` | `PrintOptions`, `StimPrint`, `Display` impls, `to_stim` |
| `tests/*.rs` | Ported integration tests + proptests |
| `tests/parity.rs` | Differential harness vs old `stim-parser` |

Swap phase modifies: `Cargo.toml` (workspace members), `crates/ppvm-stim/Cargo.toml`, `crates/ppvm-stim/src/{lib,executor,validate}.rs`, `crates/ppvm-stim/benches/tableau-msd-stim.rs`.

---

## Phase 0 — Crate scaffold

### Task 0: Create the crate and register it in the workspace

**Files:**
- Create: `crates/stim-parser-2/Cargo.toml`
- Create: `crates/stim-parser-2/src/lib.rs`
- Modify: `Cargo.toml` (workspace `members`)

- [ ] **Step 1: Write `Cargo.toml`**

```toml
[package]
name = "stim-parser-2"
version = "0.1.0"
edition = "2024"

[dependencies]
chumsky = "0.12.0"
thiserror = "2.0.18"

[dev-dependencies]
proptest = "1"
stim-parser = { version = "0.1.0", path = "../stim-parser" }   # parity harness only; removed at swap
```

- [ ] **Step 2: Write a placeholder `src/lib.rs`**

```rust
// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0
```

- [ ] **Step 3: Register in workspace members**

In the root `Cargo.toml`, change the line `"crates/ppvm-stim", "crates/stim-parser",` to add the new crate:

```toml
    "crates/ppvm-stim", "crates/stim-parser", "crates/stim-parser-2",
```

- [ ] **Step 4: Verify it builds**

Run: `cargo build -p stim-parser-2`
Expected: `Finished` (empty crate compiles).

- [ ] **Step 5: Commit**

```bash
git add Cargo.toml crates/stim-parser-2/Cargo.toml crates/stim-parser-2/src/lib.rs
git commit -m "feat(stim-parser-2): scaffold new crate"
```

---

## Phase 1 — `diagnostics/` (Span, LineMap, sink, Diagnostics)

### Task 1: `Span` and `LineMap`

**Files:**
- Create: `crates/stim-parser-2/src/diagnostics/span.rs`
- Test: same file (`#[cfg(test)]`)
- Reference: `crates/stim-parser/src/line_map.rs`

- [ ] **Step 1: Write the failing tests**

```rust
// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

#[cfg(test)]
mod tests {
    use super::*;

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
    fn span_resolves_against_line_map() {
        let m = LineMap::new("X 0\nH 0");
        let span = Span::new(4, 5);
        assert_eq!(span.line_col(&m), (2, 1));
        assert_eq!(span.line(&m), 2);
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test -p stim-parser-2 --lib span`
Expected: FAIL — `LineMap`/`Span` not found.

- [ ] **Step 3: Implement `Span` and `LineMap`**

Port `LineMap` verbatim from `crates/stim-parser/src/line_map.rs` (the `new`/`line_of`/`line_col`/`starts_at` methods and the `Debug` impl are unchanged), and add `Span` above it:

```rust
// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Byte spans and the line/column map shared by every diagnostic.

/// Half-open byte range `[start, end)` into the source string.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}

impl Span {
    pub fn new(start: usize, end: usize) -> Self {
        Span { start, end }
    }

    /// 1-indexed `(line, col)` of the span start.
    pub fn line_col(&self, line_map: &LineMap) -> (usize, usize) {
        line_map.line_col(self.start)
    }

    /// 1-indexed line of the span start.
    pub fn line(&self, line_map: &LineMap) -> usize {
        line_map.line_of(self.start)
    }
}

impl From<chumsky::span::SimpleSpan<usize>> for Span {
    fn from(s: chumsky::span::SimpleSpan<usize>) -> Self {
        Span::new(s.start, s.end)
    }
}

// ---- LineMap: ported verbatim from crates/stim-parser/src/line_map.rs ----
```

Then paste the entire `LineMap` struct, its `Debug` impl, and its `impl LineMap { new, line_of, line_col, starts_at }` block from the reference file (drop that file's own `#[cfg(test)]` block — the tests live above).

- [ ] **Step 4: Wire the module**

Create `crates/stim-parser-2/src/diagnostics/mod.rs` with:

```rust
// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod span;
pub use span::{LineMap, Span};
```

And in `src/lib.rs` add:

```rust
pub mod diagnostics;
```

- [ ] **Step 5: Run to verify pass**

Run: `cargo test -p stim-parser-2 --lib span`
Expected: PASS (3 tests).

- [ ] **Step 6: Commit**

```bash
git add crates/stim-parser-2/src
git commit -m "feat(stim-parser-2): add Span and LineMap"
```

### Task 2: `Diagnostic`, `Severity`, `Flow`, `DiagnosticSink`, `Aborted`, `Diagnostics`

**Files:**
- Modify: `crates/stim-parser-2/src/diagnostics/mod.rs`
- Test: `crates/stim-parser-2/src/diagnostics/mod.rs` (`#[cfg(test)]`)

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn diagnostics_display_renders_line_col() {
        let line_map = Arc::new(LineMap::new("X 0\nBADINSTR 0"));
        let diag = Diagnostic::error(Span::new(4, 12), "unknown-instruction", "unknown instruction 'BADINSTR'");
        let diags = Diagnostics::new(vec![diag], line_map);
        assert_eq!(
            diags.to_string(),
            "error at line 2, col 1: unknown instruction 'BADINSTR'"
        );
    }

    #[test]
    fn diagnostics_is_empty_when_no_items() {
        let line_map = Arc::new(LineMap::new(""));
        assert!(Diagnostics::new(vec![], line_map).is_empty());
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test -p stim-parser-2 --lib diagnostics::tests`
Expected: FAIL — types not found.

- [ ] **Step 3: Implement the types**

Add to `src/diagnostics/mod.rs` (above the `mod span;` line keep module decls together; types below):

```rust
use std::fmt;
use std::sync::Arc;

/// Severity of a diagnostic. Only `Error` aborts a `FailFast` run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    Error,
    Warning,
}

/// A sink's continuation decision: keep processing the current stage, or
/// abort it as soon as possible. The handler returning `Flow` is how the
/// effect model "handles errors as continuations".
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Flow {
    Continue,
    Abort,
}

/// One diagnostic, carrying a span so every message can render `line:col`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnostic {
    pub severity: Severity,
    pub span: Span,
    pub message: String,
    /// Stable short kind tag (e.g. "unknown-instruction") for matching
    /// without string-sniffing.
    pub code: Option<&'static str>,
}

impl Diagnostic {
    pub fn error(span: Span, code: &'static str, message: impl Into<String>) -> Self {
        Diagnostic { severity: Severity::Error, span, message: message.into(), code: Some(code) }
    }

    pub fn warning(span: Span, code: &'static str, message: impl Into<String>) -> Self {
        Diagnostic { severity: Severity::Warning, span, message: message.into(), code: Some(code) }
    }
}

/// A handler the pipeline emits diagnostics to. The returned `Flow` tells
/// the emitting stage whether to continue (recover) or abort.
pub trait DiagnosticSink {
    fn emit(&mut self, diagnostic: Diagnostic) -> Flow;
}

/// Marker returned by a pipeline transition that was told to `Abort`.
/// The diagnostics themselves live in the caller-owned sink.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Aborted;

/// Aggregate returned by the Tier-1 `parse`/`parse_extended` functions on
/// failure. Owns a `LineMap` so `Display` can render `line:col`.
#[derive(Debug, Clone)]
pub struct Diagnostics {
    items: Vec<Diagnostic>,
    line_map: Arc<LineMap>,
}

impl Diagnostics {
    pub fn new(items: Vec<Diagnostic>, line_map: Arc<LineMap>) -> Self {
        Diagnostics { items, line_map }
    }

    pub fn iter(&self) -> impl Iterator<Item = &Diagnostic> {
        self.items.iter()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
}

impl fmt::Display for Diagnostics {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, d) in self.items.iter().enumerate() {
            if i > 0 {
                writeln!(f)?;
            }
            let (line, col) = d.span.line_col(&self.line_map);
            let sev = match d.severity {
                Severity::Error => "error",
                Severity::Warning => "warning",
            };
            write!(f, "{sev} at line {line}, col {col}: {}", d.message)?;
        }
        Ok(())
    }
}

impl std::error::Error for Diagnostics {}
```

- [ ] **Step 4: Re-export the new types**

Update the top of `src/diagnostics/mod.rs` so the public surface is:

```rust
mod span;
pub use span::{LineMap, Span};
```

(The new types are already `pub` at module root, so `pub mod diagnostics;` in `lib.rs` exposes them.)

- [ ] **Step 5: Run to verify pass**

Run: `cargo test -p stim-parser-2 --lib diagnostics::tests`
Expected: PASS (2 tests).

- [ ] **Step 6: Commit**

```bash
git add crates/stim-parser-2/src/diagnostics/mod.rs
git commit -m "feat(stim-parser-2): add Diagnostic, DiagnosticSink, Diagnostics"
```

### Task 3: `FailFast` and `Collect` sinks

**Files:**
- Create: `crates/stim-parser-2/src/diagnostics/sinks.rs`
- Modify: `crates/stim-parser-2/src/diagnostics/mod.rs`
- Test: `crates/stim-parser-2/src/diagnostics/sinks.rs` (`#[cfg(test)]`)

- [ ] **Step 1: Write the failing tests**

```rust
// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagnostics::{Diagnostic, DiagnosticSink, Flow, Severity, Span};

    fn err() -> Diagnostic { Diagnostic::error(Span::new(0, 1), "x", "boom") }
    fn warn() -> Diagnostic { Diagnostic::warning(Span::new(0, 1), "x", "heads up") }

    #[test]
    fn fail_fast_aborts_on_first_error_and_keeps_it() {
        let mut s = FailFast::new();
        assert_eq!(s.emit(warn()), Flow::Continue);   // warnings don't abort
        assert_eq!(s.emit(err()), Flow::Abort);
        let items = s.into_items();
        assert_eq!(items.len(), 2);
        assert_eq!(items[1].severity, Severity::Error);
    }

    #[test]
    fn collect_never_aborts_and_gathers_all() {
        let mut s = Collect::new();
        assert_eq!(s.emit(err()), Flow::Continue);
        assert_eq!(s.emit(err()), Flow::Continue);
        assert_eq!(s.into_items().len(), 2);
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test -p stim-parser-2 --lib sinks`
Expected: FAIL — `FailFast`/`Collect` not found.

- [ ] **Step 3: Implement the sinks**

```rust
// (license header above)

//! Two provided `DiagnosticSink` policies.

use crate::diagnostics::{Diagnostic, DiagnosticSink, Flow, Severity};

/// Aborts the current stage on the first `Error` (warnings pass through).
/// Used by the Tier-1 `parse`/`parse_extended` functions.
#[derive(Debug, Default)]
pub struct FailFast {
    items: Vec<Diagnostic>,
    saw_error: bool,
}

impl FailFast {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn saw_error(&self) -> bool {
        self.saw_error
    }

    pub fn into_items(self) -> Vec<Diagnostic> {
        self.items
    }
}

impl DiagnosticSink for FailFast {
    fn emit(&mut self, diagnostic: Diagnostic) -> Flow {
        let is_error = diagnostic.severity == Severity::Error;
        self.items.push(diagnostic);
        if is_error {
            self.saw_error = true;
            Flow::Abort
        } else {
            Flow::Continue
        }
    }
}

/// Never aborts; accumulates every diagnostic for one-pass reporting.
#[derive(Debug, Default)]
pub struct Collect {
    items: Vec<Diagnostic>,
}

impl Collect {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn into_items(self) -> Vec<Diagnostic> {
        self.items
    }
}

impl DiagnosticSink for Collect {
    fn emit(&mut self, diagnostic: Diagnostic) -> Flow {
        self.items.push(diagnostic);
        Flow::Continue
    }
}
```

- [ ] **Step 4: Wire the module**

In `src/diagnostics/mod.rs` add `mod sinks;` and `pub use sinks::{Collect, FailFast};`.

- [ ] **Step 5: Run to verify pass**

Run: `cargo test -p stim-parser-2 --lib sinks`
Expected: PASS (2 tests).

- [ ] **Step 6: Commit**

```bash
git add crates/stim-parser-2/src/diagnostics
git commit -m "feat(stim-parser-2): add FailFast and Collect sinks"
```

---

## Phase 2 — `instructions/` (single-source table)

### Task 4: Name enums, arity enums, table types

**Files:**
- Create: `crates/stim-parser-2/src/instructions/mod.rs`
- Modify: `crates/stim-parser-2/src/lib.rs`
- Reference: `crates/stim-parser/src/ast.rs` (enums + `canonical_name`), `crates/stim-parser/src/table.rs`

- [ ] **Step 1: Port the name enums and arity enums**

Into `src/instructions/mod.rs`, port verbatim from `crates/stim-parser/src/ast.rs`:
- `GateName`, `NoiseName`, `MeasureName`, `AnnotationKind` (the enum definitions only — *not* their `canonical_name` impls yet),
- `ArgCount`, `TargetArity`.

Then port from `crates/stim-parser/src/table.rs`: `TableEntry`, `EntryKind` (keep `MPad` variant), and the `gate`/`noise`/`measure`/`measure_pairs`/`annotation` const constructors. Add the license header. Do **not** port `canonical_name` or `TABLE` yet (next steps).

- [ ] **Step 2: Wire the module and build**

In `src/lib.rs` add `pub mod instructions;`.
Run: `cargo build -p stim-parser-2`
Expected: `Finished`.

- [ ] **Step 3: Commit**

```bash
git add crates/stim-parser-2/src
git commit -m "feat(stim-parser-2): add instruction name/arity enums and table types"
```

### Task 5: The `TABLE` with a canonical column, `lookup`, `canonical_name`

**Files:**
- Modify: `crates/stim-parser-2/src/instructions/mod.rs`
- Test: same file
- Reference: `crates/stim-parser/src/table.rs` (`TABLE`), `crates/stim-parser/src/ast.rs` (`canonical_name` impls)

- [ ] **Step 1: Add a canonical column to `TableEntry`**

Extend `TableEntry` with a `canonical: &'static str` field:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TableEntry {
    pub kind: EntryKind,
    pub args: ArgCount,
    pub targets: TargetArity,
    pub canonical: &'static str,
}
```

Update the const constructors to take and store it, e.g.:

```rust
const fn gate(name: GateName, args: ArgCount, targets: TargetArity, canonical: &'static str) -> TableEntry {
    TableEntry { kind: EntryKind::Gate(name), args, targets, canonical }
}
// …same shape for noise/measure/measure_pairs/annotation…
```

- [ ] **Step 2: Port `TABLE` with the canonical spelling per row**

Port the full `TABLE` from `crates/stim-parser/src/table.rs`, adding the canonical spelling as the last constructor arg. The canonical spelling for each row is the spelling produced by the old `canonical_name()` for that variant (open `crates/stim-parser/src/ast.rs` and read the matching arm). Examples (note the alias rows `R`/`RZ`, `CX`/`ZCX`/`CNOT`, `S`/`SQRT_Z` each carry their *own* canonical spelling, matching the old behavior where every variant is distinct):

```rust
const TABLE: &[(&str, TableEntry)] = &[
    ("R",   gate(G::Reset,  NoArgs, AtLeastOne, "R")),
    ("RZ",  gate(G::ResetZ, NoArgs, AtLeastOne, "RZ")),
    ("RX",  gate(G::ResetX, NoArgs, AtLeastOne, "RX")),
    ("RY",  gate(G::ResetY, NoArgs, AtLeastOne, "RY")),
    ("X",   gate(G::X, NoArgs, AtLeastOne, "X")),
    // … port every remaining row from the reference table, canonical = old canonical_name(variant) …
    ("CNOT", gate(G::CNot, NoArgs, Pairs, "CNOT")),
    ("CX",   gate(G::CX,   NoArgs, Pairs, "CX")),
    ("ZCX",  gate(G::ZCX,  NoArgs, Pairs, "ZCX")),
    // … noise rows (I_ERROR keeps ArgCount::Deferred) …
    // … measure rows via measure()/measure_pairs() …
    // … annotation rows + MPAD (EntryKind::MPad, Optional(1), AtLeastOne, "MPAD") + TICK (None args, Any targets, "TICK") …
];
```

Keep the `use` aliases from the reference (`use GateName as G;` etc.) and add nothing new.

- [ ] **Step 3: Implement `lookup` and `canonical_name`**

```rust
/// Look up a Stim instruction name. `None` means unknown.
pub fn lookup(name: &str) -> Option<TableEntry> {
    TABLE.iter().find(|(n, _)| *n == name).map(|(_, e)| *e)
}

/// Canonical spelling for a decoded instruction kind, derived from the
/// same TABLE rows that drive parsing (single source of truth).
pub fn canonical_name(kind: EntryKind) -> &'static str {
    TABLE
        .iter()
        .find(|(_, e)| e.kind == kind)
        .map(|(_, e)| e.canonical)
        .expect("every EntryKind has a TABLE row (enforced by completeness test)")
}
```

- [ ] **Step 4: Write the drift-proofing tests**

```rust
#[cfg(test)]
mod table_tests {
    use super::*;

    #[test]
    fn every_table_key_is_unique() {
        use std::collections::HashSet;
        let mut seen = HashSet::new();
        for (key, _) in TABLE {
            assert!(seen.insert(*key), "duplicate key {key:?} in TABLE");
        }
    }

    #[test]
    fn canonical_round_trips_through_lookup() {
        // For every row, looking up its canonical spelling yields the same kind.
        for (_, entry) in TABLE {
            let via_canonical = lookup(entry.canonical)
                .unwrap_or_else(|| panic!("canonical {:?} not in TABLE", entry.canonical));
            assert_eq!(via_canonical.kind, entry.kind, "canonical mismatch for {:?}", entry.canonical);
        }
    }

    #[test]
    fn every_variant_has_exactly_one_row() {
        // Build the set of kinds present and assert each expected variant appears.
        // Gate/Noise/Measure/Annotation variant lists mirror the enums; if a
        // variant is added without a row, canonical_name() will panic and this
        // test (which calls it for every kind) will fail.
        for (_, entry) in TABLE {
            let _ = canonical_name(entry.kind); // must not panic
        }
    }
}
```

- [ ] **Step 5: Run to verify pass**

Run: `cargo test -p stim-parser-2 --lib table_tests`
Expected: PASS (3 tests).

- [ ] **Step 6: Commit**

```bash
git add crates/stim-parser-2/src/instructions/mod.rs
git commit -m "feat(stim-parser-2): single-source instruction TABLE with canonical column"
```

---

## Phase 3 — `ast/` (shared payloads + two enums)

### Task 6: Shared payload structs and leaf types

**Files:**
- Create: `crates/stim-parser-2/src/ast/shared.rs`, `crates/stim-parser-2/src/ast/mod.rs`
- Modify: `crates/stim-parser-2/src/lib.rs`
- Reference: `crates/stim-parser/src/ast.rs` (`Target`, `PauliAxis`, `PauliFactor`, `Tag`, `TagParam`), `crates/stim-parser/src/extended/ast.rs` (`Axis`)

- [ ] **Step 1: Write the failing test**

```rust
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
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test -p stim-parser-2 --lib ast::shared`
Expected: FAIL — types not found.

- [ ] **Step 3: Implement `shared.rs`**

Port the leaf types verbatim from the reference:
- `Target` (with `as_qubit` and the `PartialEq<usize>` impl) — from `ast.rs`.
- `PauliAxis` (with `as_char`) and `PauliFactor` — from `ast.rs`.
- `Tag`, `TagParam` — from `ast.rs`.
- `Axis` — from `extended/ast.rs`.

Then add the shared family payload structs (each carries `span: Span` instead of `line: usize`):

```rust
use crate::diagnostics::Span;
use crate::instructions::{AnnotationKind, GateName, MeasureName, NoiseName};

#[derive(Debug, Clone, PartialEq)]
pub struct GateOp {
    pub name: GateName,
    pub tags: Vec<Tag>,
    pub args: Vec<f64>,
    pub targets: Vec<Target>,
    pub span: Span,
}

#[derive(Debug, Clone, PartialEq)]
pub struct NoiseOp {
    pub name: NoiseName,
    pub tags: Vec<Tag>,
    pub args: Vec<f64>,
    pub targets: Vec<usize>,
    pub span: Span,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MeasureOp {
    pub name: MeasureName,
    pub tags: Vec<Tag>,
    pub args: Vec<f64>,
    pub targets: Vec<usize>,
    pub span: Span,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AnnotationOp {
    pub kind: AnnotationKind,
    pub args: Vec<f64>,
    pub targets: Vec<usize>,
    pub span: Span,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MppOp {
    pub tags: Vec<Tag>,
    pub args: Vec<f64>,
    pub products: Vec<Vec<PauliFactor>>,
    pub span: Span,
}
```

- [ ] **Step 4: Wire `ast/mod.rs` and `lib.rs`**

`src/ast/mod.rs`:

```rust
// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

pub mod shared;
pub use shared::*;
```

In `src/lib.rs` add `pub mod ast;`.

- [ ] **Step 5: Run to verify pass**

Run: `cargo test -p stim-parser-2 --lib ast::shared`
Expected: PASS (2 tests).

- [ ] **Step 6: Commit**

```bash
git add crates/stim-parser-2/src
git commit -m "feat(stim-parser-2): shared AST payloads and leaf types"
```

### Task 7: `Instruction` and `Program` (vanilla)

**Files:**
- Create: `crates/stim-parser-2/src/ast/vanilla.rs`
- Modify: `crates/stim-parser-2/src/ast/mod.rs`

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ast::shared::GateOp;
    use crate::diagnostics::{LineMap, Span};
    use crate::instructions::GateName;
    use std::sync::Arc;

    #[test]
    fn program_holds_instructions() {
        let p = Program {
            instructions: vec![Instruction::Gate(GateOp {
                name: GateName::H,
                tags: vec![],
                args: vec![],
                targets: vec![],
                span: Span::new(0, 1),
            })],
            line_map: Arc::new(LineMap::new("H 0")),
        };
        assert_eq!(p.instructions.len(), 1);
    }

    #[test]
    fn program_eq_ignores_line_map() {
        let g = || Instruction::Gate(GateOp {
            name: GateName::H, tags: vec![], args: vec![], targets: vec![], span: Span::new(0, 1),
        });
        let a = Program { instructions: vec![g()], line_map: Arc::new(LineMap::new("H 0")) };
        let b = Program { instructions: vec![g()], line_map: Arc::new(LineMap::new("\n\nH 0")) };
        assert_eq!(a, b); // equality is by instructions only
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test -p stim-parser-2 --lib ast::vanilla`
Expected: FAIL — `Instruction`/`Program` not found.

- [ ] **Step 3: Implement `vanilla.rs`**

```rust
// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Vanilla Stim AST. Tags are preserved verbatim; the parser does not
//! resolve the Stim dialect — that is the consumer's responsibility.

use std::sync::Arc;

use crate::ast::shared::{AnnotationOp, GateOp, MeasureOp, MppOp, NoiseOp, Tag};
use crate::diagnostics::{LineMap, Span};

#[derive(Debug, Clone, PartialEq)]
pub enum Instruction {
    Gate(GateOp),
    Noise(NoiseOp),
    Measure(MeasureOp),
    Annotation(AnnotationOp),
    Mpp(MppOp),
    MPad {
        tags: Vec<Tag>,
        prob: Option<f64>,
        bits: Vec<usize>,
        span: Span,
    },
    Repeat {
        count: u64,
        body: Vec<Instruction>,
        span: Span,
    },
}

/// Vanilla program. Owns the `LineMap` so consumers can resolve any node's
/// `span` to `line:col` (spec §6.3). Equality is by `instructions` only —
/// the line map is positional metadata, not identity.
#[derive(Debug, Clone)]
pub struct Program {
    pub instructions: Vec<Instruction>,
    pub line_map: Arc<LineMap>,
}

impl PartialEq for Program {
    fn eq(&self, other: &Self) -> bool {
        self.instructions == other.instructions
    }
}
```

- [ ] **Step 4: Wire `ast/mod.rs`**

Add `pub mod vanilla;` and `pub use vanilla::{Instruction, Program};`.

- [ ] **Step 5: Run to verify pass**

Run: `cargo test -p stim-parser-2 --lib ast::vanilla`
Expected: PASS (1 test).

- [ ] **Step 6: Commit**

```bash
git add crates/stim-parser-2/src/ast
git commit -m "feat(stim-parser-2): vanilla Instruction/Program AST"
```

### Task 8: `ExtendedInstruction`, `ExtendedProgram`, `measurement_count`

**Files:**
- Create: `crates/stim-parser-2/src/ast/extended.rs`
- Modify: `crates/stim-parser-2/src/ast/mod.rs`
- Reference: `crates/stim-parser/src/extended/ast.rs` (`measurement_count`/`count_in_slice`)

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ast::shared::{MeasureOp, GateOp};
    use crate::diagnostics::{LineMap, Span};
    use crate::instructions::{GateName, MeasureName};
    use std::sync::Arc;

    fn span() -> Span { Span::new(0, 1) }

    #[test]
    fn measurement_count_scales_with_repeat() {
        let m = ExtendedInstruction::Measure(MeasureOp {
            name: MeasureName::M, tags: vec![], args: vec![], targets: vec![0, 1], span: span(),
        });
        let prog = ExtendedProgram {
            instructions: vec![ExtendedInstruction::Repeat {
                count: 3,
                body: vec![m],
                span: span(),
            }],
            line_map: Arc::new(LineMap::new("")),
        };
        assert_eq!(prog.measurement_count(), 6);
    }

    #[test]
    fn gate_op_is_shared_with_vanilla() {
        // The same GateOp struct constructs an ExtendedInstruction::Gate.
        let _ = ExtendedInstruction::Gate(GateOp {
            name: GateName::H, tags: vec![], args: vec![], targets: vec![], span: span(),
        });
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test -p stim-parser-2 --lib ast::extended`
Expected: FAIL — types not found.

- [ ] **Step 3: Implement `extended.rs`**

```rust
// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Typed AST for Stim with PPVM tag-based extensions promoted to
//! first-class instruction variants.

use std::sync::Arc;

use crate::ast::shared::{AnnotationOp, Axis, GateOp, MeasureOp, MppOp, NoiseOp, Tag};
use crate::diagnostics::{LineMap, Span};

#[derive(Debug, Clone, PartialEq)]
pub enum ExtendedInstruction {
    // Pass-through families — the SAME structs as the vanilla AST.
    Gate(GateOp),
    Noise(NoiseOp),
    Measure(MeasureOp),
    Annotation(AnnotationOp),
    Mpp(MppOp),

    // Promoted sugar.
    T { targets: Vec<usize>, span: Span },
    TDag { targets: Vec<usize>, span: Span },
    Rotation { axis: Axis, theta: f64, targets: Vec<usize>, span: Span },
    U3 { theta: f64, phi: f64, lambda: f64, targets: Vec<usize>, span: Span },
    Loss { p: f64, targets: Vec<usize>, span: Span },
    CorrelatedLoss { ps: [f64; 3], targets: Vec<(usize, usize)>, span: Span },
    MPad { tags: Vec<Tag>, prob: Option<f64>, bits: Vec<bool>, span: Span },
    Repeat { count: u64, body: Vec<ExtendedInstruction>, span: Span },
}

/// Extended program. Owns the `LineMap` (spec §6.3) so `ppvm-stim` can
/// resolve a node's `span` to a line for `ExecError`. Equality by
/// `instructions` only.
#[derive(Debug, Clone)]
pub struct ExtendedProgram {
    pub instructions: Vec<ExtendedInstruction>,
    pub line_map: Arc<LineMap>,
}

impl PartialEq for ExtendedProgram {
    fn eq(&self, other: &Self) -> bool {
        self.instructions == other.instructions
    }
}

impl ExtendedProgram {
    /// Total recorded bits the program produces, accounting for REPEAT
    /// factors. Pure AST property; backend-agnostic.
    pub fn measurement_count(&self) -> usize {
        count_in_slice(&self.instructions, 1)
    }
}

fn count_in_slice(instructions: &[ExtendedInstruction], factor: u64) -> usize {
    let mut total = 0usize;
    let factor_usize = usize::try_from(factor).unwrap_or(usize::MAX);
    for instr in instructions {
        match instr {
            ExtendedInstruction::Measure(op) => {
                total = total.saturating_add(op.targets.len().saturating_mul(factor_usize));
            }
            ExtendedInstruction::MPad { bits, .. } => {
                total = total.saturating_add(bits.len().saturating_mul(factor_usize));
            }
            ExtendedInstruction::Mpp(op) => {
                total = total.saturating_add(op.products.len().saturating_mul(factor_usize));
            }
            ExtendedInstruction::Repeat { count, body, .. } => {
                total = total.saturating_add(count_in_slice(body, factor.saturating_mul(*count)));
            }
            ExtendedInstruction::Gate(_)
            | ExtendedInstruction::Noise(_)
            | ExtendedInstruction::Annotation(_)
            | ExtendedInstruction::T { .. }
            | ExtendedInstruction::TDag { .. }
            | ExtendedInstruction::Rotation { .. }
            | ExtendedInstruction::U3 { .. }
            | ExtendedInstruction::Loss { .. }
            | ExtendedInstruction::CorrelatedLoss { .. } => {}
        }
    }
    total
}
```

Note: unlike the reference, `Measure` is now `ExtendedInstruction::Measure(MeasureOp)` (not `Raw(RawPassthrough::Measure{..})`), so the `count_in_slice` match arm reads `op.targets`.

- [ ] **Step 4: Wire `ast/mod.rs`**

Add `pub mod extended;` and `pub use extended::{ExtendedInstruction, ExtendedProgram};`.

- [ ] **Step 5: Run to verify pass**

Run: `cargo test -p stim-parser-2 --lib ast::extended`
Expected: PASS (2 tests).

- [ ] **Step 6: Commit**

```bash
git add crates/stim-parser-2/src/ast
git commit -m "feat(stim-parser-2): extended AST with shared payloads"
```

---

## Phase 4 — `syntax/` (chumsky grammar)

### Task 9: `RawSyntaxNode`, `RawTarget`, parser-stack thread

**Files:**
- Create: `crates/stim-parser-2/src/syntax/raw.rs`, `crates/stim-parser-2/src/syntax/mod.rs`
- Modify: `crates/stim-parser-2/src/lib.rs`
- Reference: `crates/stim-parser/src/parser.rs` (`RawSyntaxNode`, `RawTarget`, `run_on_parser_stack`, `PARSER_STACK_SIZE`)

- [ ] **Step 1: Implement `raw.rs`**

Port `RawSyntaxNode` and `RawTarget` from `crates/stim-parser/src/parser.rs`, changing the `Tag` import to `crate::ast::shared::Tag` and keeping `chumsky::span::SimpleSpan<usize>` for spans (the grammar works in chumsky spans; conversion to `Span` happens in validate). Add a type alias:

```rust
// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use chumsky::span::SimpleSpan;
use crate::ast::shared::Tag;

pub(crate) type RawSyntaxTree = Vec<RawSyntaxNode>;

#[derive(Debug, Clone)]
pub(crate) enum RawSyntaxNode {
    Instruction {
        name: String,
        tags: Vec<Tag>,
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
pub(crate) struct RawTarget {
    pub text: String,
    pub span: SimpleSpan<usize>,
}
```

- [ ] **Step 2: Implement `syntax/mod.rs` with `run_on_parser_stack`**

Port `run_on_parser_stack` and `PARSER_STACK_SIZE` verbatim from `crates/stim-parser/src/parser.rs` (the `#[cfg(target_arch = "wasm32")]` carve-out included):

```rust
// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod grammar;
mod raw;

pub(crate) use raw::{RawSyntaxNode, RawSyntaxTree, RawTarget};

// … paste PARSER_STACK_SIZE const and run_on_parser_stack<R, F> fn here …
pub(crate) use self::stack::run_on_parser_stack;

mod stack {
    // (paste the const + fn; keep them in a submodule or inline — either is fine)
}
```

Simpler: inline the const + `pub(crate) fn run_on_parser_stack` directly in `mod.rs` (no `stack` submodule). Use whichever compiles cleanly.

In `src/lib.rs` add `pub(crate) mod syntax;`.

- [ ] **Step 3: Build**

Run: `cargo build -p stim-parser-2`
Expected: FAIL — `grammar` module is empty/missing (next task). If you stubbed `grammar.rs` empty, expected `Finished`. Create an empty `grammar.rs` with just the license header for now.

- [ ] **Step 4: Commit**

```bash
git add crates/stim-parser-2/src
git commit -m "feat(stim-parser-2): raw syntax nodes and parser-stack thread"
```

### Task 10: Port the chumsky grammar

**Files:**
- Modify: `crates/stim-parser-2/src/syntax/grammar.rs`
- Test: `crates/stim-parser-2/src/syntax/grammar.rs` (`#[cfg(test)]` — port the reference grammar tests)
- Reference: `crates/stim-parser/src/grammar.rs` (entire file)

- [ ] **Step 1: Port the grammar verbatim with two import changes**

Copy the entire body of `crates/stim-parser/src/grammar.rs` into the new `grammar.rs`. Apply exactly these changes:
- `use crate::ast::{Tag, TagParam};` → `use crate::ast::shared::{Tag, TagParam};`
- `use crate::parser::{RawSyntaxNode, RawTarget};` → `use crate::syntax::raw::{RawSyntaxNode, RawTarget};`
- Keep `program_parser`, `instruction_line`, `repeat_block`, all combinators, and the entire `#[cfg(test)] mod tests` unchanged (the test `use crate::parser::RawSyntaxNode;` becomes `use crate::syntax::raw::RawSyntaxNode;`).

Mark `program_parser` `pub(crate)`.

- [ ] **Step 2: Run the ported grammar tests to verify they pass**

Run: `cargo test -p stim-parser-2 --lib grammar`
Expected: PASS (same ~30 grammar tests as the reference: `ident_matches…`, `signed_float_parses…`, `program_parses_repeat_block`, `program_rejects_oversized_repeat_count_without_panicking`, etc.).

- [ ] **Step 3: Commit**

```bash
git add crates/stim-parser-2/src/syntax/grammar.rs
git commit -m "feat(stim-parser-2): port chumsky grammar"
```

---

## Phase 5 — `pipeline/` (typestate + validate + lower)

### Task 11: Pipeline skeleton and the `parse` transition

**Files:**
- Create: `crates/stim-parser-2/src/pipeline/mod.rs`
- Modify: `crates/stim-parser-2/src/lib.rs`
- Reference: `crates/stim-parser/src/parser.rs` (`parse_impl`'s chumsky-driving + error-forwarding)

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagnostics::FailFast;

    #[test]
    fn parse_transition_produces_parsed_state() {
        let mut sink = FailFast::new();
        let parsed = Pipeline::new("H 0\n").parse(&mut sink);
        assert!(parsed.is_ok());
    }

    #[test]
    fn parse_transition_emits_diagnostic_and_aborts_on_syntax_error() {
        let mut sink = FailFast::new();
        let res = Pipeline::new("REPEAT 2 {\nH 0\n").parse(&mut sink); // unclosed
        assert!(res.is_err());
        assert!(sink.saw_error());
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test -p stim-parser-2 --lib pipeline`
Expected: FAIL — `Pipeline` not found.

- [ ] **Step 3: Implement the pipeline skeleton + `parse`**

```rust
// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Typestate lowering pipeline: Source → Parsed → Validated → Lowered.

mod lower;
mod validate;

use std::sync::Arc;

use chumsky::Parser;

use crate::ast::{ExtendedProgram, Program};
use crate::diagnostics::{Diagnostic, DiagnosticSink, Flow, LineMap, Severity, Span};
use crate::diagnostics::Aborted;
use crate::syntax::{self, RawSyntaxTree};

pub struct Pipeline<S> {
    state: S,
}

pub struct Source<'a> {
    src: &'a str,
}
pub struct Parsed {
    pub(crate) tree: RawSyntaxTree,
    pub(crate) line_map: Arc<LineMap>,
}
// Once a Program exists it owns the LineMap (spec §6.3), so the later
// states need only hold the program.
pub struct Validated {
    pub(crate) program: Program,
}
pub struct Lowered {
    pub(crate) program: ExtendedProgram,
}

impl<'a> Pipeline<Source<'a>> {
    pub fn new(src: &'a str) -> Self {
        Pipeline { state: Source { src } }
    }

    /// Stage 1: pure syntax. Forwards every chumsky error into the sink.
    pub fn parse(self, sink: &mut dyn DiagnosticSink) -> Result<Pipeline<Parsed>, Aborted> {
        let src = self.state.src;
        let line_map = Arc::new(LineMap::new(src));
        let result = crate::syntax::program_parser().parse(src);
        match result.into_result() {
            Ok(tree) => Ok(Pipeline { state: Parsed { tree, line_map } }),
            Err(errors) => {
                for err in errors {
                    let span: Span = (*err.span()).into();
                    let flow = sink.emit(Diagnostic::error(span, "syntax", err.to_string()));
                    if flow == Flow::Abort {
                        return Err(Aborted);
                    }
                }
                // All syntax errors forwarded; with a non-aborting sink we still
                // cannot produce a tree, so abort the stage.
                Err(Aborted)
            }
        }
    }
}

impl Pipeline<Validated> {
    pub fn finish(self) -> Program {
        self.state.program
    }
}

impl Pipeline<Lowered> {
    pub fn finish(self) -> ExtendedProgram {
        self.state.program
    }
}
```

Add to `src/lib.rs`: `pub mod pipeline;` and make `syntax::program_parser` reachable — in `src/syntax/mod.rs` add `pub(crate) use grammar::program_parser;`.

Note: `Severity` import is used by later transitions; if the compiler warns it's unused now, remove it and re-add in Task 12/13.

- [ ] **Step 4: Run to verify pass**

Run: `cargo test -p stim-parser-2 --lib pipeline`
Expected: PASS (2 tests). The `validate`/`lower` modules are empty stubs for now — create `src/pipeline/validate.rs` and `src/pipeline/lower.rs` each with just the license header so the `mod` lines compile.

- [ ] **Step 5: Commit**

```bash
git add crates/stim-parser-2/src
git commit -m "feat(stim-parser-2): typestate pipeline skeleton + parse transition"
```

### Task 12: `validate` transition (RawSyntaxTree → Program)

**Files:**
- Modify: `crates/stim-parser-2/src/pipeline/validate.rs`, `crates/stim-parser-2/src/pipeline/mod.rs`
- Test: `crates/stim-parser-2/src/pipeline/validate.rs` (port the reference `validate_tests`)
- Reference: `crates/stim-parser/src/parser.rs` (`validate_program`, `validate_node`, `check_arg_count`, `parse_mpp_products`, `parse_rec`, `qubit_indices`, `build_instruction`)

- [ ] **Step 1: Implement `validate.rs`**

Port the validation logic from `crates/stim-parser/src/parser.rs`, with these transformations:
- The entry point becomes `pub(crate) fn validate(tree: RawSyntaxTree, line_map: &Arc<LineMap>, sink: &mut dyn DiagnosticSink) -> Result<Program, Aborted>`. It builds `Program { instructions, line_map: Arc::clone(line_map) }` so the program owns its line map (spec §6.3).
- Instead of `return Err(ParseError::…)`, emit a `Diagnostic::error(span, code, message)` to `sink` and honor the returned `Flow` (on `Abort` return `Err(Aborted)`; on `Continue` skip the offending instruction — push nothing and proceed). Use these `code`s and spans:
  - unknown instruction → `code = "unknown-instruction"`, span = the instruction-name span. Message: `format!("unknown instruction '{name}'")`.
  - arg-count → `code = "arg-count"`, span = name span. Message mirrors the reference `ParseError::ArgCount` text.
  - target-count → `code = "target-count"`, span = name span. Message mirrors `ParseError::TargetCount`.
  - invalid target → `code = "invalid-target"`, span = the target's span. Message: `format!("invalid target {:?}", t.text)`.
  - invalid MPP target → `code = "invalid-mpp-target"`, span = target span. Message: `format!("invalid MPP target {:?}", t.text)`.
- Convert chumsky `SimpleSpan` → `Span` via `.into()` when constructing diagnostics and when storing spans on AST nodes.
- `build_instruction` now constructs the shared payload structs (`GateOp{..}`, `NoiseOp{..}`, etc.) with `span` set to the instruction-name span (`name_span.into()`), and returns `Instruction`. MPP builds `Instruction::Mpp(MppOp{..})`; MPAD builds `Instruction::MPad{..}`.
- `lookup`/`canonical_name`/`EntryKind`/`ArgCount`/`TargetArity` come from `crate::instructions`.
- Keep `parse_rec` and `parse_mpp_products` and `qubit_indices` as private helpers (the recursion into `Repeat` bodies calls `validate` recursively).

The transformation from "first-error `Result`" to "emit-and-maybe-continue" is mechanical: each `return Err(e)` becomes:

```rust
if sink.emit(diagnostic) == Flow::Abort {
    return Err(Aborted);
}
// Continue: for a per-instruction error, skip pushing this instruction; for a
// per-target error inside the target loop, the instruction is abandoned — emit,
// then on Continue `continue` the outer instruction loop without pushing.
```

- [ ] **Step 2: Add the `validate` transition to `mod.rs`**

```rust
impl Pipeline<Parsed> {
    pub fn validate(self, sink: &mut dyn DiagnosticSink) -> Result<Pipeline<Validated>, Aborted> {
        let Parsed { tree, line_map } = self.state;
        let program = validate::validate(tree, &line_map, sink)?;
        Ok(Pipeline { state: Validated { program } })
    }
}
```

(No `finish` on `Pipeline<Parsed>`: `RawSyntaxTree` is crate-private, so the raw tree is never exposed.)

- [ ] **Step 3: Port the validate tests**

Port `validate_tests` from `crates/stim-parser/src/parser.rs`, adapting:
- Build a `Collect` or `FailFast` sink and call `validate(nodes, &lm(), &mut sink)`.
- Assert on the returned `Program`'s instructions (now `Instruction::Gate(GateOp{ name: GateName::H, .. })` etc.) and, for error cases, on `sink.into_items()[0].code` (e.g. `"unknown-instruction"`, `"arg-count"`, `"target-count"`, `"invalid-target"`).
- The reference test `invalid_target_uses_target_span_for_line_col` becomes: emit to a `Collect` sink, then resolve `items[0].span.line_col(&line_map)` and assert `(2, 5)`.

Example adapted test:

```rust
#[test]
fn unknown_instruction_emits_diagnostic() {
    let mut sink = Collect::new();
    let nodes = vec![instr("FROBNICATE", vec![], vec!["0"], (0, 10))];
    let prog = validate(nodes, &lm(), &mut sink).unwrap(); // Collect never aborts
    assert!(prog.instructions.is_empty());
    let items = sink.into_items();
    assert_eq!(items[0].code, Some("unknown-instruction"));
}
```

Port the `instr`, `instr_with_target_spans`, `instr_with_tags`, `lm` helpers from the reference (adjusting `RawTarget`/`RawSyntaxNode` import to `crate::syntax::raw`).

- [ ] **Step 4: Run to verify pass**

Run: `cargo test -p stim-parser-2 --lib validate`
Expected: PASS (the ported ~13 validate tests).

- [ ] **Step 5: Commit**

```bash
git add crates/stim-parser-2/src/pipeline
git commit -m "feat(stim-parser-2): validate transition (RawSyntaxTree -> Program)"
```

### Task 13: `lower` transition (Program → ExtendedProgram)

**Files:**
- Modify: `crates/stim-parser-2/src/pipeline/lower.rs`, `crates/stim-parser-2/src/pipeline/mod.rs`
- Test: `crates/stim-parser-2/src/pipeline/lower.rs`
- Reference: `crates/stim-parser/src/extended/interpret.rs` (entire file)

- [ ] **Step 1: Implement `lower.rs`**

Port `interpret.rs` with these transformations:
- Entry point: `pub(crate) fn lower(program: Program, sink: &mut dyn DiagnosticSink) -> Result<ExtendedProgram, Aborted>`. Capture `let line_map = Arc::clone(&program.line_map);` up front and build `ExtendedProgram { instructions, line_map }` at the end so the extended program also owns the line map (spec §6.3).
- Input is now the new `Instruction` enum: match `Instruction::Gate(op)` → `interpret_gate(op, sink)`, `Instruction::Noise(op)` → `interpret_noise(op, sink)`, `Instruction::Measure(op)` → `Ok(ExtendedInstruction::Measure(op))` (move the shared struct straight through — no `RawPassthrough`), `Instruction::Annotation(op)` → `Ok(ExtendedInstruction::Annotation(op))`, `Instruction::Mpp(op)` → `Ok(ExtendedInstruction::Mpp(op))`, `Instruction::MPad{..}` → convert bits, `Instruction::Repeat{..}` → recurse.
- `interpret_gate` takes a `GateOp` and returns `ExtendedInstruction`: the no-tag/other-gate fall-throughs return `ExtendedInstruction::Gate(op)` (reusing the same struct); the `T`/`TDag`/`S[T]`/`I[..]` paths build the sugar variants. Replace `qubit_targets(targets, name, line)` with a helper that, on a `Target::Rec`, emits `Diagnostic::error(op.span, "record-target-not-allowed", msg)` to the sink and returns `Err(Aborted)` (or skips on `Continue` — but record-target-on-sugar is fatal for that instruction, so emit then `Err(Aborted)` is the faithful behavior; on `Continue` from a permissive sink, skip the instruction).
- All the `ExtendedParseError::{InvalidTag, InvalidMPadBit, RecordTargetNotAllowed}` returns become `sink.emit(Diagnostic::error(span, code, message))` + `Flow` handling, with codes:
  - invalid tag → `"invalid-tag"` (message: keep the reference's `format!`-built text).
  - invalid MPAD bit → `"invalid-mpad-bit"`.
  - record target not allowed → `"record-target-not-allowed"`.
  - Span for these is the instruction's `op.span` (or the MPad/Mpp `span`).
- `exact_named_params`, `require_no_params`, `convert_mpad_bits`, `pair_targets`, `interpret_identity_tag` port across with the `line: usize` parameter replaced by `span: Span` (used only to build diagnostics).

- [ ] **Step 2: Add the `lower` transition to `mod.rs`**

```rust
impl Pipeline<Validated> {
    pub fn lower(self, sink: &mut dyn DiagnosticSink) -> Result<Pipeline<Lowered>, Aborted> {
        let program = lower::lower(self.state.program, sink)?;
        Ok(Pipeline { state: Lowered { program } })
    }
}
```

- [ ] **Step 3: Write/port lower tests**

Port the gate/noise interpretation assertions exercised indirectly by `tests/extended.rs` as focused unit tests here. Minimum set:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagnostics::{Collect, FailFast};
    use crate::pipeline::Pipeline;

    fn lower_extended(src: &str) -> Result<ExtendedProgram, Vec<crate::diagnostics::Diagnostic>> {
        let mut sink = FailFast::new();
        Pipeline::new(src)
            .parse(&mut sink)
            .and_then(|p| p.validate(&mut sink))
            .and_then(|p| p.lower(&mut sink))
            .map(|p| p.finish())
            .map_err(|_| sink.into_items())
    }

    #[test]
    fn s_t_tag_lowers_to_t() {
        let prog = lower_extended("S[T] 0\n").unwrap();
        assert!(matches!(prog.instructions[0], ExtendedInstruction::T { .. }));
    }

    #[test]
    fn native_t_lowers_to_t() {
        let prog = lower_extended("T 0\n").unwrap();
        assert!(matches!(prog.instructions[0], ExtendedInstruction::T { .. }));
    }

    #[test]
    fn i_rx_lowers_to_rotation() {
        let prog = lower_extended("I[R_X(theta=0.5)] 0\n").unwrap();
        assert!(matches!(prog.instructions[0], ExtendedInstruction::Rotation { .. }));
    }

    #[test]
    fn i_error_loss_lowers_to_loss() {
        let prog = lower_extended("I_ERROR[loss](0.01) 0\n").unwrap();
        assert!(matches!(prog.instructions[0], ExtendedInstruction::Loss { .. }));
    }

    #[test]
    fn t_with_record_target_is_rejected() {
        let err = lower_extended("T rec[-1]\n").unwrap_err();
        assert_eq!(err.last().unwrap().code, Some("record-target-not-allowed"));
    }

    #[test]
    fn i_error_without_tag_is_rejected() {
        let err = lower_extended("I_ERROR(0.01) 0\n").unwrap_err();
        assert_eq!(err.last().unwrap().code, Some("invalid-tag"));
    }
}
```

- [ ] **Step 4: Run to verify pass**

Run: `cargo test -p stim-parser-2 --lib lower`
Expected: PASS (6 tests).

- [ ] **Step 5: Commit**

```bash
git add crates/stim-parser-2/src/pipeline
git commit -m "feat(stim-parser-2): lower transition (Program -> ExtendedProgram)"
```

### Task 14: Tier-1 `parse` / `parse_extended` + `prelude`

**Files:**
- Modify: `crates/stim-parser-2/src/lib.rs`
- Test: `crates/stim-parser-2/src/lib.rs` (`#[cfg(test)]`)

- [ ] **Step 1: Write the failing test**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_returns_program() {
        let prog = parse("H 0\nM 0\n").unwrap();
        assert_eq!(prog.instructions.len(), 2);
    }

    #[test]
    fn parse_extended_returns_extended_program() {
        let prog = parse_extended("S[T] 0\n").unwrap();
        assert_eq!(prog.instructions.len(), 1);
    }

    #[test]
    fn parse_error_renders_line_col() {
        let err = parse("REPEAT 2 {\nH 0\n").unwrap_err();
        assert!(err.to_string().starts_with("error at line"));
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test -p stim-parser-2 --lib tests::parse`
Expected: FAIL — `parse` not found.

- [ ] **Step 3: Implement the Tier-1 functions**

```rust
// in src/lib.rs, after the `pub mod` declarations
use std::sync::Arc;

use crate::ast::{ExtendedProgram, Program};
use crate::diagnostics::{Diagnostics, FailFast, LineMap};
use crate::pipeline::Pipeline;
use crate::syntax::run_on_parser_stack;

/// Parse Stim source into the vanilla [`Program`] AST. Uses a fail-fast
/// policy; the returned [`Diagnostics`] holds the first error.
pub fn parse(src: &str) -> Result<Program, Diagnostics> {
    run_on_parser_stack(|| {
        let mut sink = FailFast::new();
        let result = Pipeline::new(src)
            .parse(&mut sink)
            .and_then(|p| p.validate(&mut sink));
        match result {
            Ok(p) => Ok(p.finish()),
            Err(_) => Err(Diagnostics::new(sink.into_items(), Arc::new(LineMap::new(src)))),
        }
    })
}

/// Parse Stim source into the extended-dialect [`ExtendedProgram`] AST.
pub fn parse_extended(src: &str) -> Result<ExtendedProgram, Diagnostics> {
    run_on_parser_stack(|| {
        let mut sink = FailFast::new();
        let result = Pipeline::new(src)
            .parse(&mut sink)
            .and_then(|p| p.validate(&mut sink))
            .and_then(|p| p.lower(&mut sink));
        match result {
            Ok(p) => Ok(p.finish()),
            Err(_) => Err(Diagnostics::new(sink.into_items(), Arc::new(LineMap::new(src)))),
        }
    })
}

pub mod prelude {
    pub use crate::ast::{
        AnnotationOp, Axis, ExtendedInstruction, ExtendedProgram, GateOp, Instruction, MeasureOp,
        MppOp, NoiseOp, PauliAxis, PauliFactor, Program, Tag, TagParam, Target,
    };
    pub use crate::diagnostics::{Diagnostic, Diagnostics, DiagnosticSink, Flow, Severity};
    pub use crate::instructions::{AnnotationKind, GateName, MeasureName, NoiseName};
    pub use crate::pipeline::Pipeline;
    pub use crate::{parse, parse_extended};
}
```

Ensure `run_on_parser_stack` is re-exported from `syntax`: in `src/syntax/mod.rs` it is `pub(crate) fn`, so `use crate::syntax::run_on_parser_stack;` resolves.

- [ ] **Step 4: Run to verify pass**

Run: `cargo test -p stim-parser-2 --lib`
Expected: PASS (all lib tests, including the 3 new ones).

- [ ] **Step 5: Commit**

```bash
git add crates/stim-parser-2/src/lib.rs
git commit -m "feat(stim-parser-2): Tier-1 parse/parse_extended + prelude"
```

---

## Phase 6 — `print/` (first-class canonical printer)

### Task 15: `StimPrint`, `PrintOptions`, `Display`, `to_stim`

**Files:**
- Create: `crates/stim-parser-2/src/print/mod.rs`
- Modify: `crates/stim-parser-2/src/lib.rs`
- Test: `crates/stim-parser-2/src/print/mod.rs`
- Reference: `crates/stim-parser/src/display.rs` (entire file — but note it has TWO printers `fmt_raw`/`fmt_raw_passthrough`; we collapse them)

- [ ] **Step 1: Write the failing test (the canonical-shape spot checks from the reference)**

```rust
#[cfg(test)]
mod tests {
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
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test -p stim-parser-2 --lib print`
Expected: FAIL — `to_stim` not found.

- [ ] **Step 3: Implement `print/mod.rs`**

Port the formatting from `crates/stim-parser/src/display.rs` with these structural changes:
- Define `PrintOptions { pub indent: std::borrow::Cow<'static, str> }` with `Default` = `Cow::Borrowed("    ")`.
- Define `pub trait StimPrint { fn print(&self, out: &mut dyn fmt::Write, opts: &PrintOptions, depth: usize) -> fmt::Result; }`.
- Implement `StimPrint` for `GateOp`, `NoiseOp`, `MeasureOp`, `AnnotationOp`, `MppOp` (the shared family writers — ported once from `fmt_raw`; this removes the `fmt_raw`/`fmt_raw_passthrough` duplication). The helpers `write_tags`, `write_args`, `write_usize_targets`, `write_targets`, `write_mpp_products`, and the `FloatLit` struct port verbatim (use `name.canonical_name()` → replace with `crate::instructions::canonical_name(EntryKind::Gate(self.name))`, or add an inherent `canonical_name()` accessor on each name enum that delegates to `crate::instructions::canonical_name` — pick one and use consistently).
- Implement `StimPrint` for `Instruction` (matches the shared-op arms by delegating to the op's `StimPrint`, plus `MPad`/`Repeat`) and for `ExtendedInstruction` (shared-op arms delegate; sugar arms ported from `fmt_ext`; `Repeat` recurses; `MPad` prints `bool` bits as `u8::from(bit)`).
- Implement `StimPrint` for `Program` and `ExtendedProgram` (iterate instructions at depth 0).
- Add inherent methods + `Display`:

```rust
impl Program {
    pub fn to_stim(&self) -> String { self.to_stim_with(&PrintOptions::default()) }
    pub fn to_stim_with(&self, opts: &PrintOptions) -> String {
        let mut s = String::new();
        let _ = self.print(&mut s, opts, 0);
        s
    }
}
impl fmt::Display for Program {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.print(f, &PrintOptions::default(), 0)
    }
}
// …identical pair for ExtendedProgram…
```

Indentation uses `opts.indent` repeated `depth` times (replaces the hard-coded `const INDENT`).

For the name-enum canonical accessor: add to `src/instructions/mod.rs`:

```rust
impl GateName { pub fn canonical_name(self) -> &'static str { canonical_name(EntryKind::Gate(self)) } }
impl NoiseName { pub fn canonical_name(self) -> &'static str { canonical_name(EntryKind::Noise(self)) } }
impl MeasureName { pub fn canonical_name(self) -> &'static str { canonical_name(EntryKind::Measure(self)) } }
impl AnnotationKind { pub fn canonical_name(self) -> &'static str { canonical_name(EntryKind::Annotation(self)) } }
```

These let the printer and any consumer get a spelling without touching the table directly, and they keep `GateName: Display` (add a `Display` impl delegating to `canonical_name`, matching the reference).

- [ ] **Step 4: Wire `lib.rs`**

Add `pub mod print;` and to the `prelude` add `pub use crate::print::{PrintOptions, StimPrint};`.

- [ ] **Step 5: Run to verify pass**

Run: `cargo test -p stim-parser-2 --lib print`
Expected: PASS (2 tests).

- [ ] **Step 6: Commit**

```bash
git add crates/stim-parser-2/src
git commit -m "feat(stim-parser-2): first-class canonical StimPrint/Display printer"
```

---

## Phase 7 — Port the integration & property tests

### Task 16: Port the integration tests

**Files:**
- Create: `crates/stim-parser-2/tests/{errors,extended,gates,measure,noise,roundtrip,syntax,tags,record_targets,mpp}.rs`
- Reference: the same-named files under `crates/stim-parser/tests/`

- [ ] **Step 1: Port each integration test file**

For each reference file, copy it into `crates/stim-parser-2/tests/` and apply the API mapping:
- `use stim_parser::prelude::parse;` → `use stim_parser_2::prelude::parse;`
- `use stim_parser::extended::parse_extended;` → `use stim_parser_2::prelude::parse_extended;`
- `use stim_parser::ast::{…};` and `use stim_parser::extended::{…};` → `use stim_parser_2::prelude::{…};`
- AST matching changes: `RawInstruction::Gate { name, .. }` → `Instruction::Gate(GateOp { name, .. })`; `RawPassthrough::Measure { .. }` → `ExtendedInstruction::Measure(MeasureOp { .. })`; etc.
- Error assertions: where a test matched `ParseError::UnknownInstruction { .. }` / `ExtendedParseError::InvalidTag { .. }`, change to inspect the returned `Diagnostics` — `err.iter().next().unwrap().code == Some("unknown-instruction")` (or `"invalid-tag"`, `"arg-count"`, `"target-count"`, `"invalid-target"`, `"invalid-mpad-bit"`, `"record-target-not-allowed"`). The `errors.rs` file is the one most affected; map each asserted error variant to its `code` per the table in Tasks 12–13.
- `format!("{ast}")` round-trip assertions are unchanged (Display is preserved).

Port `roundtrip.rs` last — it should pass unmodified except the imports, since the canonical output format is preserved byte-for-byte.

- [ ] **Step 2: Run the full ported suite**

Run: `cargo test -p stim-parser-2 --tests`
Expected: PASS. Fix any `code`/shape mismatches until green. If a reference test asserted an exact `ParseError` `Display` string, assert the new `Diagnostics` `Display` string instead (format: `error at line L, col C: <message>`).

- [ ] **Step 3: Commit**

```bash
git add crates/stim-parser-2/tests
git commit -m "test(stim-parser-2): port integration tests"
```

### Task 17: Port the proptests

**Files:**
- Create: `crates/stim-parser-2/tests/{proptest_ast,proptest_parse,proptest_roundtrip}.rs`
- Reference: same-named files under `crates/stim-parser/tests/`

- [ ] **Step 1: Port the three proptest files**

Copy each, applying the same import mapping as Task 16. `proptest_roundtrip.rs` uses only `parse`/`parse_extended` + `format!` and the `program_source()` strategy — port unchanged but for imports. `proptest_ast.rs` matches on AST node shapes — apply the `Instruction::Gate(GateOp{..})` mapping. `proptest_parse.rs` likewise.

- [ ] **Step 2: Copy the proptest regressions seed file**

```bash
cp crates/stim-parser/tests/proptest_ast.proptest-regressions crates/stim-parser-2/tests/proptest_ast.proptest-regressions
```

(If the seeds reference shapes that no longer reproduce, proptest will simply not find a failure; that is acceptable. Do not delete the file — keeping it preserves known historical counterexamples.)

- [ ] **Step 3: Run the proptests**

Run: `cargo test -p stim-parser-2 --test proptest_roundtrip --test proptest_ast --test proptest_parse`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add crates/stim-parser-2/tests
git commit -m "test(stim-parser-2): port proptest suites"
```

---

## Phase 8 — Differential parity harness

### Task 18: Parity harness vs old `stim-parser`

**Files:**
- Create: `crates/stim-parser-2/tests/parity.rs`
- Reference: `crates/stim-parser/tests/proptest_roundtrip.rs` (`program_source` strategy), `crates/stim-parser/tests/roundtrip.rs` (corpora)

- [ ] **Step 1: Write the parity harness**

```rust
// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential parity: old `stim_parser` vs new `stim_parser_2`.
//! Asserts (1) same accept/reject and (2) byte-identical canonical print
//! output, over a hand corpus and proptest-generated programs.

use proptest::prelude::*;

// Same fragment strategy as tests/proptest_roundtrip.rs — keep in sync.
fn instruction_fragment() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("H 0\n".to_string()),
        Just("CX 0 1\n".to_string()),
        Just("S[T] 0\n".to_string()),
        Just("I[R_X(theta=0.5)] 0\n".to_string()),
        Just("I_ERROR[loss](0.01) 0\n".to_string()),
        Just("DEPOLARIZE1(0.05) 0\n".to_string()),
        Just("M(0.001) 0\n".to_string()),
        Just("MPAD 0 1 0\n".to_string()),
        Just("DETECTOR rec[-1]\n".to_string()),
        Just("REPEAT 3 {\n    H 0\n    M 0\n}\n".to_string()),
        Just("# leading\n".to_string()),
        Just("H 0  # trail\n".to_string()),
    ]
}

fn program_source() -> impl Strategy<Value = String> {
    prop::collection::vec(instruction_fragment(), 0..16).prop_map(|f| f.concat())
}

fn assert_parity(src: &str) {
    let old = stim_parser::extended::parse_extended(src);
    let new = stim_parser_2::prelude::parse_extended(src);
    match (old, new) {
        (Ok(o), Ok(n)) => {
            assert_eq!(format!("{o}"), n.to_stim(), "print mismatch for:\n{src}");
        }
        (Err(_), Err(_)) => {}
        (o, n) => panic!(
            "accept/reject mismatch for:\n{src}\n  old_ok={}, new_ok={}",
            o.is_ok(),
            n.is_ok()
        ),
    }
}

#[test]
fn parity_on_hand_corpus() {
    for src in [
        "H 0\nCX 0 1\nM 0 1\n",
        "S[T] 0\nI[R_X(theta=0.25)] 1\nI_ERROR[loss](0.01) 2\n",
        "REPEAT 2 {\n    REPEAT 3 {\n        H 0\n        M 0\n    }\n}\n",
        "MPAD 0 1 0\nMPAD(0.01) 1 1 0 0\n",
        "I_ERROR[correlated_loss](0.1, 0.05, 0.05) 0 1 2 3\n",
        "R 0 1\nMR 0\nDETECTOR rec[-1]\nOBSERVABLE_INCLUDE(0) rec[-1]\nTICK\n",
    ] {
        assert_parity(src);
    }
}

proptest! {
    #[test]
    fn parity_on_generated_programs(src in program_source()) {
        assert_parity(&src);
    }
}
```

- [ ] **Step 2: Run the parity harness**

Run: `cargo test -p stim-parser-2 --test parity`
Expected: PASS. If a print mismatch appears, the new printer diverged from the old `Display`; fix `print/mod.rs` until byte-identical. If an accept/reject mismatch appears, the new validate/lower diverged; reconcile against the reference logic.

- [ ] **Step 3: Commit**

```bash
git add crates/stim-parser-2/tests/parity.rs
git commit -m "test(stim-parser-2): differential parity harness vs stim-parser"
```

### Task 19: Full-workspace gate before swap

- [ ] **Step 1: Run the whole workspace**

Run: `cargo test --workspace`
Expected: PASS (old crate, new crate, and every other crate still green; `ppvm-stim` still on the old `stim-parser`).

- [ ] **Step 2: Run clippy on the new crate**

Run: `cargo clippy -p stim-parser-2 --all-targets -- -D warnings`
Expected: no warnings. Fix any.

- [ ] **Step 3: Commit any clippy fixes**

```bash
git add crates/stim-parser-2
git commit -m "chore(stim-parser-2): clippy clean"
```

---

## Phase 9 — Swap

> Only start this phase once Task 18 and Task 19 are green. This is the point of no return for `ppvm-stim`'s dependency.

### Task 20: Repoint `ppvm-stim` at `stim-parser-2`

**Files:**
- Modify: `crates/ppvm-stim/Cargo.toml`, `crates/ppvm-stim/src/lib.rs`, `crates/ppvm-stim/src/executor.rs`, `crates/ppvm-stim/src/validate.rs`, `crates/ppvm-stim/benches/tableau-msd-stim.rs`
- Reference: the import/match sites found via grep (Task 0 of the review found them)

- [ ] **Step 1: Switch the dependency**

In `crates/ppvm-stim/Cargo.toml` change:

```toml
stim-parser = { version = "0.1.0", path = "../stim-parser" }
```
to:
```toml
stim-parser-2 = { version = "0.1.0", path = "../stim-parser-2" }
```

- [ ] **Step 2: Update imports and AST match sites**

In each `ppvm-stim` source file, apply:
- `use stim_parser::…` → `use stim_parser_2::…` (path within is `prelude`, or `ast`/`extended` re-exports; use `stim_parser_2::prelude::{…}`).
- `RawPassthrough::Gate { name, targets, .. }` → `ExtendedInstruction::Gate(GateOp { name, targets, .. })`; same for `Noise`/`Measure`/`Annotation` (these are now top-level `ExtendedInstruction` variants wrapping the shared `*Op` structs, **not** nested under `Raw`).
- Remove `ExtendedInstruction::Raw(...)` wrapping entirely — match the family variants directly.
- `*line` field reads (e.g. in `validate.rs`'s `ExecError` construction): `ExtendedProgram` already owns `line_map` (Task 8), so resolve `op.span.line(&program.line_map)` to recover the `usize` line. `ExecError`'s `line: usize` fields and `Display` strings stay exactly as they are — only the *source* of the number changes from a stored `line` field to `span.line(&line_map)`. Thread `&program.line_map` into `validate_slice` alongside the existing `measurements` accumulator.
- `ExtendedInstruction::Mpp { products, .. }` → `ExtendedInstruction::Mpp(MppOp { products, .. })`.

- [ ] **Step 3: Build and test ppvm-stim**

Run: `cargo test -p ppvm-stim`
Expected: PASS. Iterate on the match sites until green.

- [ ] **Step 4: Full workspace**

Run: `cargo test --workspace`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add crates/ppvm-stim
git commit -m "refactor(ppvm-stim): switch to stim-parser-2"
```

### Task 21: Delete old crate, rename new → `stim-parser`

**Files:**
- Delete: `crates/stim-parser/`
- Rename: `crates/stim-parser-2/` → `crates/stim-parser/`
- Modify: root `Cargo.toml`, `crates/stim-parser/Cargo.toml`, `crates/ppvm-stim/Cargo.toml`, all `ppvm-stim` imports, `crates/stim-parser-2/tests/parity.rs` (delete)

- [ ] **Step 1: Remove the parity harness (it dev-depends on the old crate)**

```bash
git rm crates/stim-parser-2/tests/parity.rs
```

- [ ] **Step 2: Delete the old crate and remove its dev-dep**

```bash
git rm -r crates/stim-parser
```

In `crates/stim-parser-2/Cargo.toml`, delete the line:
```toml
stim-parser = { version = "0.1.0", path = "../stim-parser" }   # parity harness only; removed at swap
```

- [ ] **Step 3: Rename the directory and the crate**

```bash
git mv crates/stim-parser-2 crates/stim-parser
```

In `crates/stim-parser/Cargo.toml` change `name = "stim-parser-2"` → `name = "stim-parser"`.

- [ ] **Step 4: Update workspace members and the consumer dependency**

In root `Cargo.toml`, change the members line back to a single entry — replace `"crates/ppvm-stim", "crates/stim-parser", "crates/stim-parser-2",` with `"crates/ppvm-stim", "crates/stim-parser",`.

In `crates/ppvm-stim/Cargo.toml` change:
```toml
stim-parser-2 = { version = "0.1.0", path = "../stim-parser-2" }
```
to:
```toml
stim-parser = { version = "0.1.0", path = "../stim-parser" }
```

- [ ] **Step 5: Rename the crate references in `ppvm-stim` source**

Replace every `stim_parser_2` with `stim_parser` across `crates/ppvm-stim/src/**` and `crates/ppvm-stim/benches/**`:

Run: `grep -rl stim_parser_2 crates/ppvm-stim` then edit each, or use your editor's project replace. Verify none remain:
Run: `grep -rn "stim_parser_2\|stim-parser-2" crates/ Cargo.toml`
Expected: no matches.

- [ ] **Step 6: Full workspace build + test**

Run: `cargo test --workspace`
Expected: PASS.

Run: `cargo clippy --workspace --all-targets -- -D warnings`
Expected: no warnings.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "refactor(stim-parser): replace with stim-parser-2 implementation"
```

---

## Self-Review (completed by plan author)

**Spec coverage:**
- Typestate pipeline → Tasks 11–14. ✓
- Diagnostic-sink effect model (Flow as continuation) → Tasks 2–3, threaded in 12–13. ✓
- Two enums sharing payloads → Tasks 6–8 (shared `*Op` structs reused). ✓
- Single-source table → Tasks 4–5 + completeness tests. ✓
- `Span` over `line` → Task 1, carried on all ops (6–8), resolved in diagnostics. ✓
- First-class canonical printer → Task 15 (`to_stim`/`to_stim_with`/`Display`/`StimPrint`). ✓
- Keep chumsky + stack thread → Tasks 9–10. ✓
- Validation-responsibility contract unchanged → validate (12) keeps dialect checks; `ppvm-stim` keeps capability checks (20). ✓
- Parity-then-swap delivery → Tasks 18–21. ✓
- Non-goals (no trivia, no hand-written parser, no streaming) → respected; printer is canonical-only. ✓

**Placeholder scan:** No "TBD"/"add error handling" left. Port tasks name the exact reference file and the exact transformation. Line-map ownership is fixed upfront — `Program` (Task 7) and `ExtendedProgram` (Task 8) own `Arc<LineMap>` from the start, so `ppvm-stim` resolves `op.span.line(&program.line_map)` at swap time (Task 20) with no open decision.

**Type consistency:** `Diagnostic::error(span, code, message)` signature is used identically in Tasks 2, 12, 13. `sink.emit(..) -> Flow` consistent. `canonical_name(EntryKind)` free fn (Task 5) + per-enum `canonical_name(self)` accessors (Task 15) are consistent. `ExtendedInstruction::Measure(MeasureOp)` (not `Raw`) consistent across Tasks 8, 13, 15, 16, 20. `Program`/`ExtendedProgram` own `Arc<LineMap>` with hand-written `PartialEq` (instructions only) — defined in Tasks 7–8, produced in Tasks 12–13, consumed in Task 20. Pipeline state structs `Validated`/`Lowered` hold only the program (which owns the line map); `Parsed` holds `tree` + `line_map`.

**Scope:** One cohesive crate reimplementation with a clean swap. Sequential layers, not independent subsystems — correctly a single plan.
