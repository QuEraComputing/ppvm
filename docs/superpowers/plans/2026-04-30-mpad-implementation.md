# MPAD Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote `MPAD` from a no-op annotation to a first-class instruction in the parse → normalize → execute pipeline of the `ppvm-stim` crate, honouring Stim's optional flip-probability argument and rejecting non-bit targets.

**Architecture:** Layered change with compileable commits: parser shape first, normalized representation plus deterministic execution second, statistical noise verification last. Each layer adds a dedicated MPAD variant alongside the existing `Annotation` / `Measure` variants -- no field overloading. Parser carries `bits: Vec<usize>` (so target-validation diagnostics keep the offending value), normalize converts to `Vec<bool>` and validates each value is `0` or `1`, executor pushes one `Some(bool)` per bit and applies optional flip noise via `tab.bernoulli`.

**Tech Stack:** Rust 2024, chumsky 0.12 (already in use), thiserror 2.

**Spec:** `docs/superpowers/specs/2026-04-30-mpad-implementation-design.md`.

**Branch:** `david/ppvm-stim-2` (existing). All work continues on this branch.

**Standing rules for every commit:**
- Run `cargo fmt --all` *before* every `git commit`.
- Run `cargo clippy --workspace --all-targets -- -D warnings` before every commit; fix warnings rather than silencing.
- Test corpus (`tests/data/`, `tests/regen-stim/`, `tests/stim_corpus.rs`) is **not** modified.
- Commit messages follow the imperative-mood convention used on this branch (e.g. `Add MPAD parser variant`).

---

## File-by-file overview

| File | Action | Responsibility after change |
|---|---|---|
| `crates/ppvm-stim/src/parser/ast.rs` | Modify | Drop `AnnotationKind::MPad` (and its `canonical_name` arm). Add `RawInstruction::MPad { tags, prob, bits, line }`. |
| `crates/ppvm-stim/src/parser/table.rs` | Modify | Add `TableEntry::MPad { args, targets }` variant. Reclassify the `MPAD` row from `Annotation` to `MPad` with `ArgCount::Optional(1)` / `TargetArity::AtLeastOne`. |
| `crates/ppvm-stim/src/parser/mod.rs` | Modify | Add 5th arm to `arity_of`. Add 5th arm to `build_instruction` returning `RawInstruction::MPad`. |
| `crates/ppvm-stim/src/tableau_program.rs` | Modify | Add `Instruction::MPad { bits: Vec<bool>, noise: f64, line }`. |
| `crates/ppvm-stim/src/normalize.rs` | Modify | Add `RawInstruction::MPad` arm in `normalize_slice`. Add `NormalizeError::InvalidMPadTarget`. Bump `expected_measurement_count` by `bits.len() * enclosing_repeat_factor`. |
| `crates/ppvm-stim/src/executor.rs` | Modify | Add `Instruction::MPad` arm in `execute_slice`: push one `Some(bit)` per pad, optionally flipped by `tab.bernoulli(noise)`. |
| `crates/ppvm-stim/src/lib.rs` | No change | Existing `pub use` for `Instruction`, `RawInstruction`, `NormalizeError` covers the new variants. |
| `crates/ppvm-stim/tests/parser_measure.rs` | Test additions | Parse-positive cases for MPAD. |
| `crates/ppvm-stim/tests/parser_errors.rs` | Test additions | Parse-error cases for MPAD (zero targets, two args, `rec[-1]` target). |
| `crates/ppvm-stim/tests/normalize.rs` | Test additions | Normalize success + `InvalidMPadTarget` failure. |
| `crates/ppvm-stim/tests/executor.rs` | Test additions | Execution: deterministic pads, in-`REPEAT` pads, noisy pad distribution. Drop stale `MPAD` mention from comment at line 305. |

No new files are created.

---

## Reference: relevant existing patterns

These are the call shapes the new code mirrors. If any of these have drifted in `main` since the spec was written, follow what is in the file at task time and keep the structure.

**1. `arity_of` in `crates/ppvm-stim/src/parser/mod.rs:128-151`** — 4-arm match returning `(ArgCount, TargetArity, &'static str)`. Adding a `TableEntry::MPad` arm requires a hardcoded `"MPAD"` literal because there is no `MPadName` enum (MPAD has only one form).

**2. `build_instruction` in `crates/ppvm-stim/src/parser/mod.rs:298-334`** — 4-arm match constructing `RawInstruction`. Adding the `MPad` arm splits `args` into an `Option<f64>` and uses `targets` directly as `Vec<usize>`:

```rust
TableEntry::MPad { .. } => RawInstruction::MPad {
    tags,
    prob: args.into_iter().next(),
    bits: targets,
    line,
},
```

**3. `Measure` arm in `crates/ppvm-stim/src/normalize.rs:68-88`** — pattern for bumping `*measure_count` by `targets.len() * enclosing_repeat_factor` (saturating). MPAD reuses this formula with `bits.len()` substituted for `targets.len()`.

**4. `MeasureKind::MR` arm in `crates/ppvm-stim/src/executor.rs:225-240`** — pattern for applying readout noise via `tab.bernoulli`. The MPAD body is the same `if *noise > 0.0 && tab.bernoulli(*noise) { !bit } else { bit }` shape, pushed as `Some(_)`.

**5. `measure_noise_distribution_within_3_sigma` in `crates/ppvm-stim/tests/executor.rs:383-403`** — pattern for shot-loop statistical tests using `ppvm_stim::sample` with `GeneralizedTableau::new_with_seed`. The MPAD noise test is structurally identical, only the input circuit changes.

---

## Task 1: Parser layer — accept `MPAD(p) 0 1 0`

Add the parser-level MPAD shape and the table entry that drives validation. After this task, `parse("MPAD 0")` returns `RawInstruction::MPad`, `parse("MPAD")` is a target-count error, and `parse("MPAD(0.1, 0.2) 0")` is an arg-count error.

**Files:**
- Modify: `crates/ppvm-stim/src/parser/ast.rs` (`AnnotationKind`, `RawInstruction`)
- Modify: `crates/ppvm-stim/src/parser/table.rs` (`TableEntry`, the `MPAD` row)
- Modify: `crates/ppvm-stim/src/parser/mod.rs` (`arity_of`, `build_instruction`)
- Modify: `crates/ppvm-stim/tests/parser_measure.rs` (positive parse tests)
- Modify: `crates/ppvm-stim/tests/parser_errors.rs` (negative parse tests)

- [ ] **Step 1: Write the failing parser tests**

Append to `crates/ppvm-stim/tests/parser_measure.rs`:

```rust
#[test]
fn parse_mpad_no_args_no_tags() {
    let p = parse("MPAD 0 1 0").unwrap();
    let RawInstruction::MPad { tags, prob, bits, .. } = &p.instructions[0] else {
        panic!("{:?}", p.instructions[0]);
    };
    assert!(tags.is_empty());
    assert_eq!(*prob, None);
    assert_eq!(bits, &[0usize, 1, 0]);
}

#[test]
fn parse_mpad_with_prob() {
    let p = parse("MPAD(0.25) 0 1").unwrap();
    let RawInstruction::MPad { prob, bits, .. } = &p.instructions[0] else {
        panic!()
    };
    assert_eq!(*prob, Some(0.25));
    assert_eq!(bits, &[0usize, 1]);
}
```

Append to `crates/ppvm-stim/tests/parser_errors.rs`:

```rust
#[test]
fn mpad_zero_targets_is_target_count_error() {
    let err = parse("MPAD").expect_err("must reject");
    match err {
        ParseError::TargetCount { name, divisor, found, .. } => {
            assert_eq!(name, "MPAD");
            assert_eq!(divisor, 1);
            assert_eq!(found, 0);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn mpad_two_args_is_arg_count_error() {
    let err = parse("MPAD(0.1, 0.2) 0").expect_err("must reject");
    match err {
        ParseError::ArgCount { name, expected, found, .. } => {
            assert_eq!(name, "MPAD");
            assert_eq!(expected, 1);
            assert_eq!(found, 2);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn mpad_with_rec_target_is_syntax_error() {
    let err = parse("MPAD rec[-1]").expect_err("must reject");
    assert!(matches!(err, ParseError::Syntax(_)));
}
```

The new tests need `RawInstruction::MPad` to compile. The import line at the top of `parser_measure.rs` is currently:

```rust
use ppvm_stim::{AnnotationKind, MeasureName, ParseError, RawInstruction, parse};
```

The import already includes `RawInstruction`, so no edit is required for that line. Compilation will fail because `RawInstruction::MPad` does not yet exist.

- [ ] **Step 2: Run the new tests to confirm they fail to compile**

Run: `cargo test -p ppvm-stim --test parser_measure parse_mpad -- --nocapture`
Expected: compile error — `RawInstruction::MPad` not found.

- [ ] **Step 3: Drop `AnnotationKind::MPad`**

In `crates/ppvm-stim/src/parser/ast.rs`:

Remove the `MPad` variant from `AnnotationKind` (currently between `Detector` and `ObservableInclude`):

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotationKind {
    Detector,
    ObservableInclude,
    QubitCoords,
    ShiftCoords,
    Tick,
}
```

Delete the comment line above the removed variant ("`MPAD` is treated as an annotation, not a measurement (matches today's `stim.rs`)"); it is stale.

Remove the `MPad => "MPAD"` arm from `AnnotationKind::canonical_name` so the enum has exactly five arms.

- [ ] **Step 4: Add `RawInstruction::MPad`**

In the same file, add a new variant to `RawInstruction` (place it between `Annotation` and `Repeat` so visited order matches lowering order):

```rust
MPad {
    tags: Vec<Tag>,
    prob: Option<f64>,
    bits: Vec<usize>,
    line: usize,
},
```

- [ ] **Step 5: Add `TableEntry::MPad`**

In `crates/ppvm-stim/src/parser/table.rs`, add a new variant to `TableEntry` (place it after `Annotation`):

```rust
MPad { args: ArgCount, targets: TargetArity },
```

Replace the `"MPAD"` row inside `TABLE` (currently `TableEntry::Annotation { kind: AnnotationKind::MPad, … }`) with:

```rust
(
    "MPAD",
    TableEntry::MPad {
        args: ArgCount::Optional(1),
        targets: TargetArity::AtLeastOne,
    },
),
```

- [ ] **Step 6: Wire `arity_of` and `build_instruction`**

In `crates/ppvm-stim/src/parser/mod.rs`, extend `arity_of` (around line 128) with a 5th arm:

```rust
TableEntry::MPad { args, targets } => (args, targets, "MPAD"),
```

Extend `build_instruction` (around line 298) with a 5th arm:

```rust
TableEntry::MPad { .. } => RawInstruction::MPad {
    tags,
    prob: args.into_iter().next(),
    bits: targets,
    line,
},
```

The existing flags upstream of `build_instruction` (`tolerate_non_numeric` at line 199, `skip_arg_validation` at line 217) match only `TableEntry::Annotation { .. }`. Because `MPAD` is no longer that variant, both flags correctly stay `false` for MPAD: targets must parse as `usize`, and `ArgCount::Optional(1)` is enforced.

- [ ] **Step 7: Run the parser tests**

Run: `cargo test -p ppvm-stim --test parser_measure parse_mpad -- --nocapture`
Run: `cargo test -p ppvm-stim --test parser_errors mpad -- --nocapture`
Expected: all five new tests pass.

- [ ] **Step 8: Run the full ppvm-stim test suite**

Run: `cargo test -p ppvm-stim`
Expected: pre-existing tests still pass; the only new passing tests are the five added in Step 1.

If any pre-existing test fails, it most likely matches `AnnotationKind::MPad` somewhere — fix at the call site by removing the `MPad` arm or replacing it with the new `RawInstruction::MPad` shape. Do not reintroduce `AnnotationKind::MPad`.

- [ ] **Step 9: Format, lint, and commit**

Run: `cargo fmt --all`
Run: `cargo clippy --workspace --all-targets -- -D warnings`

Commit:

```bash
git add crates/ppvm-stim/src/parser/ast.rs \
        crates/ppvm-stim/src/parser/table.rs \
        crates/ppvm-stim/src/parser/mod.rs \
        crates/ppvm-stim/tests/parser_measure.rs \
        crates/ppvm-stim/tests/parser_errors.rs
git commit -m "Promote MPAD parser variant out of AnnotationKind"
```

---

## Task 2: Normalize + deterministic executor layer

Add the normalized variant, wire `RawInstruction::MPad` -> `Instruction::MPad` with 0/1 validation, and handle deterministic execution. Rust enum exhaustiveness couples the normalized `Instruction` enum to the executor match, so these changes are one commit. After this task, `normalize::to_tableau` rejects `MPAD 2`, accepts `MPAD 0 1 0`, and `execute(&prog, &mut tab)` for `MPAD 0 1` returns `[Some(false), Some(true)]`.

**Files:**
- Modify: `crates/ppvm-stim/src/tableau_program.rs` (add `Instruction::MPad`)
- Modify: `crates/ppvm-stim/src/normalize.rs` (add error, add normalize arm)
- Modify: `crates/ppvm-stim/src/executor.rs` (add `Instruction::MPad` arm)
- Modify: `crates/ppvm-stim/tests/normalize.rs` (add tests)
- Modify: `crates/ppvm-stim/tests/executor.rs` (deterministic execution tests)

- [ ] **Step 1: Write the failing normalize tests**

Append to `crates/ppvm-stim/tests/normalize.rs`:

```rust
#[test]
fn mpad_normalize_zero_one_succeeds() {
    let p = norm("MPAD 0 1 0");
    let Instruction::MPad { bits, noise, .. } = &p.instructions[0] else {
        panic!("{:?}", p.instructions[0]);
    };
    assert_eq!(bits, &[false, true, false]);
    assert_eq!(*noise, 0.0);
    assert_eq!(p.expected_measurement_count, 3);
}

#[test]
fn mpad_normalize_with_prob() {
    let p = norm("MPAD(0.25) 1 0");
    let Instruction::MPad { bits, noise, .. } = &p.instructions[0] else {
        panic!()
    };
    assert_eq!(bits, &[true, false]);
    assert_eq!(*noise, 0.25);
    assert_eq!(p.expected_measurement_count, 2);
}

#[test]
fn mpad_inside_repeat_block_multiplies_count() {
    let p = norm("REPEAT 3 {\n    MPAD 1\n}");
    assert_eq!(p.expected_measurement_count, 3);
}

#[test]
fn mpad_target_two_is_invalid_target() {
    let err = norm_err("MPAD 0 2 1");
    match err {
        NormalizeError::InvalidMPadTarget { line, index, value } => {
            assert_eq!(line, 1);
            assert_eq!(index, 1);
            assert_eq!(value, 2);
        }
        other => panic!("{other:?}"),
    }
}
```

These tests need `Instruction::MPad` and `NormalizeError::InvalidMPadTarget` to compile.

- [ ] **Step 2: Write the failing deterministic executor tests**

Append to `crates/ppvm-stim/tests/executor.rs`:

```rust
#[test]
fn mpad_single_zero_appends_some_false() {
    let (results, _) = run("MPAD 0", 1);
    assert_eq!(results, vec![Some(false)]);
}

#[test]
fn mpad_single_one_appends_some_true() {
    let (results, _) = run("MPAD 1", 1);
    assert_eq!(results, vec![Some(true)]);
}

#[test]
fn mpad_multi_bit_in_order() {
    let (results, _) = run("MPAD 0 1 0 1", 1);
    assert_eq!(
        results,
        vec![Some(false), Some(true), Some(false), Some(true)]
    );
}

#[test]
fn mpad_interleaved_with_measurement() {
    let (results, _) = run("X 0\nMPAD 1\nM 0\nMPAD 0", 1);
    assert_eq!(results, vec![Some(true), Some(true), Some(false)]);
}

#[test]
fn mpad_inside_repeat_block_executes_each_iteration() {
    let (results, _) = run("REPEAT 3 {\n    MPAD 1\n}", 1);
    assert_eq!(results, vec![Some(true), Some(true), Some(true)]);
}
```

These tests need `RawInstruction::MPad` from Task 1, then `Instruction::MPad` and its executor arm from this task.

- [ ] **Step 3: Run the new tests to confirm they fail to compile**

Run: `cargo test -p ppvm-stim --test normalize mpad -- --nocapture`
Expected: compile error -- `Instruction::MPad` and `NormalizeError::InvalidMPadTarget` not found.

Run: `cargo test -p ppvm-stim --test executor mpad -- --nocapture`
Expected: compile error -- `Instruction::MPad` not found.

- [ ] **Step 4: Add `Instruction::MPad`**

In `crates/ppvm-stim/src/tableau_program.rs`, add a new variant inside the `Instruction` enum (place it after `Measure`, before `Annotation`, so the file's variant order matches the lowering order):

```rust
MPad {
    /// Bits to append to the measurement record, in order.
    bits: Vec<bool>,
    /// Probability of the *recorded* bit being flipped (mirrors `Measure::noise`).
    /// 0.0 means deterministic.
    noise: f64,
    line: usize,
},
```

Update the doc comment on `expected_measurement_count` (currently `crates/ppvm-stim/src/tableau_program.rs:7-9`) to mention MPAD:

```rust
/// Sum over `M` / `MZ` / `MR` / `MPAD` target counts, multiplied by enclosing
/// `REPEAT` counts. Used to pre-size shot result buffers.
pub expected_measurement_count: usize,
```

- [ ] **Step 5: Add `NormalizeError::InvalidMPadTarget`**

In `crates/ppvm-stim/src/normalize.rs`, append a new variant to the `NormalizeError` enum:

```rust
#[error("'MPAD' at line {line} target #{index} = {value}, must be 0 or 1")]
InvalidMPadTarget {
    line: usize,
    index: usize,
    value: usize,
},
```

- [ ] **Step 6: Add the normalize arm**

In `normalize_slice`, add a new arm to the `match raw` block (place it after `RawInstruction::Annotation`, before `RawInstruction::Repeat`):

```rust
RawInstruction::MPad {
    tags: _,
    prob,
    bits,
    line,
} => {
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
        noise: (*prob).unwrap_or(0.0),
        line: *line,
    });
}
```

Borrow note: in this match the destructured fields are references — `prob: &Option<f64>`, `bits: &Vec<usize>`, `line: &usize`. `Option<f64>: Copy` (because `f64: Copy`), so `(*prob)` is a cheap copy, mirroring how the existing `Measure` arm uses `*name` and `*line`. Calling `.unwrap_or(0.0)` directly on `&Option<f64>` does not compile — `Option::unwrap_or` takes `self` by value and method resolution does not auto-deref into a by-value receiver.

The `tags` field is currently unused in normalize: MPAD has no defined tag semantics today. Leaving it parsed-and-discarded mirrors how `RawInstruction::Measure` ignores tags too. If a future tag (e.g. dialect extension) is needed, it would be added here.

- [ ] **Step 7: Add the executor arm**

In `crates/ppvm-stim/src/executor.rs`, inside `execute_slice` (around line 119, the `match instr` block), add a new arm — place it after `Instruction::Measure { … }` and before `Instruction::Annotation { .. }`:

```rust
Instruction::MPad { bits, noise, .. } => {
    for &bit in bits {
        let recorded = if *noise > 0.0 && tab.bernoulli(*noise) {
            !bit
        } else {
            bit
        };
        results.push(Some(recorded));
    }
}
```

Borrow shapes: in this context `bits: &Vec<bool>`, so `for &bit in bits` deref-binds `bit: bool`; `noise: &f64`, so `*noise > 0.0` and `tab.bernoulli(*noise)` are correct.

- [ ] **Step 8: Run the normalize tests**

Run: `cargo test -p ppvm-stim --test normalize mpad -- --nocapture`
Expected: all four new tests pass.

- [ ] **Step 9: Run the deterministic executor tests**

Run: `cargo test -p ppvm-stim --test executor mpad -- --nocapture`
Expected: all five new tests pass.

- [ ] **Step 10: Run the full ppvm-stim test suite**

Run: `cargo test -p ppvm-stim`
Expected: every test passes.

- [ ] **Step 11: Format, lint, and commit**

Run: `cargo fmt --all`
Run: `cargo clippy --workspace --all-targets -- -D warnings`

```bash
git add crates/ppvm-stim/src/tableau_program.rs \
        crates/ppvm-stim/src/normalize.rs \
        crates/ppvm-stim/src/executor.rs \
        crates/ppvm-stim/tests/normalize.rs \
        crates/ppvm-stim/tests/executor.rs
git commit -m "Add MPAD normalize and executor arms"
```

---

## Task 3: Noise distribution test + comment cleanup

Add a 3σ statistical check that MPAD's noise argument actually flips with the requested probability, and remove the stale `MPAD` mention from the existing `test_stim_noop_instructions` comment.

**Files:**
- Modify: `crates/ppvm-stim/tests/executor.rs` (new test, comment fix)

- [ ] **Step 1: Add the distribution test**

Append to `crates/ppvm-stim/tests/executor.rs`:

```rust
#[test]
fn mpad_noise_distribution_within_3_sigma() {
    use ppvm_stim::sample;
    // MPAD(0.3) 0 — pad value is 0; recorded bit flips to 1 with prob 0.3.
    let prog = parse("MPAD(0.3) 0").unwrap();
    let tprog = normalize::to_tableau(&prog).unwrap();
    let n = 4096usize;
    let mut seed_counter: u64 = 0;
    let shots = sample::<_, _, _, _>(&tprog, n, || {
        seed_counter += 1;
        GeneralizedTableau::<ByteFxHashF64<1>, usize>::new_with_seed(1, 1e-10, seed_counter)
    })
    .unwrap();
    let ones = shots.iter().filter(|s| s[0] == Some(true)).count();
    let mean = (n as f64) * 0.3;
    let std = ((n as f64) * 0.3 * 0.7).sqrt();
    assert!(
        ((ones as f64) - mean).abs() < 3.0 * std,
        "got {ones} ones, expected mean {mean} +/- 3*{std}"
    );
}
```

This test should pass directly because Task 2 already wired up `tab.bernoulli`. Run it as a verification, not a fail-then-fix step.

- [ ] **Step 2: Run the distribution test**

Run: `cargo test -p ppvm-stim --test executor mpad_noise_distribution -- --nocapture`
Expected: PASS.

If it fails statistically (rare — ≈0.27% probability under correct behaviour), re-run once. If it fails twice in a row, the noise plumbing is wrong; debug `tab.bernoulli` invocation in the MPAD arm.

- [ ] **Step 3: Drop the stale MPAD mention from the no-op comment**

In `crates/ppvm-stim/tests/executor.rs`, the existing comment at line 305 reads:

```rust
    // TICK, DETECTOR, QUBIT_COORDS, SHIFT_COORDS, MPAD, OBSERVABLE_INCLUDE should not crash
```

The actual source string in that test does not include MPAD. Drop `MPAD,` from the comment:

```rust
    // TICK, DETECTOR, QUBIT_COORDS, SHIFT_COORDS, OBSERVABLE_INCLUDE should not crash
```

Leave the test body unchanged.

- [ ] **Step 4: Run the full ppvm-stim test suite one last time**

Run: `cargo test -p ppvm-stim`
Expected: every test passes, including the new MPAD distribution test.

- [ ] **Step 5: Format, lint, and commit**

Run: `cargo fmt --all`
Run: `cargo clippy --workspace --all-targets -- -D warnings`

```bash
git add crates/ppvm-stim/tests/executor.rs
git commit -m "Verify MPAD noise distribution; clean stale annotation comment"
```

---

## Self-review checklist

After completing all three tasks, run a final review to confirm spec coverage:

- [ ] **Spec section "Decisions" item 1** (dedicated AST variant) — covered by Task 1 (parser), Task 2 (normalized).
- [ ] **Spec item 2** (probability arg implemented) — covered by Task 1 (parser accepts `Optional(1)`), Task 2 (normalize threads `prob` into `noise` and executor applies bernoulli flip).
- [ ] **Spec item 3** (non-0/1 caught in normalize) — covered by Task 2's `InvalidMPadTarget` and the `mpad_target_two_is_invalid_target` test.
- [ ] **Spec item 4** (parser `Vec<usize>` → normalized `Vec<bool>`) — covered by Task 1's `RawInstruction::MPad` and Task 2's normalize arm.
- [ ] **Spec item 5** (`expected_measurement_count` bumped with repeat factor) — covered by Task 2's `mpad_inside_repeat_block_multiplies_count` test.
- [ ] **Spec "Compatibility" — `AnnotationKind::MPad` removed** — covered by Task 1 Step 3.
- [ ] **Spec "Compatibility" — stale comment cleaned up** — covered by Task 3 Step 3.
- [ ] **Spec "Testing" table — every row has a corresponding step** — verified in tasks 1, 2, 3 collectively.

If any item is unchecked, return to the relevant task and add the missing piece.
