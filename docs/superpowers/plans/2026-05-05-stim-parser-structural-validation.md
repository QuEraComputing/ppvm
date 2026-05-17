# Move Stim Structural Validation Into stim-parser Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Tighten `stim-parser`'s extended AST so `MPAD` bits and `I_ERROR[correlated_loss]` target pairs are structurally valid before `ppvm-stim` sees them.

**Architecture:** `stim-parser` converts `MPAD` targets into `Vec<bool>` and `CorrelatedLoss` targets into `Vec<(usize, usize)>` during extended interpretation. `ppvm-stim::prepare` then validates only ppvm-tableau backend support and returns `Result<(), ExecError>`, while measurement counting stays on `ExtendedProgram::measurement_count()`.

**Tech Stack:** Rust 2024 edition, `thiserror`, existing `stim-parser`, `ppvm-stim`, `ppvm-python-native`, and Python tests through `uv`.

---

## File Structure

**`crates/stim-parser/src/extended/ast.rs`**: change extended AST field types:
- `ExtendedInstruction::MPad { bits: Vec<bool>, .. }`
- `ExtendedInstruction::CorrelatedLoss { targets: Vec<(usize, usize)>, .. }`

**`crates/stim-parser/src/extended/parser.rs`**: add `ExtendedParseError::InvalidMPadBit`.

**`crates/stim-parser/src/extended/interpret.rs`**: validate and convert `MPAD` bits; pair `CorrelatedLoss` targets after the existing nonzero-even check.

**`crates/stim-parser/tests/extended.rs`**: assert bool MPAD bits, paired correlated-loss targets, and parse-time MPAD bit rejection.

**`crates/ppvm-stim/src/prepare.rs`**: remove structural error variants and structural checks; `prepare` returns `Result<(), ExecError>`.

**`crates/ppvm-stim/src/executor.rs`**: call `prepare(program)?` then `program.measurement_count()`; consume bool MPAD bits and correlated-loss target pairs directly.

**`crates/ppvm-stim/tests/prepare.rs`**: keep only backend unsupported-instruction tests; add one supported-structural smoke test.

**`crates/ppvm-python-native/src/stim_program.rs`** and **`crates/ppvm-python-native/Cargo.toml`**: adapt `prepare` call and use `ppvm_stim::ExtendedProgram` re-export so the native Python crate does not need a direct `stim-parser` dependency.

---

## Task 1: Tighten `stim-parser` Extended AST

**Files:**
- Modify: `crates/stim-parser/src/extended/ast.rs`
- Modify: `crates/stim-parser/src/extended/parser.rs`
- Modify: `crates/stim-parser/src/extended/interpret.rs`
- Modify: `crates/stim-parser/tests/extended.rs`

- [ ] **Step 1: Update parser tests first**

In `crates/stim-parser/tests/extended.rs`, replace `vanilla_mpad_passes_through` with:

```rust
#[test]
fn vanilla_mpad_passes_through_as_bool_bits() {
    let p = parse_ok("MPAD 0 1\n");
    match &p.instructions[0] {
        ExtendedInstruction::MPad { bits, prob, .. } => {
            assert_eq!(bits, &vec![false, true]);
            assert!(prob.is_none());
        }
        other => panic!("{other:?}"),
    }
}
```

Add this test near the MPAD test:

```rust
#[test]
fn mpad_non_bit_target_errors_in_extended_parser() {
    let err = parse_err("MPAD 0 2 1\n");
    match err {
        ExtendedParseError::InvalidMPadBit { line, index, value } => {
            assert_eq!(line, 1);
            assert_eq!(index, 1);
            assert_eq!(value, 2);
        }
        other => panic!("{other:?}"),
    }
}
```

Replace `i_error_correlated_loss_one_arg_expands` with:

```rust
#[test]
fn i_error_correlated_loss_one_arg_expands_and_pairs_targets() {
    let p = parse_ok("I_ERROR[correlated_loss](0.5) 0 1 2 3\n");
    match &p.instructions[0] {
        ExtendedInstruction::CorrelatedLoss { ps, targets, line } => {
            approx_eq(ps[0], 0.5);
            approx_eq(ps[1], 0.0);
            approx_eq(ps[2], 0.0);
            assert_eq!(targets, &vec![(0, 1), (2, 3)]);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}
```

- [ ] **Step 2: Run the focused parser tests and confirm failure**

Run:

```bash
cargo test -p stim-parser --test extended mpad_non_bit_target_errors_in_extended_parser
cargo test -p stim-parser --test extended vanilla_mpad_passes_through_as_bool_bits
cargo test -p stim-parser --test extended i_error_correlated_loss_one_arg_expands_and_pairs_targets
```

Expected: FAIL. The first test fails because `ExtendedParseError::InvalidMPadBit` does not exist. The other tests fail because the AST still uses `Vec<usize>`.

- [ ] **Step 3: Add the parser error variant**

In `crates/stim-parser/src/extended/parser.rs`, change `ExtendedParseError` to:

```rust
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
    #[error("'MPAD' at line {line} target #{index} = {value}, must be 0 or 1")]
    InvalidMPadBit {
        line: usize,
        index: usize,
        value: usize,
    },
}
```

- [ ] **Step 4: Tighten the extended AST field types**

In `crates/stim-parser/src/extended/ast.rs`, change the relevant enum fields to:

```rust
MPad {
    tags: Vec<Tag>,
    prob: Option<f64>,
    bits: Vec<bool>,
    line: usize,
},
```

and:

```rust
CorrelatedLoss {
    ps: [f64; 3],
    targets: Vec<(usize, usize)>,
    line: usize,
},
```

Leave `measurement_count()` unchanged; it only needs `bits.len()`.

- [ ] **Step 5: Convert MPAD bits during extended interpretation**

In `crates/stim-parser/src/extended/interpret.rs`, replace the `RawInstruction::MPad` arm with:

```rust
RawInstruction::MPad {
    tags,
    prob,
    bits,
    line,
} => Ok(ExtendedInstruction::MPad {
    tags: tags.clone(),
    prob: *prob,
    bits: convert_mpad_bits(bits, *line)?,
    line: *line,
}),
```

Add this helper near `require_no_params`:

```rust
fn convert_mpad_bits(bits: &[usize], line: usize) -> Result<Vec<bool>, ExtendedParseError> {
    let mut out = Vec::with_capacity(bits.len());
    for (index, value) in bits.iter().copied().enumerate() {
        match value {
            0 => out.push(false),
            1 => out.push(true),
            _ => return Err(ExtendedParseError::InvalidMPadBit { line, index, value }),
        }
    }
    Ok(out)
}
```

- [ ] **Step 6: Pair correlated-loss targets during extended interpretation**

In the `IError[correlated_loss]` arm in `crates/stim-parser/src/extended/interpret.rs`, replace the output construction with:

```rust
Ok(ExtendedInstruction::CorrelatedLoss {
    ps,
    targets: pair_targets(targets),
    line,
})
```

Add this helper near `convert_mpad_bits`:

```rust
fn pair_targets(targets: &[usize]) -> Vec<(usize, usize)> {
    targets
        .chunks_exact(2)
        .map(|pair| (pair[0], pair[1]))
        .collect()
}
```

The existing nonzero-even validation stays before this helper is called.

- [ ] **Step 7: Run parser tests**

Run:

```bash
cargo test -p stim-parser --test extended
cargo test -p stim-parser
```

Expected: PASS.

- [ ] **Step 8: Format and lint parser crate**

Run:

```bash
cargo fmt -p stim-parser
cargo clippy -p stim-parser --all-targets -- -D warnings
```

Expected: PASS with no diagnostics.

- [ ] **Step 9: Commit parser AST tightening**

```bash
git add crates/stim-parser/src/extended/ast.rs crates/stim-parser/src/extended/parser.rs crates/stim-parser/src/extended/interpret.rs crates/stim-parser/tests/extended.rs
git commit -m "refactor(stim-parser): validate extended structural invariants"
```

---

## Task 2: Collapse `ppvm-stim::prepare` to Backend Support Validation

**Files:**
- Modify: `crates/ppvm-stim/src/prepare.rs`
- Modify: `crates/ppvm-stim/src/executor.rs`
- Modify: `crates/ppvm-stim/tests/prepare.rs`
- Modify: `crates/ppvm-stim/tests/executor.rs` if compiler errors reveal type assumptions

- [ ] **Step 1: Replace prepare tests**

Replace `crates/ppvm-stim/tests/prepare.rs` with:

```rust
use ppvm_stim::{ExecError, parse_extended, prepare};

fn err_from_src(src: &str) -> ExecError {
    let prog = parse_extended(src).expect("parse_extended");
    prepare(&prog).expect_err("must reject")
}

#[test]
fn unsupported_swap_rejected() {
    let e = err_from_src("SWAP 0 1");
    match e {
        ExecError::Unsupported { name, line } => {
            assert_eq!(name, "SWAP");
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn unsupported_mx_rejected() {
    let e = err_from_src("MX 0");
    assert!(matches!(e, ExecError::Unsupported { .. }));
}

#[test]
fn unsupported_heralded_erase_rejected() {
    let e = err_from_src("HERALDED_ERASE(0.1) 0");
    assert!(matches!(e, ExecError::Unsupported { .. }));
}

#[test]
fn supported_structural_instructions_are_not_rejected_by_prepare() {
    let prog = parse_extended("MPAD 0 1\nI_ERROR[correlated_loss](0.5) 0 1\n")
        .expect("parse_extended");
    assert_eq!(prepare(&prog), Ok(()));
}
```

- [ ] **Step 2: Run prepare tests and confirm failure**

Run:

```bash
cargo test -p ppvm-stim --test prepare
```

Expected: FAIL because `prepare` still returns `Result<usize, ExecError>` and `ExecError` still has structural variants.

- [ ] **Step 3: Simplify `ExecError` and `prepare`**

In `crates/ppvm-stim/src/prepare.rs`, change `ExecError` to:

```rust
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum ExecError {
    #[error("unsupported instruction '{name}' at line {line}")]
    Unsupported { name: String, line: usize },
}
```

Change `prepare` to:

```rust
pub fn prepare(program: &ExtendedProgram) -> Result<(), ExecError> {
    validate_slice(&program.instructions)
}
```

In `validate_slice`, remove the `MPad` bit loop and `CorrelatedLoss` arity check. The match tail should become:

```rust
ExtendedInstruction::Repeat { body, .. } => validate_slice(body)?,
ExtendedInstruction::Annotation { .. }
| ExtendedInstruction::MPad { .. }
| ExtendedInstruction::T { .. }
| ExtendedInstruction::TDag { .. }
| ExtendedInstruction::Rotation { .. }
| ExtendedInstruction::U3 { .. }
| ExtendedInstruction::Loss { .. }
| ExtendedInstruction::CorrelatedLoss { .. } => {}
_ => unreachable!("ExtendedInstruction variant added but not handled in prepare"),
```

- [ ] **Step 4: Update executor pre-sizing and dispatch**

In `crates/ppvm-stim/src/executor.rs`, change `execute` from:

```rust
let count = prepare(program)?;
let mut results = Vec::with_capacity(count);
```

to:

```rust
prepare(program)?;
let mut results = Vec::with_capacity(program.measurement_count());
```

Change `sample` from:

```rust
let count = prepare(program)?;
Ok((0..num_shots)
    .map(|_| {
        let mut tab = make_tableau();
        let mut results = Vec::with_capacity(count);
        execute_slice(&program.instructions, &mut tab, &mut results);
        results
    })
    .collect())
```

to:

```rust
prepare(program)?;
let count = program.measurement_count();
Ok((0..num_shots)
    .map(|_| {
        let mut tab = make_tableau();
        let mut results = Vec::with_capacity(count);
        execute_slice(&program.instructions, &mut tab, &mut results);
        results
    })
    .collect())
```

Change correlated-loss dispatch from:

```rust
ExtendedInstruction::CorrelatedLoss { ps, targets, .. } => {
    let ps: [T::Coeff; 3] = [ps[0].into(), ps[1].into(), ps[2].into()];
    for (a, b) in targets.iter().copied().tuples() {
        tab.correlated_loss_channel(a, b, ps.clone());
    }
}
```

to:

```rust
ExtendedInstruction::CorrelatedLoss { ps, targets, .. } => {
    let ps: [T::Coeff; 3] = [ps[0].into(), ps[1].into(), ps[2].into()];
    for &(a, b) in targets {
        tab.correlated_loss_channel(a, b, ps.clone());
    }
}
```

Change MPAD dispatch from:

```rust
ExtendedInstruction::MPad { bits, prob, .. } => {
    let noise = prob.unwrap_or(0.0);
    for &bit in bits {
        let bit_bool = bit != 0;
        let recorded = if noise > 0.0 && tab.bernoulli(noise) {
            !bit_bool
        } else {
            bit_bool
        };
        results.push(Some(recorded));
    }
}
```

to:

```rust
ExtendedInstruction::MPad { bits, prob, .. } => {
    let noise = prob.unwrap_or(0.0);
    for &bit in bits {
        let recorded = if noise > 0.0 && tab.bernoulli(noise) {
            !bit
        } else {
            bit
        };
        results.push(Some(recorded));
    }
}
```

- [ ] **Step 5: Run ppvm-stim tests**

Run:

```bash
cargo test -p ppvm-stim --test prepare
cargo test -p ppvm-stim
```

Expected: PASS.

- [ ] **Step 6: Format and lint ppvm-stim**

Run:

```bash
cargo fmt -p ppvm-stim
cargo clippy -p ppvm-stim --all-targets -- -D warnings
```

Expected: PASS with no diagnostics.

- [ ] **Step 7: Commit ppvm-stim boundary cleanup**

```bash
git add crates/ppvm-stim/src/prepare.rs crates/ppvm-stim/src/executor.rs crates/ppvm-stim/tests/prepare.rs crates/ppvm-stim/tests/executor.rs
git commit -m "refactor(ppvm-stim): keep prepare backend-specific"
```

---

## Task 3: Update Python Native Consumer and Remove Direct Parser Dependency

**Files:**
- Modify: `crates/ppvm-python-native/src/stim_program.rs`
- Modify: `crates/ppvm-python-native/Cargo.toml`
- Modify: `Cargo.lock`
- Optional if needed: `ppvm-python/src/ppvm/stim_program.py`

- [ ] **Step 1: Update native StimProgram parse flow**

In `crates/ppvm-python-native/src/stim_program.rs`, replace:

```rust
use stim_parser::extended::ExtendedProgram;

use ppvm_stim::{parse_extended, prepare};
```

with:

```rust
use ppvm_stim::{ExtendedProgram, parse_extended, prepare};
```

Change `parse` from:

```rust
let program = parse_extended(src).map_err(stim_to_pyerr)?;
let measurement_count = prepare(&program).map_err(stim_to_pyerr_exec)?;
Ok(Self {
    program,
    measurement_count,
})
```

to:

```rust
let program = parse_extended(src).map_err(stim_to_pyerr)?;
prepare(&program).map_err(stim_to_pyerr_exec)?;
let measurement_count = program.measurement_count();
Ok(Self {
    program,
    measurement_count,
})
```

- [ ] **Step 2: Remove direct `stim-parser` dependency**

In `crates/ppvm-python-native/Cargo.toml`, delete:

```toml
stim-parser = { version = "0.1.0", path = "../stim-parser" }
```

Run:

```bash
cargo check -p ppvm-python-native
```

Expected: PASS and `Cargo.lock` removes `stim-parser` from the direct dependency list for `ppvm-python-native`.

- [ ] **Step 3: Run Python-facing tests**

Run:

```bash
uv run --project ppvm-python --group dev pytest ppvm-python/test/generalized_tableau/test_stim.py
```

Expected: PASS.

- [ ] **Step 4: Format and lint native crate**

Run:

```bash
cargo fmt -p ppvm-python-native
cargo clippy -p ppvm-python-native --all-targets -- -D warnings
```

Expected: PASS with no diagnostics.

- [ ] **Step 5: Commit Python native cleanup**

```bash
git add crates/ppvm-python-native/src/stim_program.rs crates/ppvm-python-native/Cargo.toml Cargo.lock ppvm-python/src/ppvm/stim_program.py
git commit -m "refactor(python-native): use prepared stim parser re-export"
```

If `ppvm-python/src/ppvm/stim_program.py` is unchanged, omit it from `git add`.

---

## Task 4: Workspace Verification

**Files:**
- No source edits expected unless verification reveals a missed call site.

- [ ] **Step 1: Search for stale structural errors and old AST assumptions**

Run:

```bash
rg -n "InvalidMPadBit|InvalidCorrelatedLossArity|CorrelatedLoss \\{ ps, targets: vec|MPad \\{[^\\n]*bits: vec|bit != 0|tuples\\(\\).*correlated_loss|prepare\\(&.*\\).*measurement" crates ppvm-python
```

Expected: no stale `ExecError::InvalidMPadBit` or `ExecError::InvalidCorrelatedLossArity` references. Any remaining `ExtendedParseError::InvalidMPadBit` references should be in `stim-parser`.

- [ ] **Step 2: Run full Rust formatting**

Run:

```bash
cargo fmt --all
```

Expected: no diff beyond already formatted files.

- [ ] **Step 3: Run full Rust lint**

Run:

```bash
cargo clippy --workspace --all-targets -- -D warnings
```

Expected: PASS.

- [ ] **Step 4: Run full Rust tests**

Run:

```bash
cargo test --workspace
```

Expected: PASS.

- [ ] **Step 5: Run Python tests**

Run:

```bash
uv run --project ppvm-python --group dev pytest ppvm-python/test/
```

Expected: PASS.

- [ ] **Step 6: Commit any verification follow-up**

If verification required fixes:

```bash
git add <changed-files>
git commit -m "test(stim): cover parser-owned structural validation"
```

If no files changed, skip this commit.

---

## Design Notes

- `stim-parser` remains backend-agnostic. It only rejects malformed extended-dialect structure: `MPAD` targets must be bits, and correlated-loss targets must form pairs.
- `ppvm-stim::prepare` remains the ppvm-tableau capability gate. Unsupported valid Stim instructions still fail there, not in `stim-parser`.
- `ExtendedProgram::measurement_count()` remains the only measurement-count API. It is called directly by `execute`, `sample`, and Python `StimProgram.parse`.
- The public `prepare` function keeps its name but changes from `Result<usize, ExecError>` to `Result<(), ExecError>` because counting is no longer part of preparation.
- `ExecError` intentionally collapses to one variant: `Unsupported { name, line }`.

## Self-Review

- Spec coverage: MPAD validation moves to parser through `Vec<bool>` and `InvalidMPadBit`; correlated-loss arity is represented by `Vec<(usize, usize)>`; `prepare` remains backend-specific; measurement count stays in `stim-parser`.
- Placeholder scan: no TBD/TODO placeholders.
- Type consistency: `MPad.bits` is `Vec<bool>` in parser tests and executor dispatch; `CorrelatedLoss.targets` is `Vec<(usize, usize)>` in parser tests and executor dispatch; `prepare` returns `Result<(), ExecError>` everywhere.
