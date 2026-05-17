# Collapse `ppvm-stim` IR Layers — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate the intermediate `TableauProgram` IR in `ppvm-stim`. Make the executor consume `stim_parser::extended::ExtendedProgram` directly. Replace `normalize::to_tableau` with two up-front pre-walks: `ExtendedProgram::measurement_count` (in stim-parser) and `ppvm_stim::prepare` (backend-capability validation).

**Architecture:** stim-parser gains a pure `measurement_count()` method on `ExtendedProgram` (counts `Measure`/`MPad` targets multiplied by `REPEAT` factors). ppvm-stim gains a `prepare(&ExtendedProgram) -> Result<usize, ExecError>` that walks the program once to enforce ppvm-tableau's supported-gate set and the structural validations (MPad bits ∈ {0,1}, CorrelatedLoss arity), then returns the measurement count. The executor's `execute`/`sample` take `&ExtendedProgram` directly, call `prepare` once at entry to validate and pre-size the result vec, then dispatch infallibly via match arms unioned by alias (e.g. `H | HXZ → tab.h(q)`).

**Tech Stack:** Rust 2024 edition, `thiserror` for errors, `criterion` for benches, `pyo3` for Python bindings. No new dependencies.

**Spec:** `docs/superpowers/specs/2026-05-05-collapse-ppvm-stim-irs-design.md`

---

## File Structure (post-implementation)

**`crates/stim-parser/src/extended/ast.rs`** (modified): existing struct/enum definitions unchanged; appends a new `impl ExtendedProgram { pub fn measurement_count(&self) -> usize }` block plus a private `count_in_slice` helper.

**`crates/ppvm-stim/src/prepare.rs`** (NEW, ~110 lines): owns `pub enum ExecError`, `pub fn prepare(&ExtendedProgram) -> Result<usize, ExecError>`, plus private `validate_slice` recursion and three small `check_*_supported` helpers.

**`crates/ppvm-stim/src/executor.rs`** (rewritten): `execute`/`sample` take `&ExtendedProgram`; `execute_slice` becomes infallible (`-> ()`), matches on `ExtendedInstruction` directly, fans out to tableau methods via alias-unioned `GateName`/`NoiseName`/`MeasureName` arms.

**`crates/ppvm-stim/src/lib.rs`** (modified): drops `pub mod normalize` and `pub mod tableau_program`; adds `pub mod prepare`. Re-exports `prepare` and `ExecError` from `prepare`. `pub enum Error` drops `Normalize(NormalizeError)`. `run_string` body becomes `parse_extended → execute`.

**Deleted files:**
- `crates/ppvm-stim/src/normalize.rs`
- `crates/ppvm-stim/src/tableau_program.rs`
- `crates/ppvm-stim/tests/normalize.rs`

**Test files:**
- `crates/ppvm-stim/tests/prepare.rs` (NEW, ~80 lines, 6 tests)
- `crates/stim-parser/tests/extended.rs` (modified — append 3 tests)
- `crates/ppvm-stim/tests/executor.rs` (modified — drop normalize step from helper + tests)
- `crates/ppvm-stim/tests/run.rs` (modified — rename one test, swap one error variant in assertion)
- `crates/ppvm-stim/tests/stim_corpus.rs` (modified — `Expect::NormalizeUnsupported` → `Expect::ExecUnsupported`, drop `normalize` import)
- `crates/ppvm-stim/benches/tableau-msd-stim.rs` (modified — bench fn signature + setup)

**`crates/ppvm-python-native/src/stim_program.rs`** (modified): `PyStimProgram { program: ExtendedProgram, measurement_count: usize }`; `parse()` calls `prepare`, errors at parse time.

**`crates/ppvm-python-native/src/interface_tableau.rs`** (modified): two call sites `&prog.inner` → `&prog.program`.

---

## Task 1: Add `ExtendedProgram::measurement_count` to stim-parser

**Files:**
- Modify: `crates/stim-parser/src/extended/ast.rs` (append impl + helper)
- Modify: `crates/stim-parser/tests/extended.rs` (append 3 tests)

- [ ] **Step 1: Append three failing tests to `crates/stim-parser/tests/extended.rs`**

Open the file and append the following at the end:

```rust
// ----------------------------------------------------------------
// measurement_count
// ----------------------------------------------------------------

#[test]
fn measurement_count_counts_m_mz_mr() {
    let p = parse_ok("X 0\nM 0 1 2\nMR 5");
    assert_eq!(p.measurement_count(), 4);
}

#[test]
fn measurement_count_includes_repeat_multiplier() {
    let p = parse_ok("REPEAT 10 {\n    X 0\n    M 0 1\n}");
    assert_eq!(p.measurement_count(), 20);
}

#[test]
fn measurement_count_mpad_inside_repeat_block_multiplies() {
    let p = parse_ok("REPEAT 3 {\n    MPAD 1\n}");
    assert_eq!(p.measurement_count(), 3);
}
```

The `parse_ok` helper already exists at the top of `tests/extended.rs:6-8`; reuse it.

- [ ] **Step 2: Run the new tests to confirm they fail**

Run: `cargo test -p stim-parser --test extended measurement_count`
Expected: FAIL — compile error `no method named measurement_count found for struct ExtendedProgram`.

- [ ] **Step 3: Implement `measurement_count` in `crates/stim-parser/src/extended/ast.rs`**

Append the following after the existing `pub enum Axis { … }` definition (after line 92):

```rust
impl ExtendedProgram {
    /// Total number of recorded bits the program will produce, accounting for
    /// `REPEAT` factors. Pure AST property; backend-agnostic.
    pub fn measurement_count(&self) -> usize {
        count_in_slice(&self.instructions, 1)
    }
}

fn count_in_slice(instructions: &[ExtendedInstruction], factor: u64) -> usize {
    let mut total = 0usize;
    for instr in instructions {
        match instr {
            ExtendedInstruction::Measure { targets, .. } => {
                total = total.saturating_add(
                    targets.len().saturating_mul(factor as usize),
                );
            }
            ExtendedInstruction::MPad { bits, .. } => {
                total = total.saturating_add(
                    bits.len().saturating_mul(factor as usize),
                );
            }
            ExtendedInstruction::Repeat { count, body, .. } => {
                total = total.saturating_add(count_in_slice(
                    body,
                    factor.saturating_mul(*count),
                ));
            }
            ExtendedInstruction::Gate { .. }
            | ExtendedInstruction::Noise { .. }
            | ExtendedInstruction::Annotation { .. }
            | ExtendedInstruction::T { .. }
            | ExtendedInstruction::TDag { .. }
            | ExtendedInstruction::Rotation { .. }
            | ExtendedInstruction::U3 { .. }
            | ExtendedInstruction::Loss { .. }
            | ExtendedInstruction::CorrelatedLoss { .. } => {}
            _ => unreachable!(
                "ExtendedInstruction variant added but not handled in measurement_count"
            ),
        }
    }
    total
}
```

- [ ] **Step 4: Run the new tests to confirm they pass**

Run: `cargo test -p stim-parser --test extended measurement_count`
Expected: PASS — all three tests green.

- [ ] **Step 5: Run the full stim-parser test suite to confirm no regression**

Run: `cargo test -p stim-parser`
Expected: PASS.

- [ ] **Step 6: Format and lint**

Run:
```
cargo fmt -p stim-parser
cargo clippy -p stim-parser --all-targets -- -D warnings
```
Expected: both succeed with no diagnostics.

- [ ] **Step 7: Commit**

```bash
git add crates/stim-parser/src/extended/ast.rs crates/stim-parser/tests/extended.rs
git commit -m "$(cat <<'EOF'
feat(stim-parser): add ExtendedProgram::measurement_count

Public AST utility for downstream consumers — counts recorded bits
across Measure/MPad targets, accounting for REPEAT factors.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Add `prepare` and `ExecError` to ppvm-stim (additive)

This task adds `prepare.rs` alongside the existing `normalize.rs` and `tableau_program.rs`. The old IR keeps working until Task 3 deletes it.

**Files:**
- Modify: `crates/ppvm-stim/src/executor.rs` (delete old empty `enum ExecError {}` and import from `prepare`)
- Modify: `crates/ppvm-stim/src/lib.rs` (add `pub mod prepare`, re-export from prepare)
- Create: `crates/ppvm-stim/src/prepare.rs`
- Create: `crates/ppvm-stim/tests/prepare.rs`

- [ ] **Step 1: Write the 6 failing prepare tests**

Create `crates/ppvm-stim/tests/prepare.rs` with the following content:

```rust
use ppvm_stim::{ExecError, parse_extended, prepare};
use stim_parser::extended::{ExtendedInstruction, ExtendedProgram};

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
fn correlated_loss_single_target_rejected() {
    let prog = ExtendedProgram {
        instructions: vec![ExtendedInstruction::CorrelatedLoss {
            ps: [0.5, 0.0, 0.0],
            targets: vec![0],
            line: 7,
        }],
    };
    let e = prepare(&prog).expect_err("must reject");
    assert_eq!(
        e,
        ExecError::InvalidCorrelatedLossArity { line: 7, found: 1 }
    );
}

#[test]
fn correlated_loss_odd_targets_rejected() {
    let prog = ExtendedProgram {
        instructions: vec![ExtendedInstruction::CorrelatedLoss {
            ps: [0.5, 0.0, 0.0],
            targets: vec![0, 1, 2],
            line: 7,
        }],
    };
    let e = prepare(&prog).expect_err("must reject");
    assert_eq!(
        e,
        ExecError::InvalidCorrelatedLossArity { line: 7, found: 3 }
    );
}

#[test]
fn mpad_target_two_rejected() {
    let e = err_from_src("MPAD 0 2 1");
    match e {
        ExecError::InvalidMPadBit { line, index, value } => {
            assert_eq!(line, 1);
            assert_eq!(index, 1);
            assert_eq!(value, 2);
        }
        other => panic!("{other:?}"),
    }
}
```

- [ ] **Step 2: Run the new tests to confirm they fail**

Run: `cargo test -p ppvm-stim --test prepare`
Expected: FAIL — compile error (`prepare`, `ExecError` variants not yet in scope).

- [ ] **Step 3: Create `crates/ppvm-stim/src/prepare.rs`**

Create the file with the following content:

```rust
//! Validate an [`ExtendedProgram`] against the ppvm-tableau backend's
//! capabilities. On success, returns the measurement count so [`execute`]
//! and [`sample`] can pre-size their result buffers.
//!
//! [`ExtendedProgram`]: stim_parser::extended::ExtendedProgram
//! [`execute`]: crate::executor::execute
//! [`sample`]: crate::executor::sample

use stim_parser::ast::{GateName, MeasureName, NoiseName};
use stim_parser::extended::{ExtendedInstruction, ExtendedProgram};

#[derive(Debug, thiserror::Error, PartialEq)]
pub enum ExecError {
    #[error("unsupported instruction '{name}' at line {line}")]
    Unsupported { name: String, line: usize },

    #[error("'MPAD' at line {line} target #{index} = {value}, must be 0 or 1")]
    InvalidMPadBit {
        line: usize,
        index: usize,
        value: usize,
    },

    #[error(
        "'I_ERROR[correlated_loss]' at line {line} expected nonzero target count divisible by 2, got {found}"
    )]
    InvalidCorrelatedLossArity { line: usize, found: usize },
}

/// Validate the program (rejecting unsupported instructions, malformed MPad
/// bits, and bad CorrelatedLoss arity) and return its measurement count.
///
/// Call once per program, before any shot loop. After this returns `Ok`,
/// the executor can dispatch infallibly.
pub fn prepare(program: &ExtendedProgram) -> Result<usize, ExecError> {
    validate_slice(&program.instructions)?;
    Ok(program.measurement_count())
}

fn validate_slice(instructions: &[ExtendedInstruction]) -> Result<(), ExecError> {
    for instr in instructions {
        match instr {
            ExtendedInstruction::Gate { name, line, .. } => {
                check_gate_supported(*name, *line)?;
            }
            ExtendedInstruction::Noise { name, line, .. } => {
                check_noise_supported(*name, *line)?;
            }
            ExtendedInstruction::Measure { name, line, .. } => {
                check_measure_supported(*name, *line)?;
            }
            ExtendedInstruction::MPad { bits, line, .. } => {
                for (index, &value) in bits.iter().enumerate() {
                    if value > 1 {
                        return Err(ExecError::InvalidMPadBit {
                            line: *line,
                            index,
                            value,
                        });
                    }
                }
            }
            ExtendedInstruction::CorrelatedLoss { targets, line, .. } => {
                if targets.is_empty() || !targets.len().is_multiple_of(2) {
                    return Err(ExecError::InvalidCorrelatedLossArity {
                        line: *line,
                        found: targets.len(),
                    });
                }
            }
            ExtendedInstruction::Repeat { body, .. } => {
                validate_slice(body)?;
            }
            ExtendedInstruction::Annotation { .. }
            | ExtendedInstruction::T { .. }
            | ExtendedInstruction::TDag { .. }
            | ExtendedInstruction::Rotation { .. }
            | ExtendedInstruction::U3 { .. }
            | ExtendedInstruction::Loss { .. } => {}
            _ => unreachable!(
                "ExtendedInstruction variant added but not handled in prepare"
            ),
        }
    }
    Ok(())
}

fn check_gate_supported(name: GateName, line: usize) -> Result<(), ExecError> {
    use GateName::*;
    match name {
        Reset | ResetZ | X | Y | Z | H | HXZ | S | SqrtZ | SDag | SqrtZDag | SqrtX
        | SqrtXDag | SqrtY | SqrtYDag | Identity | CX | ZCX | CNot | CY | ZCY | CZ
        | ZCZ => Ok(()),
        Swap | ISwap | ISwapDag | SqrtXX | SqrtYY | SqrtZZ | CXSwap | SwapCX | XCX
        | XCY | XCZ | YCX | YCY | YCZ | CXYZ | CZYX | HXY | HYZ => {
            Err(ExecError::Unsupported {
                name: name.canonical_name().to_string(),
                line,
            })
        }
    }
}

fn check_noise_supported(name: NoiseName, line: usize) -> Result<(), ExecError> {
    use NoiseName::*;
    match name {
        Depolarize1 | Depolarize2 | PauliChannel1 | PauliChannel2 | XError | YError
        | ZError => Ok(()),
        IError | HeraldedErase | HeraldedPauliChannel1 | CorrelatedError
        | ElseCorrelatedError => Err(ExecError::Unsupported {
            name: name.canonical_name().to_string(),
            line,
        }),
    }
}

fn check_measure_supported(name: MeasureName, line: usize) -> Result<(), ExecError> {
    use MeasureName::*;
    match name {
        M | MZ | MR => Ok(()),
        other => Err(ExecError::Unsupported {
            name: other.canonical_name().to_string(),
            line,
        }),
    }
}
```

- [ ] **Step 4: Remove the old empty `ExecError` from `crates/ppvm-stim/src/executor.rs`**

In `executor.rs`, delete lines 14-16:

```rust
#[derive(Debug, thiserror::Error, PartialEq)]
#[non_exhaustive]
pub enum ExecError {}
```

Then add this import near the top of `executor.rs` (after the existing `use crate::tableau_program::…` line on line 12):

```rust
use crate::prepare::ExecError;
```

The `execute` and `sample` signatures still reference `ExecError` in their `Result<…, ExecError>` return type; the import resolves them to the new `prepare::ExecError` (whose new variants are extra-permissive — current code returns no `Err` from these functions yet, so behavior is unchanged).

- [ ] **Step 5: Update `crates/ppvm-stim/src/lib.rs` to expose `prepare` and re-export the new `ExecError`**

In `lib.rs`, find the existing module/use lines (lines 36-46):

```rust
pub mod executor;
pub mod normalize;
pub mod tableau_program;

pub use stim_parser::prelude::*;

pub use tableau_program::{GateKind, Instruction, MeasureKind, NoiseKind, TableauProgram};

pub use normalize::NormalizeError;

pub use executor::{ExecError, execute, sample};
```

Add `pub mod prepare;` right after `pub mod executor;` and change the `pub use executor::{ExecError, execute, sample};` line to `pub use executor::{execute, sample};`. Add `pub use prepare::{ExecError, prepare};` after the executor re-export. The block becomes:

```rust
pub mod executor;
pub mod normalize;
pub mod prepare;
pub mod tableau_program;

pub use stim_parser::prelude::*;

pub use tableau_program::{GateKind, Instruction, MeasureKind, NoiseKind, TableauProgram};

pub use normalize::NormalizeError;

pub use executor::{execute, sample};
pub use prepare::{ExecError, prepare};
```

- [ ] **Step 6: Run the new tests to confirm they pass**

Run: `cargo test -p ppvm-stim --test prepare`
Expected: PASS — all 6 tests green.

- [ ] **Step 7: Run the full ppvm-stim test suite to confirm no regression**

Run: `cargo test -p ppvm-stim`
Expected: PASS — the existing `normalize.rs` / `executor.rs` / `run.rs` / `stim_corpus.rs` tests are still wired through `normalize::to_tableau`, which still works.

- [ ] **Step 8: Format and lint**

Run:
```
cargo fmt -p ppvm-stim
cargo clippy -p ppvm-stim --all-targets -- -D warnings
```
Expected: both succeed with no diagnostics.

- [ ] **Step 9: Commit**

```bash
git add crates/ppvm-stim/src/prepare.rs crates/ppvm-stim/src/lib.rs crates/ppvm-stim/src/executor.rs crates/ppvm-stim/tests/prepare.rs
git commit -m "$(cat <<'EOF'
feat(ppvm-stim): add prepare() and ExecError variants

Adds the validate-and-count walk that will replace normalize::to_tableau.
Validation tests live in tests/prepare.rs. The old normalize path is
still wired up; Task 3 removes it.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Wire executor against `&ExtendedProgram`; remove old IR

This task is one logical change committed as a single unit. Intermediate file states between steps will not compile — that's expected. Verification (cargo fmt + clippy + test) runs only at the end before the commit.

**Files:**
- Modify: `crates/ppvm-stim/src/executor.rs` (rewrite)
- Modify: `crates/ppvm-stim/src/lib.rs` (drop modules + Error variant; update run_string)
- Modify: `crates/ppvm-stim/tests/executor.rs` (call sites)
- Modify: `crates/ppvm-stim/tests/run.rs` (test rename + assertion)
- Modify: `crates/ppvm-stim/tests/stim_corpus.rs` (Expect variant rename + imports)
- Modify: `crates/ppvm-stim/benches/tableau-msd-stim.rs` (signature + setup)
- Modify: `crates/ppvm-python-native/src/stim_program.rs` (PyStimProgram shape)
- Modify: `crates/ppvm-python-native/src/interface_tableau.rs` (call sites: `&prog.inner` → `&prog.program`)
- Delete: `crates/ppvm-stim/src/normalize.rs`
- Delete: `crates/ppvm-stim/src/tableau_program.rs`
- Delete: `crates/ppvm-stim/tests/normalize.rs`

- [ ] **Step 1: Rewrite `crates/ppvm-stim/src/executor.rs` to consume `&ExtendedProgram`**

Replace the entire file with the following content:

```rust
use bitvec::view::BitView;
use itertools::Itertools;
use num::Integer;
use num::PrimInt;
use num::complex::{Complex64, ComplexFloat};
use num::{Complex, One, ToPrimitive, Zero};
use std::fmt::Debug;

use ppvm_runtime::prelude::*;
use ppvm_tableau::prelude::*;
use stim_parser::ast::{GateName, MeasureName, NoiseName};
use stim_parser::extended::{Axis, ExtendedInstruction, ExtendedProgram};

use crate::prepare::{ExecError, prepare};

/// Execute a program against a tableau, returning the per-measurement
/// results in circuit order. Validates the program once via [`prepare`]
/// before dispatch.
pub fn execute<T, I, C>(
    program: &ExtendedProgram,
    tab: &mut GeneralizedTableau<T, I, C>,
) -> Result<Vec<Option<bool>>, ExecError>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Debug + Send + Sync,
{
    let count = prepare(program)?;
    let mut results = Vec::with_capacity(count);
    execute_slice(&program.instructions, tab, &mut results);
    Ok(results)
}

/// Execute many shots, building a fresh tableau per shot via `make_tableau`.
/// Validates the program once up front; per-shot execution is infallible.
pub fn sample<T, I, C, F>(
    program: &ExtendedProgram,
    num_shots: usize,
    mut make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Debug + Send + Sync,
    F: FnMut() -> GeneralizedTableau<T, I, C>,
{
    let count = prepare(program)?;
    Ok((0..num_shots)
        .map(|_| {
            let mut tab = make_tableau();
            let mut results = Vec::with_capacity(count);
            execute_slice(&program.instructions, &mut tab, &mut results);
            results
        })
        .collect())
}

fn execute_slice<T, I, C>(
    instructions: &[ExtendedInstruction],
    tab: &mut GeneralizedTableau<T, I, C>,
    results: &mut Vec<Option<bool>>,
) where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Debug + Send + Sync,
{
    for instr in instructions {
        match instr {
            ExtendedInstruction::Gate { name, targets, .. } => {
                use GateName::*;
                match name {
                    Reset | ResetZ => targets.iter().for_each(|&q| tab.reset(q)),
                    X => targets.iter().for_each(|&q| tab.x(q)),
                    Y => targets.iter().for_each(|&q| tab.y(q)),
                    Z => targets.iter().for_each(|&q| tab.z(q)),
                    H | HXZ => targets.iter().for_each(|&q| tab.h(q)),
                    S | SqrtZ => targets.iter().for_each(|&q| tab.s(q)),
                    SDag | SqrtZDag => targets.iter().for_each(|&q| tab.s_adj(q)),
                    SqrtX => targets.iter().for_each(|&q| tab.sqrt_x(q)),
                    SqrtXDag => targets.iter().for_each(|&q| tab.sqrt_x_adj(q)),
                    SqrtY => targets.iter().for_each(|&q| tab.sqrt_y(q)),
                    SqrtYDag => targets.iter().for_each(|&q| tab.sqrt_y_adj(q)),
                    Identity => {}
                    CX | ZCX | CNot => targets
                        .chunks_exact(2)
                        .for_each(|p| tab.cnot(p[0], p[1])),
                    CY | ZCY => targets
                        .chunks_exact(2)
                        .for_each(|p| tab.cy(p[0], p[1])),
                    CZ | ZCZ => targets
                        .chunks_exact(2)
                        .for_each(|p| tab.cz(p[0], p[1])),
                    Swap | ISwap | ISwapDag | SqrtXX | SqrtYY | SqrtZZ | CXSwap
                    | SwapCX | XCX | XCY | XCZ | YCX | YCY | YCZ | CXYZ | CZYX
                    | HXY | HYZ => unreachable!(
                        "unsupported gate {name:?} should have been rejected by prepare"
                    ),
                }
            }
            ExtendedInstruction::T { targets, .. } => {
                targets.iter().for_each(|&q| tab.t(q))
            }
            ExtendedInstruction::TDag { targets, .. } => {
                targets.iter().for_each(|&q| tab.t_adj(q))
            }
            ExtendedInstruction::Rotation {
                axis, theta, targets, ..
            } => match axis {
                Axis::X => targets.iter().for_each(|&q| tab.rx(q, *theta)),
                Axis::Y => targets.iter().for_each(|&q| tab.ry(q, *theta)),
                Axis::Z => targets.iter().for_each(|&q| tab.rz(q, *theta)),
                _ => unreachable!("Axis variant {axis:?} not handled in execute"),
            },
            ExtendedInstruction::U3 {
                theta,
                phi,
                lambda,
                targets,
                ..
            } => targets.iter().for_each(|&q| {
                tab.u3(q, (*theta).into(), (*phi).into(), (*lambda).into())
            }),
            ExtendedInstruction::Noise {
                name,
                args,
                targets,
                ..
            } => {
                use NoiseName::*;
                match name {
                    Depolarize1 => {
                        debug_assert_eq!(args.len(), 1);
                        let p = args[0];
                        for &q in targets {
                            tab.depolarize(q, p.into());
                        }
                    }
                    Depolarize2 => {
                        debug_assert_eq!(args.len(), 1);
                        let p = args[0];
                        for (a, b) in targets.iter().copied().tuples() {
                            tab.depolarize2(a, b, p.into());
                        }
                    }
                    PauliChannel1 => {
                        debug_assert_eq!(args.len(), 3);
                        let ps: [T::Coeff; 3] =
                            [args[0].into(), args[1].into(), args[2].into()];
                        for &q in targets {
                            tab.pauli_error(q, ps.clone());
                        }
                    }
                    PauliChannel2 => {
                        debug_assert_eq!(args.len(), 15);
                        let ps: [T::Coeff; 15] = std::array::from_fn(|i| args[i].into());
                        debug_assert!(targets.len().is_even());
                        for (a, b) in targets.iter().copied().tuples() {
                            tab.two_qubit_pauli_error(a, b, ps.clone());
                        }
                    }
                    XError => {
                        debug_assert_eq!(args.len(), 1);
                        let ps: [T::Coeff; 3] = [
                            args[0].into(),
                            T::Coeff::zero(),
                            T::Coeff::zero(),
                        ];
                        for &q in targets {
                            tab.pauli_error(q, ps.clone());
                        }
                    }
                    YError => {
                        debug_assert_eq!(args.len(), 1);
                        let ps: [T::Coeff; 3] = [
                            T::Coeff::zero(),
                            args[0].into(),
                            T::Coeff::zero(),
                        ];
                        for &q in targets {
                            tab.pauli_error(q, ps.clone());
                        }
                    }
                    ZError => {
                        debug_assert_eq!(args.len(), 1);
                        let ps: [T::Coeff; 3] = [
                            T::Coeff::zero(),
                            T::Coeff::zero(),
                            args[0].into(),
                        ];
                        for &q in targets {
                            tab.pauli_error(q, ps.clone());
                        }
                    }
                    IError | HeraldedErase | HeraldedPauliChannel1 | CorrelatedError
                    | ElseCorrelatedError => unreachable!(
                        "unsupported noise {name:?} should have been rejected by prepare"
                    ),
                }
            }
            ExtendedInstruction::Loss { p, targets, .. } => {
                for &q in targets {
                    tab.loss_channel(q, (*p).into());
                }
            }
            ExtendedInstruction::CorrelatedLoss { ps, targets, .. } => {
                let pps: [T::Coeff; 3] = [ps[0].into(), ps[1].into(), ps[2].into()];
                for (a, b) in targets.iter().copied().tuples() {
                    tab.correlated_loss_channel(a, b, pps.clone());
                }
            }
            ExtendedInstruction::Measure {
                name,
                args,
                targets,
                ..
            } => {
                use MeasureName::*;
                let noise = args.first().copied().unwrap_or(0.0);
                match name {
                    M | MZ => {
                        for &q in targets {
                            results.push(tab.measure_noisy(q, noise));
                        }
                    }
                    MR => {
                        for &q in targets {
                            // Use the true outcome to decide whether to reset, then
                            // apply measurement noise only to the *recorded* bit.
                            let true_outcome = tab.measure(q);
                            if true_outcome == Some(true) {
                                tab.x(q);
                            }
                            let recorded = match true_outcome {
                                Some(b) if noise > 0.0 && tab.bernoulli(noise) => Some(!b),
                                other => other,
                            };
                            results.push(recorded);
                        }
                    }
                    other => unreachable!(
                        "unsupported measure {other:?} should have been rejected by prepare"
                    ),
                }
            }
            ExtendedInstruction::MPad { bits, prob, .. } => {
                let noise = prob.unwrap_or(0.0);
                for &bit in bits {
                    let bit_bool = bit != 0; // safe: prepare validated bits ∈ {0,1}
                    let recorded = if noise > 0.0 && tab.bernoulli(noise) {
                        !bit_bool
                    } else {
                        bit_bool
                    };
                    results.push(Some(recorded));
                }
            }
            ExtendedInstruction::Annotation { .. } => { /* phase-1 no-op */ }
            ExtendedInstruction::Repeat { count, body, .. } => {
                for _ in 0..*count {
                    execute_slice(body, tab, results);
                }
            }
            _ => unreachable!(
                "ExtendedInstruction variant added but not handled in execute"
            ),
        }
    }
}
```

- [ ] **Step 2: Rewrite `crates/ppvm-stim/src/lib.rs`**

Replace the entire file with the following content:

```rust
//! Validate and execute Stim circuits against a [`GeneralizedTableau`].
//!
//! Two-stage pipeline:
//!
//! 1. [`parse_extended`] — `&str` → [`ExtendedProgram`] (re-exported from
//!    [`stim_parser`]).
//! 2. [`execute`] / [`sample`] — apply an [`ExtendedProgram`] to a
//!    [`GeneralizedTableau`]. Both call [`prepare`] internally to validate
//!    and pre-size the result vec.
//!
//! Multi-shot usage should call [`parse_extended`] once and call [`sample`]
//! for the shot loop. The [`run_string`] / [`run_file`] convenience helpers
//! re-parse on every call and are intended for single-shot demos only.
//!
//! # Multi-shot pattern (recommended)
//!
//! ```ignore
//! use ppvm_stim::{parse_extended, sample};
//! use ppvm_tableau::prelude::*;
//!
//! let prog = parse_extended(circuit_src)?;
//! let shots = sample(&prog, 10_000, || {
//!     GeneralizedTableau::<_, usize, _>::new(n_qubits, 1e-10)
//! })?;
//! # Ok::<(), ppvm_stim::Error>(())
//! ```
//!
//! [`run_string`] / [`run_file`] re-parse on every call and exist only for
//! single-shot demos — never call them from a shot loop.
//!
//! [`ExtendedProgram`]: stim_parser::extended::ExtendedProgram
//! [`GeneralizedTableau`]: ppvm_tableau::prelude::GeneralizedTableau

pub mod executor;
pub mod prepare;

pub use stim_parser::prelude::*;

pub use executor::{execute, sample};
pub use prepare::{ExecError, prepare};

use std::path::{Path, PathBuf};

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    #[error(transparent)]
    Parse(#[from] ExtendedParseError),
    #[error(transparent)]
    Exec(#[from] ExecError),
    #[error("failed to read stim file {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
}

/// Parse → execute in one shot. Re-parses each call; do **not** use in shot
/// loops — use [`parse_extended`] + [`sample`] instead.
pub fn run_string<T, I, C>(
    src: &str,
    tab: &mut ppvm_tableau::prelude::GeneralizedTableau<T, I, C>,
) -> Result<Vec<Option<bool>>, Error>
where
    T: ppvm_runtime::prelude::Config,
    <<T as ppvm_runtime::prelude::Config>::Storage as bitvec::view::BitView>::Store: num::PrimInt,
    C: ppvm_tableau::prelude::SparseVector<num::Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: num::One
        + num::Zero
        + Clone
        + num::Num
        + num::ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>
        + Send
        + Sync,
    num::Complex<T::Coeff>: std::ops::Mul<Output = num::Complex<T::Coeff>>
        + From<num::complex::Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + num::One
        + num::complex::ComplexFloat
        + Copy,
    I: ppvm_tableau::prelude::TableauIndex + std::fmt::Debug + Send + Sync,
{
    let prog = parse_extended(src)?;
    let results = execute(&prog, tab)?;
    Ok(results)
}

pub fn run_file<T, I, C>(
    path: &Path,
    tab: &mut ppvm_tableau::prelude::GeneralizedTableau<T, I, C>,
) -> Result<Vec<Option<bool>>, Error>
where
    T: ppvm_runtime::prelude::Config,
    <<T as ppvm_runtime::prelude::Config>::Storage as bitvec::view::BitView>::Store: num::PrimInt,
    C: ppvm_tableau::prelude::SparseVector<num::Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: num::One
        + num::Zero
        + Clone
        + num::Num
        + num::ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>
        + Send
        + Sync,
    num::Complex<T::Coeff>: std::ops::Mul<Output = num::Complex<T::Coeff>>
        + From<num::complex::Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + num::One
        + num::complex::ComplexFloat
        + Copy,
    I: ppvm_tableau::prelude::TableauIndex + std::fmt::Debug + Send + Sync,
{
    let src = std::fs::read_to_string(path).map_err(|source| Error::Io {
        path: path.to_path_buf(),
        source,
    })?;
    run_string(&src, tab)
}
```

- [ ] **Step 3: Delete `crates/ppvm-stim/src/normalize.rs`**

```bash
rm crates/ppvm-stim/src/normalize.rs
```

- [ ] **Step 4: Delete `crates/ppvm-stim/src/tableau_program.rs`**

```bash
rm crates/ppvm-stim/src/tableau_program.rs
```

- [ ] **Step 5: Delete `crates/ppvm-stim/tests/normalize.rs`**

```bash
rm crates/ppvm-stim/tests/normalize.rs
```

(All its content is either covered by `tests/prepare.rs`, by `tests/executor.rs`, or moved to `crates/stim-parser/tests/extended.rs` in Task 1.)

- [ ] **Step 6: Update `crates/ppvm-stim/tests/executor.rs`**

Replace the use line at the top (line 2):

```rust
use ppvm_stim::{execute, normalize, parse_extended};
```

with:

```rust
use ppvm_stim::{execute, parse_extended};
```

Then update the `run` helper at lines 7-13:

```rust
fn run(src: &str, n_qubits: usize) -> (Vec<Option<bool>>, Tab) {
    let prog = parse_extended(src).expect("parse_extended");
    let tprog = normalize::to_tableau(&prog).expect("normalize");
    let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
    let results = execute(&tprog, &mut tab).expect("execute");
    (results, tab)
}
```

becomes:

```rust
fn run(src: &str, n_qubits: usize) -> (Vec<Option<bool>>, Tab) {
    let prog = parse_extended(src).expect("parse_extended");
    let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
    let results = execute(&prog, &mut tab).expect("execute");
    (results, tab)
}
```

The same transform applies to every other test in the file that does the `parse_extended → normalize::to_tableau → execute` sequence. Replace each occurrence of:

```rust
let prog = parse_extended(SRC).unwrap();
let tprog = normalize::to_tableau(&prog).unwrap();
```

with:

```rust
let prog = parse_extended(SRC).unwrap();
```

and replace each `&tprog` with `&prog` in the following `execute(...)` / `sample(...)` call. There are about 10 such sites in this file. Specifically:

- `loss_channel_with_p1_marks_qubit_lost` (lines 76-83)
- `measurement_buffer_is_pre_sized` (lines 100-107)
- `sample_runs_n_shots_each_with_fresh_tableau` (lines 110-122)
- `sample_zero_shots_returns_empty` (lines 125-134)
- `sample_random_h_distribution_within_3_sigma` (lines 137-156)
- `measure_noise_distribution_within_3_sigma` (lines 383-403)
- `mpad_noise_distribution_within_3_sigma` (lines 439-458)

In `measurement_buffer_is_pre_sized`, change:

```rust
assert_eq!(tprog.expected_measurement_count, 5);
```

to:

```rust
assert_eq!(prog.measurement_count(), 5);
```

- [ ] **Step 7: Update `crates/ppvm-stim/tests/run.rs`**

Replace the use line at the top (line 2):

```rust
use ppvm_stim::{Error, ExtendedParseError, NormalizeError, ParseError, run_file, run_string};
```

with:

```rust
use ppvm_stim::{Error, ExecError, ExtendedParseError, ParseError, run_file, run_string};
```

Then rename the test `run_string_propagates_normalize_error` to `run_string_propagates_exec_error`, and change its body from:

```rust
#[test]
fn run_string_propagates_normalize_error() {
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-10);
    let err = run_string("SWAP 0 1", &mut tab).unwrap_err();
    assert!(matches!(
        err,
        Error::Normalize(NormalizeError::Unsupported { .. })
    ));
}
```

to:

```rust
#[test]
fn run_string_propagates_exec_error() {
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-10);
    let err = run_string("SWAP 0 1", &mut tab).unwrap_err();
    assert!(matches!(
        err,
        Error::Exec(ExecError::Unsupported { .. })
    ));
}
```

- [ ] **Step 8: Update `crates/ppvm-stim/tests/stim_corpus.rs`**

Replace lines 4 and lines 7-15 of the file.

Replace the use line:

```rust
use ppvm_stim::{NormalizeError, execute, normalize, parse_extended};
```

with:

```rust
use ppvm_stim::{ExecError, execute, parse_extended};
```

Replace the `Expect` enum:

```rust
#[derive(Debug, Clone, Copy)]
enum Expect {
    /// File parses, normalizes, executes.
    Ok,
    /// File parses, but normalize must fail with `Unsupported(name)`.
    NormalizeUnsupported(&'static str),
    /// File should fail at parse time (e.g. uses `rec[-k]` targets).
    ParseFails,
}
```

with:

```rust
#[derive(Debug, Clone, Copy)]
enum Expect {
    /// File parses and executes.
    Ok,
    /// File parses, but execute must fail with `Unsupported(name)` from prepare.
    ExecUnsupported(&'static str),
    /// File should fail at parse time (e.g. uses `rec[-k]` targets).
    ParseFails,
}
```

Update the two table entries on lines 24 and 27:

```rust
    (
        "swap_unsupported.stim",
        Expect::NormalizeUnsupported("SWAP"),
    ),
    ("mx_unsupported.stim", Expect::NormalizeUnsupported("MX")),
```

become:

```rust
    (
        "swap_unsupported.stim",
        Expect::ExecUnsupported("SWAP"),
    ),
    ("mx_unsupported.stim", Expect::ExecUnsupported("MX")),
```

Replace the `corpus_obeys_expectations` test body (lines 73-103). The new body:

```rust
#[test]
fn corpus_obeys_expectations() {
    type Tab = GeneralizedTableau<ByteFxHashF64<8>, usize>;

    for (name, expect) in CASES {
        let src = read(name);
        let parsed = parse_extended(&src);
        match (expect, parsed) {
            (Expect::ParseFails, Ok(_)) => {
                panic!("{name}: expected parse failure, but parse succeeded");
            }
            (Expect::ParseFails, Err(_)) => continue,
            (Expect::Ok, Err(e)) | (Expect::ExecUnsupported(_), Err(e)) => {
                panic!("{name}: parse failed unexpectedly: {e}");
            }
            (Expect::Ok, Ok(prog)) => {
                let mut tab: Tab = GeneralizedTableau::new(64, 1e-10);
                execute(&prog, &mut tab)
                    .unwrap_or_else(|e| panic!("{name}: execute failed: {e}"));
            }
            (Expect::ExecUnsupported(expected_name), Ok(prog)) => {
                let mut tab: Tab = GeneralizedTableau::new(64, 1e-10);
                match execute(&prog, &mut tab) {
                    Err(ExecError::Unsupported { name: n, .. }) => {
                        assert_eq!(n, *expected_name, "{name}: wrong unsupported name");
                    }
                    Err(other) => panic!("{name}: expected Unsupported, got {other:?}"),
                    Ok(_) => panic!("{name}: expected Unsupported, but execute succeeded"),
                }
            }
        }
    }
}
```

- [ ] **Step 9: Update `crates/ppvm-stim/benches/tableau-msd-stim.rs`**

Replace line 5:

```rust
use ppvm_stim::{execute, normalize, parse_extended};
```

with:

```rust
use ppvm_stim::{execute, parse_extended};
use stim_parser::extended::ExtendedProgram;
```

Replace `fn msd_stim_func` (lines 124-128):

```rust
fn msd_stim_func(prog: &ppvm_stim::TableauProgram) {
    let n_qubits = 17 * 5;
    let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
    execute(prog, &mut tab).expect("execute");
}
```

with:

```rust
fn msd_stim_func(prog: &ExtendedProgram) {
    let n_qubits = 17 * 5;
    let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
    execute(prog, &mut tab).expect("execute");
}
```

Replace the bench setup (lines 131-133):

```rust
    let circuit = msd_stim_string();
    let parsed = parse_extended(&circuit).expect("parse_extended");
    let prog = normalize::to_tableau(&parsed).expect("normalize");
```

with:

```rust
    let circuit = msd_stim_string();
    let prog = parse_extended(&circuit).expect("parse_extended");
```

Note: `stim-parser` is already a workspace member; `crates/ppvm-stim/Cargo.toml` already has `stim-parser` under `[dependencies]`, so no manifest change is needed for the new `use stim_parser::extended::ExtendedProgram;` line.

- [ ] **Step 10: Update `crates/ppvm-python-native/src/stim_program.rs`**

Replace the entire file with the following content:

```rust
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;

use ppvm_stim::{ExecError, ExtendedParseError, parse_extended, prepare};
use stim_parser::extended::ExtendedProgram;

/// Python-facing wrapper around a parsed and validated Stim program.
#[pyclass(name = "StimProgram", module = "ppvm_python_native")]
pub struct PyStimProgram {
    pub(crate) program: ExtendedProgram,
    pub(crate) measurement_count: usize,
}

#[pymethods]
impl PyStimProgram {
    /// Parse and validate a Stim circuit string. Validation errors surface here.
    #[staticmethod]
    pub fn parse(src: &str) -> PyResult<Self> {
        let program = parse_extended(src).map_err(stim_to_pyerr_parse)?;
        let measurement_count = prepare(&program).map_err(stim_to_pyerr_exec)?;
        Ok(Self {
            program,
            measurement_count,
        })
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
            self.program.instructions.len(),
            self.measurement_count
        )
    }
}

fn stim_to_pyerr_parse(e: ExtendedParseError) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}

fn stim_to_pyerr_exec(e: ExecError) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}
```

Note: `crates/ppvm-python-native/Cargo.toml` may or may not already have `stim-parser` as a direct dependency. Check by reading the manifest.

- [ ] **Step 11: Confirm `stim-parser` is a dependency of `ppvm-python-native`**

Run: `grep -n stim-parser crates/ppvm-python-native/Cargo.toml`

If the output shows a `stim-parser = …` line under `[dependencies]`, no change needed.

If the output is empty (stim-parser is not a direct dep), add it. Open `crates/ppvm-python-native/Cargo.toml` and under `[dependencies]` add a line matching the pattern of the existing `ppvm-stim` line (which today reads `ppvm-stim = { version = "0.1.0", path = "../ppvm-stim" }`):

```toml
stim-parser = { version = "0.1.0", path = "../stim-parser" }
```

- [ ] **Step 12: Update `crates/ppvm-python-native/src/interface_tableau.rs`**

Replace `&prog.inner` with `&prog.program` at the two call sites:

Line 175:

```rust
                ppvm_stim::execute(&prog.inner, &mut self.inner)
```

becomes:

```rust
                ppvm_stim::execute(&prog.program, &mut self.inner)
```

Line 190:

```rust
                ppvm_stim::sample(&prog.inner, num_shots, || {
```

becomes:

```rust
                ppvm_stim::sample(&prog.program, num_shots, || {
```

- [ ] **Step 13: Run `cargo fmt` across the workspace**

Run: `cargo fmt --all`
Expected: succeeds with no diagnostics.

- [ ] **Step 14: Run `cargo clippy` with warnings-as-errors**

Run: `cargo clippy --workspace --all-targets -- -D warnings`
Expected: PASS — no warnings, no errors. If clippy complains about the unused fields `tags` / `args` / `prob` in some `ExtendedInstruction` patterns where the executor uses `..`, leave them — `..` already discards the unused fields. Any new clippy warning that is not pre-existing must be fixed before proceeding.

- [ ] **Step 15: Run the full workspace test suite**

Run: `cargo test --workspace`
Expected: PASS — all tests across `stim-parser`, `ppvm-stim`, `ppvm-python-native`, `ppvm-tableau`, etc.

Specifically check that these test files all pass:
- `cargo test -p stim-parser --test extended` (Task 1's three count tests)
- `cargo test -p ppvm-stim --test prepare` (Task 2's six validation tests)
- `cargo test -p ppvm-stim --test executor` (the migrated executor tests)
- `cargo test -p ppvm-stim --test run` (with renamed `run_string_propagates_exec_error`)
- `cargo test -p ppvm-stim --test stim_corpus` (with renamed `Expect::ExecUnsupported`)

- [ ] **Step 16: Sanity-check the bench compiles**

Run: `cargo bench -p ppvm-stim --bench tableau-msd-stim --no-run`
Expected: compiles successfully (we don't need to run the bench, just confirm it builds).

- [ ] **Step 17: Commit**

Stage every modified or deleted file and commit:

```bash
git add -A crates/ppvm-stim/src crates/ppvm-stim/tests crates/ppvm-stim/benches crates/ppvm-python-native/src crates/ppvm-python-native/Cargo.toml
git status   # confirm normalize.rs / tableau_program.rs / tests/normalize.rs are staged for deletion
git commit -m "$(cat <<'EOF'
refactor(ppvm-stim): collapse TableauProgram IR into ExtendedProgram

Executor now takes &ExtendedProgram directly; prepare() validates
and returns the measurement count up front. Deletes normalize.rs
(282 lines), tableau_program.rs (95 lines), and the redundant
tests/normalize.rs. Updates all consumers (executor tests, run
tests, stim_corpus, bench, ppvm-python-native).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

(If the `git add -A` glob misses the deleted `Cargo.toml` change because it wasn't modified, that's fine — only stage it if Step 11 changed it.)

---

## Self-Review

**Spec coverage** (skimmed each section against the tasks):

- Goal / pipeline shape: Task 3 Step 1 (executor.rs) + Step 2 (lib.rs run_string) implement.
- `ExtendedProgram::measurement_count` (stim-parser): Task 1.
- `prepare` + `ExecError` (ppvm-stim): Task 2.
- Backend-capability filter location (ppvm-stim, not stim-parser): Task 2 Step 3 (`check_*_supported` in `prepare.rs`).
- Validation rules (Unsupported / InvalidMPadBit / InvalidCorrelatedLossArity): Task 2 Step 3.
- Executor dispatch table with alias unions: Task 3 Step 1.
- Public `Error` drops `Normalize`: Task 3 Step 2.
- `run_string`/`run_file` keep signatures, body becomes parse → execute: Task 3 Step 2.
- Test migration (delete redundant, rename `tests/normalize.rs` → `tests/prepare.rs`, move 3 to stim-parser): Tasks 1, 2 (create), 3 Steps 5/6/7/8 (delete + update).
- Bench update: Task 3 Step 9.
- ppvm-python-native (`PyStimProgram` shape, two call sites): Task 3 Steps 10–12.
- Verification gate (cargo fmt / clippy / test): Task 3 Steps 13–15.
- Phase-2 deferral (AST tightening): out of plan scope, called out in spec.

**Placeholder scan:** No "TBD" / "TODO" / "implement later" / vague handling notes anywhere. Every code step shows the exact code; every test step shows the exact assertion. The clippy step explicitly tells the engineer how to handle the `tags`/`args`/`prob` discard pattern.

**Type consistency:**
- `ExecError` variants (`Unsupported`, `InvalidMPadBit`, `InvalidCorrelatedLossArity`) are introduced in Task 2 Step 3 and consumed identically in Task 2 Step 1 (tests), Task 3 Step 7 (run.rs), Task 3 Step 8 (stim_corpus).
- `prepare` signature `fn prepare(&ExtendedProgram) -> Result<usize, ExecError>` is consistent across Task 2 (declaration) and Task 3 (callers in `executor.rs` and `stim_program.rs`).
- `ExtendedProgram::measurement_count` returns `usize` everywhere (Task 1 declaration, Task 3 Step 6 caller in `executor.rs` test, Task 3 Step 10 Python `__repr__`).
- `PyStimProgram` field rename `inner` → `program` is consistent between Task 3 Step 10 (declaration) and Task 3 Step 12 (call sites).
