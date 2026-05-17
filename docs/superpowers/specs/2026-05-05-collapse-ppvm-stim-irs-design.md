# Collapse `ppvm-stim`'s `TableauProgram` IR into the parser AST

**Date:** 2026-05-05
**Crates touched:** `ppvm-stim`, `stim-parser`, `ppvm-python-native`

## Goal

Eliminate the intermediate `TableauProgram` IR in `ppvm-stim`. The executor consumes `stim_parser::extended::ExtendedProgram` directly. Two up-front pre-walks — `ExtendedProgram::measurement_count` (in `stim-parser`, pure AST) and `ppvm_stim::prepare` (backend-capability validation) — replace the old normalize pass; both run once per program before any shot loop, returning the measurement count for `Vec::with_capacity` pre-sizing.

## Why

`ppvm-stim/src/normalize.rs` (282 lines) and `ppvm-stim/src/tableau_program.rs` (95 lines) define an intermediate IR that, after the recent extended-dialect refactor (`49cb98e`, `a2c283b`), is mostly a structural copy of `ExtendedProgram` with three filters applied:

1. Alias collapsing (e.g. `H | HXZ → GateKind::H`).
2. Backend-capability filter (rejecting `Swap`, `XCY`, `MX`, `HERALDED_ERASE`, …).
3. Three structural validations: `MPad` bits ∈ {0, 1}, `CorrelatedLoss` target count nonzero and even, plus the count walk for `expected_measurement_count`.

(1) is expressible at the executor's match site via `match name { GateName::H | GateName::HXZ => … }`. (2) and (3) are cheap up-front walks. The IR layer is not earning its keep — collapsing it removes ≈280 net lines, eliminates one of the two error types, and leaves `ppvm-stim` with three source files instead of four.

## Non-goals (deferred to phase 2)

- **AST tightening in `stim-parser`.** Moving `ExtendedInstruction::MPad::bits` from `Vec<usize>` to `Vec<bool>` and `CorrelatedLoss::targets` from `Vec<usize>` to `Vec<(usize, usize)>` is structurally cleaner and would push validations (1)+(2) above into the parser. Out of scope here. Phase 2 will handle it.
- **Moving the backend-capability filter (unsupported gates) into `stim-parser`.** This stays in `ppvm-stim` permanently — it encodes ppvm-tableau's capability set, not a property of the Stim language. A future state-vector backend would have a different list. `stim-parser` deliberately stays backend-agnostic.

## Architecture

### Two new public surfaces

**`stim_parser::extended::ExtendedProgram::measurement_count`** — pure AST utility:

```rust
impl ExtendedProgram {
    /// Total number of recorded bits the program will produce, accounting
    /// for `REPEAT` factors. Pure AST property; safe for any backend to use.
    pub fn measurement_count(&self) -> usize { … }
}
```

Walks once, sums `Measure { targets }.len()` and `MPad { bits }.len()` per occurrence multiplied by enclosing `Repeat::count` factors. Always succeeds.

**`ppvm_stim::prepare`** — backend-capability + structural validation:

```rust
pub fn prepare(program: &ExtendedProgram) -> Result<usize, ExecError> {
    validate_slice(&program.instructions)?;
    Ok(program.measurement_count())
}
```

Walks the program once to enforce `ppvm-tableau`'s contract; on success returns the cached count. Two walks total per program (one for count via the parser helper, one for validation in `ppvm-stim`). Programs are typically hundreds of instructions; this runs once before a shot loop, never per shot.

### Deleted

- `crates/ppvm-stim/src/tableau_program.rs` (95 lines): `TableauProgram`, `Instruction`, `GateKind`, `NoiseKind`, `MeasureKind`.
- `crates/ppvm-stim/src/normalize.rs` (282 lines): `to_tableau`, `NormalizeError`, all helpers.

### Added

- `crates/stim-parser/src/extended/ast.rs`: a fresh `impl ExtendedProgram { pub fn measurement_count(&self) -> usize }` block appended after the existing struct/enum definitions, plus a private recursive helper (`fn count_in_slice(instructions: &[ExtendedInstruction], factor: u64) -> usize`) at module scope.
- `crates/ppvm-stim/src/prepare.rs`: `pub fn prepare(&ExtendedProgram) -> Result<usize, ExecError>`, plus the `ExecError` enum and the `validate_slice` recursive helper.

### Modified

- `crates/ppvm-stim/src/executor.rs`: `execute` and `sample` take `&ExtendedProgram`; both call `prepare(program)?` at entry. `execute_slice` matches on `ExtendedInstruction` directly; `Gate`/`Noise`/`Measure` arms inner-match on `GateName`/`NoiseName`/`MeasureName` with alias unions; promoted variants (`T`, `TDag`, `Rotation`, `U3`, `Loss`, `CorrelatedLoss`, `MPad`) dispatch directly. `execute_slice` becomes infallible — `unreachable!()` arms guard the prepare contract for unsupported names.
- `crates/ppvm-stim/src/lib.rs`: drops `pub mod normalize` and `pub mod tableau_program`; adds `pub mod prepare`. Re-exports `prepare` and `ExecError`. `pub enum Error` drops the `Normalize(NormalizeError)` variant; final shape is `Parse(ExtendedParseError) | Exec(ExecError) | Io { … }`. `run_string` body becomes `parse_extended(src)?; execute(&prog, tab)?` — the normalize step is gone.

## Error types

`ExecError` (not `#[non_exhaustive]` — pre-1.0, workspace-internal, exhaustive matches catch new variants at compile time):

```rust
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum ExecError {
    #[error("unsupported instruction '{name}' at line {line}")]
    Unsupported { name: String, line: usize },

    #[error("'MPAD' at line {line} target #{index} = {value}, must be 0 or 1")]
    InvalidMPadBit { line: usize, index: usize, value: usize },

    #[error(
        "'I_ERROR[correlated_loss]' at line {line} expected nonzero target count divisible by 2, got {found}"
    )]
    InvalidCorrelatedLossArity { line: usize, found: usize },
}
```

The three variants are 1:1 with today's `NormalizeError::{Unsupported, InvalidMPadTarget, InvalidCorrelatedLossTargetCount}`; only the names change (and `InvalidMPadTarget` → `InvalidMPadBit` for clarity, since the field is bit value not target index).

## Validation rules (in `ppvm-stim::prepare`)

`validate_slice` recurses into `Repeat::body` and runs:

| `ExtendedInstruction` variant | Rejection rule |
|---|---|
| `Gate { name, line, .. }` | reject if `name ∈ {Swap, ISwap, ISwapDag, SqrtXX, SqrtYY, SqrtZZ, CXSwap, SwapCX, XCX, XCY, XCZ, YCX, YCY, YCZ, CXYZ, CZYX, HXY, HYZ}` → `ExecError::Unsupported` |
| `Noise { name, line, .. }` | reject if `name ∈ {IError, HeraldedErase, HeraldedPauliChannel1, CorrelatedError, ElseCorrelatedError}` → `ExecError::Unsupported` |
| `Measure { name, line, .. }` | reject if `name ∉ {M, MZ, MR}` → `ExecError::Unsupported` |
| `MPad { bits, line, .. }` | for each `(index, value)`, reject if `value > 1` → `ExecError::InvalidMPadBit` |
| `CorrelatedLoss { targets, line, .. }` | reject if `targets.is_empty() || !targets.len().is_multiple_of(2)` → `ExecError::InvalidCorrelatedLossArity` |
| `Repeat { body, .. }` | recurse into `body` |
| `T`, `TDag`, `Rotation`, `U3`, `Loss`, `Annotation` | listed explicitly with no checks |
| `_` | `unreachable!("ExtendedInstruction variant added but not handled in prepare")` |

The supported/unsupported gate sets are 1:1 with today's `vanilla_gate_to_kind`/`vanilla_noise_to_kind`/`measure_to_kind` in `normalize.rs`.

## Executor dispatch (in `ppvm-stim::executor`)

Top-level match on `ExtendedInstruction`. Concrete dispatch:

| Variant | Dispatch |
|---|---|
| `Gate { name, targets, .. }` | inner match on `GateName`. Supported set unioned by alias: `Reset \| ResetZ → tab.reset(q)`; `H \| HXZ → tab.h(q)`; `S \| SqrtZ → tab.s(q)`; `SDag \| SqrtZDag → tab.s_adj(q)`; `CX \| ZCX \| CNot → tab.cnot(c, t)`; `CY \| ZCY → tab.cy(c, t)`; `CZ \| ZCZ → tab.cz(c, t)`. Singletons: `X`, `Y`, `Z`, `SqrtX`, `SqrtXDag`, `SqrtY`, `SqrtYDag`. `Identity → ()`. Unsupported names fall through to `unreachable!()`. |
| `Noise { name, args, targets, .. }` | inner match on `NoiseName` for `Depolarize1`, `Depolarize2`, `PauliChannel1`, `PauliChannel2`, `XError`, `YError`, `ZError`. Argument extraction (`debug_assert_eq!(args.len(), …)`) and tableau dispatch identical to today's executor. Unsupported → `unreachable!()`. |
| `Measure { name, args, targets, .. }` | extract `noise = args.first().copied().unwrap_or(0.0)`. Inner match: `M \| MZ → tab.measure_noisy(q, noise)`; `MR → measure-then-reset-then-noise-on-recorded-bit` (logic identical to today's executor). Other names → `unreachable!()`. |
| `T { targets, .. }` | `targets.iter().for_each(\|&q\| tab.t(q))` |
| `TDag { targets, .. }` | `targets.iter().for_each(\|&q\| tab.t_adj(q))` |
| `Rotation { axis, theta, targets, .. }` | match on `Axis` → `tab.rx/ry/rz(q, theta)` |
| `U3 { theta, phi, lambda, targets, .. }` | `tab.u3(q, theta, phi, lambda)` |
| `Loss { p, targets, .. }` | `tab.loss_channel(q, p)` |
| `CorrelatedLoss { ps, targets, .. }` | `for (a, b) in tuples → tab.correlated_loss_channel(a, b, ps)` |
| `MPad { bits, prob, .. }` | `noise = prob.unwrap_or(0.0)`; for each `&bit in bits`, treat as `bit != 0` (safe — `prepare` validated ∈ {0,1}); apply readout noise; push `Some(recorded)` |
| `Annotation { .. }` | no-op |
| `Repeat { count, body, .. }` | recurse `count` times |
| `_` | `unreachable!("ExtendedInstruction variant added but not handled in execute")` |

`execute_slice` is infallible (`-> ()`). `execute` and `sample` validate via `prepare` once and then loop infallibly.

### Public signatures

```rust
pub fn execute<T, I, C>(
    program: &ExtendedProgram,
    tab: &mut GeneralizedTableau<T, I, C>,
) -> Result<Vec<Option<bool>>, ExecError>
where
    /* same trait bounds as today */;

pub fn sample<T, I, C, F>(
    program: &ExtendedProgram,
    num_shots: usize,
    mut make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
where
    /* same trait bounds as today */,
    F: FnMut() -> GeneralizedTableau<T, I, C>;
```

`run_string` and `run_file` keep their public signatures unchanged.

## Test migration

### `ppvm-stim/tests/normalize.rs` → `ppvm-stim/tests/prepare.rs` (≈6 tests)

**Keep** (reframed against `prepare(&prog)` returning `ExecError`):
- `unsupported_swap_rejected`
- `unsupported_mx_rejected`
- `unsupported_heralded_erase_rejected`
- `correlated_loss_single_target_rejected`
- `correlated_loss_odd_targets_rejected`
- `mpad_target_two_rejected`

**Delete** (already covered behaviorally by `tests/executor.rs`):
- `promoted_extended_variants_map_directly_to_tableau_instructions` — the IR being mapped to no longer exists; behavioral coverage in `rx_pi_flips_qubit`, `u3_pi_flip_via_y_axis`, `t_gate_via_s_t_tag_no_op_on_zero`, `loss_channel_with_p1_marks_qubit_lost`.
- `h_maps_to_gate_h`, `cnot_alias_maps_to_cx`, `h_xz_alias_maps_to_h`, `sqrt_z_alias_maps_to_s`, `r_and_rz_both_map_to_reset` — alias-collapse coverage in `cnot_alias_equivalents`, `test_stim_zcx_alias`, `test_stim_zcy_alias`, `test_stim_zcz_alias`, `test_stim_sqrt_z_is_s`, `test_stim_sqrt_z_dag_is_s_adj`.
- `x_error_y_error_z_error_supported`, `measurements_m_mz_map_to_m`, `measurement_mr_maps_to_mr`, `annotations_become_no_op_annotations` — IR translation tests with no IR.
- `measure_noise_arg_passes_through_normalize`, `measure_no_noise_arg_defaults_to_zero`, `mr_noise_passes_through` — covered by `measure_noise_*` family in `executor.rs`.
- `mpad_normalize_zero_one_succeeds`, `mpad_normalize_with_prob` — covered by `mpad_*` family in `executor.rs`.

**Move to `crates/stim-parser/tests/extended.rs`** (against `prog.measurement_count()`):
- `expected_measurement_count_counts_m_mz_mr`
- `expected_measurement_count_includes_repeat_multiplier`
- `mpad_inside_repeat_block_multiplies_count`

### `ppvm-stim/tests/executor.rs`

Mechanical: drop `normalize` from imports; the `let tprog = normalize::to_tableau(&prog).unwrap()` line goes away; `execute(&tprog, …)` / `sample(&tprog, …)` become `execute(&prog, …)` / `sample(&prog, …)`. The `measurement_buffer_is_pre_sized` assertion changes from `tprog.expected_measurement_count` to `prog.measurement_count()`.

### `ppvm-stim/tests/run.rs`

Rename `run_string_propagates_normalize_error` → `run_string_propagates_exec_error`. Update the assertion from `Error::Normalize(NormalizeError::Unsupported { .. })` to `Error::Exec(ExecError::Unsupported { .. })`.

### `ppvm-stim/tests/stim_corpus.rs`

Drop `normalize`/`NormalizeError` from imports; replace with `ExecError`. Update call sites mechanically.

### `ppvm-stim/benches/tableau-msd-stim.rs`

`fn msd_stim_func(prog: &TableauProgram)` → `fn msd_stim_func(prog: &ExtendedProgram)`. Drop the `normalize::to_tableau` call site.

## `ppvm-python-native` migration

### `crates/ppvm-python-native/src/stim_program.rs`

```rust
use ppvm_stim::{ExecError, ExtendedParseError, parse_extended, prepare};
use stim_parser::extended::ExtendedProgram;

#[pyclass(name = "StimProgram", module = "ppvm_python_native")]
pub struct PyStimProgram {
    pub(crate) program: ExtendedProgram,
    pub(crate) measurement_count: usize,
}

#[pymethods]
impl PyStimProgram {
    #[staticmethod]
    pub fn parse(src: &str) -> PyResult<Self> {
        let program = parse_extended(src).map_err(stim_to_pyerr_parse)?;
        let measurement_count = prepare(&program).map_err(stim_to_pyerr_exec)?;
        Ok(Self { program, measurement_count })
    }

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

fn stim_to_pyerr_parse(e: ExtendedParseError) -> PyErr { PyValueError::new_err(format!("{e}")) }
fn stim_to_pyerr_exec(e: ExecError) -> PyErr { PyValueError::new_err(format!("{e}")) }
```

Validation errors continue to surface at `StimProgram.parse(...)` time (matches today's UX where `to_tableau` ran at parse).

### `crates/ppvm-python-native/src/interface_tableau.rs`

`ppvm_stim::execute(&prog.inner, …)` → `ppvm_stim::execute(&prog.program, …)`. Same for `sample`. The `inner` field name is replaced by `program`.

## Performance

The hot path (per-instruction gate dispatch) goes through one match layer (`ExtendedInstruction::Gate` outer + `GateName` inner) instead of two (`Instruction::Gate` outer + `GateKind` inner). Same number of comparisons; same inlining behavior. No regression expected.

`ExtendedInstruction::Gate` carries two extra empty `Vec`s (`tags`, `args`) that the executor ignores — 48 bytes per `Gate` instruction more than today's `Instruction::Gate`. For million-instruction programs, that's tens of MB of extra static program memory. For typical workloads, negligible.

`tableau-msd-stim` bench is the regression check. Run it before and after; expected delta within noise.

## Verification gate

Before committing the implementation:

```
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

All three must pass.

## Phase-2 follow-up (separate change)

After this lands:

1. Tighten `stim_parser::extended::ExtendedInstruction::MPad::bits` to `Vec<bool>`. Validation moves into the parser; `ExecError::InvalidMPadBit` goes away.
2. Tighten `ExtendedInstruction::CorrelatedLoss::targets` to `Vec<(usize, usize)>` (or assert arity at parse construction). `ExecError::InvalidCorrelatedLossArity` goes away.
3. After (1) and (2), `ExecError` collapses to a single `Unsupported` variant.

Phase 2 is structurally nicer but ripples through stim-parser AST + parser code + parser tests. Bundling it here would triple the diff. Leave for a follow-up PR.
