# Plan: Python Bindings for ppvm-timeevolve (Task 31+)

## Context

The time evolution solver (`ppvm-timeevolve`) is currently Rust-only. This plan adds Python bindings so users can drive Lindblad ODE solves directly from Python — constructing Hamiltonians and dissipators with the existing `PauliSum` API, calling `solve()`, and obtaining state snapshots or scalar observables.

This work spans three packages:
- `crates/ppvm-python-native` — Rust PyO3 bindings (new `interface_timeevolve.rs`, Cargo.toml update, lib.rs update)
- `ppvm-python` — Python wrapper package (new `timeevolve.py`, `__init__.py` update, new tests)
- `crates/ppvm-timeevolve/plans/` — a copy of this plan is written there at implementation start

The developer guidelines from `agents/perf-developer.md` apply: allocations are suspects, prefer clear+reuse over new, no `unwrap()` in production, short focused functions, `pub(crate)` for non-public internals.

---

## Key Design Decisions

### 1. Observable / Callback Design

**Two mutually exclusive modes. Full state cloning is disabled when a scalar observable is specified.**

| Option | Verdict |
|---|---|
| Python callable passed to Rust | **Never** — GIL re-entry per save point, dangling reference risk |
| State snapshots | **Default** — `observable=None` returns `list[PauliSum]` |
| Pre-defined scalar trace | **Fast path** — `observable="trace:<pattern>"` returns `list[float]`, no clone |

Two native Rust functions (not three — multi-observable is unified):
- `solve_timeevolve_states(...)` → `(list[float], list[PauliSumNative])` — clone callback
- `solve_timeevolve_observables(state, ..., patterns: list[str])` → `(list[float], list[list[float]])` — Rust computes all traces per save point in a single pass, no cloning

Python `solve()` dispatches based on `observable` kwarg:
- `observable=None` → `solve_timeevolve_states`
- `observable="trace:Z0"` → `solve_timeevolve_observables(patterns=["Z0"])`, unwraps inner list → `list[float]`
- `observable=["trace:Z0", "trace:Z1"]` → `solve_timeevolve_observables(patterns=["Z0", "Z1"])` → `list[list[float]]`

> **Phase 2**: other observable types (`overlap`, `n_terms`, etc.)

### 2. Generic Type Handling

**Single Python `solve()` entry point, 16-arm runtime dispatch in each Rust function.**

- Python `solve()` extracts `state._interface` (raw Rust-native object) and passes it to the native function
- Rust function does `downcast` on the concrete type → calls monomorphized `ppvm_timeevolve::solve::<ConcreteType, _, _>`
- Hamiltonian (if provided) must downcast to the **same** concrete type; mismatch → `TypeError` with clear message
- **Phase 1**: non-loss types only (N = 0..15, 16 arms).

### 3. LindbladOp API

**Pure Python data container. Rust `LindbladOp<T>` built inside native solve function, not at construction.**

- `LadderOp` — plain Python dataclass: `LadderOp(qubit: int, direction: "raise"|"lower")`
- `LindbladOp` — plain Python dataclass: `LindbladOp(jump_ops: list[LadderOp], rates: list[float] | list[list[float]])`
- Native function receives `Vec<(usize, String)>` (qubit, direction) + rates as `PyAny`
- Rates dispatched in Rust: `list[float]` → `RateMatrix::Vector`, `list[list[float]]` → `RateMatrix::Dense`
- **Phase 1**: ladder operators only. `CollapseOp` deferred to Phase 2.
- **Performance note**: `LindbladOp<T>::new` runs O(n²) expansion every call. Memoization is Phase 2.

### 4. SolverConfig

Plain Python dataclass matching Rust `SolverConfig` exactly. Individual fields passed as primitives to native function.

---

## Python-Facing API

```python
from ppvm.paulisum import PauliSum
from ppvm.timeevolve import LadderOp, LindbladOp, SolverConfig, solve

ham = PauliSum.new(4, [("ZIII", 0.5), ("IZII", 0.5)])
lindblad = LindbladOp(
    jump_ops=[LadderOp(qubit=0, direction="lower"),
              LadderOp(qubit=1, direction="lower")],
    rates=[1.0, 1.0],               # diagonal: one rate per jump op
    # rates=[[g11,g12],[g21,g22]]   # or dense rate matrix
)
state = PauliSum.new(4, "ZIII")

# Mode 1: state snapshots (default, observable=None)
times, states = solve(
    state=state, lindblad=lindblad, t_span=(0.0, 5.0),
    save_at=[1.0, 2.0, 3.0, 5.0],
    hamiltonian=ham,         # optional, default None
    config=SolverConfig(),   # optional
)
# states: list[PauliSum] — fully usable evolved states
states[0].trace("ZIZI")

# Mode 2: single scalar observable (no state clone, fast path)
# Pattern is passed directly to p.trace(pattern) — full Pauli string with ? and * wildcards
times, values = solve(
    state=state, lindblad=lindblad, t_span=(0.0, 5.0),
    save_at=[1.0, 2.0, 3.0, 5.0],
    observable="trace:ZIII",   # single full pattern → list[float]
)

# Mode 3: multiple observables in one pass (no state clone)
times, values = solve(
    state=state, lindblad=lindblad, t_span=(0.0, 5.0),
    save_at=[1.0, 2.0, 3.0, 5.0],
    observable=["trace:ZIII", "trace:IZII"],  # multi-pattern → list[list[float]]
)
# values[i] = [trace_ZIII_at_save_i, trace_IZII_at_save_i]
```

> **Observable format**: patterns in `"trace:<pat>"` are passed verbatim to `p.trace(pat)`. Use full Pauli strings (`"ZIII"`, `"Z?*"`) — compact qubit-index notation (`"Z0"`) is NOT supported here.

---

## Known Pitfalls

1. **Type mismatch silent at construction, noisy at solve.** Must emit `TypeError("hamiltonian and state must use the same native type (same qubit-count N)")`.
2. **`save_at` validation.** Python `solve()` must validate: non-empty, sorted, within `t_span`. Otherwise Rust silently produces wrong output.
3. **`_interface` coupling.** If `paulisum.py` renames `_interface`, bindings break. Cover with tests.
4. **Compile time.** 16 monomorphized instantiations of the full DOPRI5 + Lindblad RHS per function (×2 functions). Expected significant increase to `ppvm-python-native` build time.
5. **Returned PauliSum stale metadata.** `initial_terms`, `n_qubits`, `coefficients` on returned snapshots reflect the *input* state, not the evolved state. Document clearly; `._interface` (and `.terms`, `.trace()`) is always correct.
6. **Rates type ambiguity.** `list[list[float]]` and `list[float]` look similar. Rust side must use `extract::<Vec<f64>>()` before `extract::<Vec<Vec<f64>>>()` (more specific first).

---

## Files to Create / Modify

| File | Action |
|---|---|
| `crates/ppvm-python-native/Cargo.toml` | **Modify** — add `ppvm-timeevolve` dependency |
| `crates/ppvm-python-native/src/interface_timeevolve.rs` | **Create** — `solve_timeevolve_states` + `solve_timeevolve_observables` |
| `crates/ppvm-python-native/src/lib.rs` | **Modify** — add `pub mod interface_timeevolve;` + exports |
| `ppvm-python/src/ppvm/timeevolve.py` | **Create** — `LadderOp`, `LindbladOp`, `SolverConfig`, `solve()` |
| `ppvm-python/src/ppvm/__init__.py` | **Modify** — add `from . import timeevolve as timeevolve` |
| `ppvm-python/test/test_timeevolve.py` | **Create** — pytest tests |
| `crates/ppvm-timeevolve/plans/PLAN_PYTHON_BINDINGS.md` | **Create** — copy of this plan in project |

### Key reuse

- `crates/ppvm-python-native/src/interface.rs` — macro pattern for 16-arm monomorphization (dispatch arms mirror the existing `create_interface!` type aliases)
- `ppvm-python/src/ppvm/paulisum.py:282-291` — `__copy__` pattern for constructing returned `PauliSum` shells from native objects
- `ppvm-python/src/ppvm/paulisum.py:167-171` — `N_interface` calculation (for verifying N matches in Python)
- `crates/ppvm-timeevolve/src/solve.rs` — `solve()` function and `SolverConfig` struct to wrap
- `crates/ppvm-timeevolve/src/lindblad.rs` — `LindbladOp::new`, `JumpOp`, `LadderOp`, `LadderDirection`, `RateMatrix`

---

## Task Breakdown

### Task 31 — Write plan file + add ppvm-timeevolve dependency
**Files:** `crates/ppvm-python-native/Cargo.toml`

1. Write `crates/ppvm-timeevolve/plans/PLAN_PYTHON_BINDINGS.md` — copy of this plan
2. Add under `[dependencies]` in `crates/ppvm-python-native/Cargo.toml`:
```toml
ppvm-timeevolve = { version = "0.1.0", path = "../ppvm-timeevolve" }
```

**Verify:** `cargo build -p ppvm-python-native` succeeds with no errors.

---

### Task 32 — Create `interface_timeevolve.rs`: state-snapshot solver
**Files:** `crates/ppvm-python-native/src/interface_timeevolve.rs` (new)

Implement `solve_timeevolve_states`:

```rust
#[pyfunction]
pub fn solve_timeevolve_states(
    py: Python<'_>,
    state: &Bound<PyAny>,
    lindblad_ops: Vec<(usize, String)>,  // (qubit, "raise"|"lower")
    rates: &Bound<PyAny>,               // list[float] or list[list[float]]
    t_span_start: f64,
    t_span_end: f64,
    save_at: Vec<f64>,
    hamiltonian: Option<&Bound<PyAny>>,
    rtol: f64,
    atol: f64,
    h0: Option<f64>,
    hmin: f64,
    hmax: f64,
) -> PyResult<(Vec<f64>, Vec<PyObject>)>
```

Structure:
1. Parse `rates`: try `extract::<Vec<f64>>()` first (diagonal), then `extract::<Vec<Vec<f64>>>()` (dense), else `PyTypeError`
2. 16-arm `if let Ok(s) = state.downcast::<PauliSumIndexMapFxHash{N}>()` match
3. Inside each arm: validate `hamiltonian` downcast to same type (if Some), build `LindbladOp<T>` from ops+rates, run `ppvm_timeevolve::solve(ham, &lindblad, s, (t0, t1), &save_at, |_, p| p.clone(), config)`, wrap each result in `PauliSumIndexMapFxHash{N}` and collect as `Vec<PyObject>`
4. Use a `macro_rules!` helper to generate the 16 arms without repetition

**Verify:** `cargo build -p ppvm-python-native` succeeds. (Python tests in Task 37.)

---

### Task 33 — Add multi-observable solver to `interface_timeevolve.rs`
**Files:** `crates/ppvm-python-native/src/interface_timeevolve.rs`

Add `solve_timeevolve_observables`:

```rust
#[pyfunction]
pub fn solve_timeevolve_observables(
    state: &Bound<PyAny>,
    lindblad_ops: Vec<(usize, String)>,
    rates: &Bound<PyAny>,
    t_span_start: f64,
    t_span_end: f64,
    save_at: Vec<f64>,
    patterns: Vec<String>,            // one or more trace patterns
    hamiltonian: Option<&Bound<PyAny>>,
    rtol: f64,
    atol: f64,
    h0: Option<f64>,
    hmin: f64,
    hmax: f64,
) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)>
// outer Vec: save points; inner Vec: one f64 per pattern
```

Same dispatch macro as Task 32. Callback: `|_, p| patterns.iter().map(|pat| p.trace(pat)).collect::<Vec<f64>>()`. No state cloning.

Python `solve()` transparently handles the `list[str]` vs `str` distinction:
- `observable="trace:Z0"` → `patterns=["Z0"]`, Python unwraps inner vec → `list[float]`
- `observable=["trace:Z0","trace:Z1"]` → `patterns=["Z0","Z1"]` → `list[list[float]]`

**Verify:** `cargo build -p ppvm-python-native` succeeds.

---

### Task 34 — Register both functions in `lib.rs`
**Files:** `crates/ppvm-python-native/src/lib.rs`

Add:
```rust
pub mod interface_timeevolve;
```

Add two `#[pymodule_export]` entries in the `ppvm_python_native` module:
```rust
#[pymodule_export]
pub use crate::interface_timeevolve::solve_timeevolve_states;
#[pymodule_export]
pub use crate::interface_timeevolve::solve_timeevolve_observables;
```

**Verify:** `cargo build -p ppvm-python-native` succeeds. Then run `maturin develop --uv` in `crates/ppvm-python-native/`. Confirm in Python: `import ppvm_python_native; ppvm_python_native.solve_timeevolve_states`.

---

### Task 35 — Create `ppvm-python/src/ppvm/timeevolve.py`
**Files:** `ppvm-python/src/ppvm/timeevolve.py` (new)

Contents:
```python
from dataclasses import dataclass, field
from typing import Sequence
import ppvm_python_native
from .paulisum import PauliSum

@dataclass
class LadderOp:
    qubit: int
    direction: str  # "raise" or "lower"

@dataclass
class LindbladOp:
    jump_ops: list[LadderOp]
    rates: list[float] | list[list[float]]

@dataclass
class SolverConfig:
    rtol: float = 1e-6
    atol: float = 1e-9
    h0: float | None = None
    hmin: float = 1e-12
    hmax: float = float("inf")

def solve(
    state: PauliSum,
    lindblad: LindbladOp,
    t_span: tuple[float, float],
    save_at: Sequence[float],
    *,
    hamiltonian: PauliSum | None = None,
    observable: str | list[str] | None = None,
    config: SolverConfig | None = None,
) -> tuple[list[float], list]:
    ...
```

Responsibilities of `solve()`:
1. Validate: `save_at` non-empty, sorted ascending, within `t_span`; `t_span[0] < t_span[1]`; `direction` values are "raise" or "lower"
2. Extract `state._interface` and (if provided) `hamiltonian._interface`
3. Convert `lindblad.jump_ops` to `list[tuple[int, str]]`
4. Build `config` kwargs (use defaults from `SolverConfig` if `config=None`)
5. Parse `observable`: if `str`, wrap in list and note it's single; dispatch to `ppvm_python_native.solve_timeevolve_states` (None) or `solve_timeevolve_observables` (str/list[str]), stripping the `"trace:"` prefix from each pattern before passing to Rust
6. For state-snapshot mode: wrap each native result in a `PauliSum` shell using `__copy__`-style construction

**Verify:** `python -c "from ppvm.timeevolve import solve, LadderOp, LindbladOp"` succeeds.

---

### Task 36 — Update `ppvm-python/src/ppvm/__init__.py`
**Files:** `ppvm-python/src/ppvm/__init__.py`

Add:
```python
from . import timeevolve as timeevolve
```

**Verify:** `python -c "from ppvm import timeevolve"` succeeds.

---

### Task 37 — Write Python tests
**Files:** `ppvm-python/test/test_timeevolve.py` (new)

Required test cases:

1. **`test_decay_state_snapshots`** — pure Lindblad (no Hamiltonian), one lowering op on a 2-qubit system, check that `state.trace("ZI")` decays monotonically across save points (state-snapshot mode)
2. **`test_decay_scalar_observable`** — same setup but using `observable="trace:ZI"`, verify `list[float]` values match those from test 1 within tolerance (cross-mode consistency)
3. **`test_no_hamiltonian`** — verify `hamiltonian=None` runs without error
4. **`test_with_hamiltonian`** — verify Hamiltonian-driven oscillation in a simple 1-qubit case (trace of Z oscillates in Rabi-like setup)
5. **`test_save_at_validation_empty`** — `solve(..., save_at=[])` raises `ValueError`
6. **`test_save_at_validation_unsorted`** — `solve(..., save_at=[2.0, 1.0])` raises `ValueError`
7. **`test_save_at_validation_out_of_bounds`** — save time outside `t_span` raises `ValueError`
8. **`test_type_mismatch`** — different-N hamiltonian and state raises `TypeError`
9. **`test_returned_states_are_pauslisum`** — verify returned objects are `PauliSum` instances with working `.trace()`, `.terms`
10. **`test_dense_rate_matrix`** — use `rates=[[g, 0.0], [0.0, g]]` (dense format), verify same result as diagonal `rates=[g, g]`

**Verify:** `cd ppvm-python && uv run pytest test/test_timeevolve.py -v` — all 10 tests pass.

---

## Phase 1 / Phase 2 Split

**Phase 1 (this plan):**
- `LadderOp`, `LindbladOp` (ladder ops only), `SolverConfig`, `solve()`
- State-snapshot mode + single/multi `trace:<pattern>` scalar mode (single Rust pass)
- Non-loss `PauliSum` types (N=0..15)
- 10 Python tests

**Phase 2 (future):**
- `CollapseOp` generic dissipators
- `SolverCache` exposure for parameter sweeps
- Additional observable types (`overlap`, `n_terms`)
- LindbladOp Rust-object memoization

---

## Verification Sequence

1. `cargo build -p ppvm-python-native` — Rust builds clean
2. `cargo clippy -p ppvm-python-native -- -D warnings` — no warnings
3. `cd crates/ppvm-python-native && maturin develop --uv` — wheel installs
4. `cd ppvm-python && uv run pytest test/ -v` — all tests (including pre-existing) pass
