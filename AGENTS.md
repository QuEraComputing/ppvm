# AGENTS.md

This file provides guidance to AI agents when working with code in this repository.

## Build and Test

```bash
# Rust
cargo test --workspace                       # Run all Rust tests
cargo test -p ppvm-tableau                   # Test a single crate
cargo test -p ppvm-runtime -- test_ghz       # Run a single test by name
cargo bench -p ppvm-tableau --bench micro    # Run benchmarks for a crate
cargo bench --bench micro -- "gates/single-qubit/h"  # Run a specific benchmark

# Python (requires uv; compiles ppvm-python-native via maturin automatically)
uv run --project ppvm-python --group dev pytest ppvm-python/test/
uv run --project ppvm-python --group dev pytest ppvm-python/test/test_basics.py  # Single file
uv run --project ppvm-python --group dev pytest ppvm-python/test/ -k test_ghz   # Single test
```

Rust edition 2024. CI sets `RUSTFLAGS="-C target-feature=+aes,+sse2"` (needed for gxhash on x86).

## Python Bindings

**Two-layer build:** `ppvm-python-native` (Rust → cdylib via maturin + PyO3 0.27) is compiled automatically when `ppvm-python` is installed. `ppvm-python` is a pure Python wrapper using `uv_build` as its build backend.

- Python ≥ 3.10 required (`.python-version` pins 3.12 for dev)
- `uv` manages the venv, deps, and triggers the maturin build
- `ppvm-python/pyproject.toml` references `ppvm-python-native` via `[tool.uv.sources]` path dependency
- The native module exports 16 PauliSum variants × 2 (with/without loss) + 32 GeneralizedTableau variants (1–32 qubits) via `create_interface!` / `create_interface_range!` macros

## Commit Messages

Use [Conventional Commits](https://www.conventionalcommits.org/): `<type>(<scope>): <description>`

Examples: `feat(tableau): add correlated loss channel`, `fix(runtime): handle zero-norm in truncation`

## Project Structure

```
ppvm (root package)
├── crates/
│   ├── ppvm-runtime     # Core: Pauli arithmetic, PauliSum, gate traits, configs
│   ├── ppvm-tableau     # Stabilizer tableau simulator (Clifford + non-Clifford)
│   ├── ppvm-sym         # Symbolic (parametric) Pauli propagation
│   └── ppvm-python-native  # PyO3 bindings
├── ppvm-python/         # Pure Python wrapper (separate build, uses uv)
└── examples/            # Rust examples (symbolic.rs, trotter.rs)
```

**Dependency graph:** `ppvm-runtime` is the foundation. `ppvm-tableau` and `ppvm-sym` depend on it. `ppvm-python-native` depends on both `ppvm-runtime` and `ppvm-tableau`.

## Architecture

ppvm implements two quantum simulation backends:

### 1. Pauli Propagation (`ppvm-runtime`)

Tracks Pauli operator evolution through circuits in the **Heisenberg picture** (circuits run backwards). The central type is `PauliSum<T: Config>`, a dictionary of Pauli strings to coefficients.

Key design patterns:
- **Config-based generics:** `Config` trait bundles Storage, Coefficient, Strategy, Map, and BuildHasher choices at compile time. Implementations live in `config/` (fxhash, indexmap, dashmap, gxhash).
- **Dual-map optimization:** `PauliSum` maintains two internal maps (main + auxiliary) and swaps between them during gate propagation to avoid repeated allocations.
- **Strategy pattern:** Truncation policies (`CoefficientThreshold`, `MaxPauliWeight`, `MaxLossWeight`, `CombinedStrategy`) control when small terms are dropped. Call `.truncate()` to apply.
- **Backward propagation:** Pauli propagation runs circuits backwards. To simulate `H(0); CNOT(0,1)`, call `state.cnot(0,1); state.h(0)` — the CNOT precedes the Hadamard.

### 2. Generalized Stabilizer Tableau (`ppvm-tableau`)

Full state simulation using stabilizer formalism, extended to handle non-Clifford gates (T, rotations) via stabilizer rank decomposition with sparse coefficient tracking.

- **`Tableau<T: Config>`:** 2n-row stabilizer/destabilizer tableau (rows 0..n = destabilizers, n..2n = stabilizers).
- **`GeneralizedTableau<T: Config, IndexType>`:** Extends Tableau with a sparse coefficient vector for non-Clifford state tracking. `IndexType` can be `usize`, `u128`, or `bnum::types::U256` for large qubit counts.
- **`SparseVector<T, I>` trait:** Stores coefficients indexed by bitstrings. Indices can be large integers (U256, U512, U1024) for simulations beyond 64 qubits.
- **Stim compatibility:** Parse Stim circuits with `ppvm-stim` (`StimProgram.parse` / `StimProgram.from_file`); execute with `tab.run(prog)` or sample many shots with `ppvm.sample_stim` / `GeneralizedTableau.sample`.

### Trait hierarchy (in `ppvm-runtime/src/traits/`)

Gate behavior is defined via traits reused across both backends:
- `Clifford` / `CliffordExtensions` — single/two-qubit Clifford gates
- `TGate`, `RotationOne`, `RotationTwo`, `U3Gate` — non-Clifford gates (branching)
- `Measure` / `LossyMeasure` — Z-basis measurement
- `Depolarizing`, `PauliError`, `LossChannel`, `CorrelatedLossChannel` — noise channels

### Python bindings

`ppvm-python-native` uses PyO3 macros (`create_interface!`, `create_interface_range!`) to generate multiple Python classes per config/qubit-count combination. `ppvm-python` wraps these with Pythonic APIs via mixins.

## `crates/ppvm-stim` test corpus

Tests under `crates/ppvm-stim/tests/data/` are committed `.stim` + `.expected.json` pairs consumed by `tests/stim_corpus.rs`. The harness asserts ppvm's output matches the committed reference bit-for-bit. Cross-check against `quantumlib/Stim` happens at regen time, not at test time:

```bash
cd crates/ppvm-stim/tests/regen-stim
uv sync
uv run regen-stim all
```

When phase-2 lifts a restriction, `uv run regen-stim refresh ../data/unsupported/<name>.stim` flips that fixture from "expected to fail normalize" to "expected to match Stim's pre-recorded distribution". Schema and category overview: `crates/ppvm-stim/tests/data/README.md`; design rationale: `docs/superpowers/specs/2026-04-28-ppvm-stim-test-corpus-design.md`.
