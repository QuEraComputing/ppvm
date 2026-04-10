# Tableau Microbenchmark Suite

## Goal

Add a focused microbenchmark suite for the generalized tableau implementation (`ppvm-tableau`) to track performance improvements and regressions. The benchmarks should execute quickly (target ~1 minute total wall-clock time) while providing statistically meaningful measurements.

## Changes

### 1. Tune Criterion settings on all existing benchmark files

Apply to `benches/tableau.rs`, `benches/tableau-msd.rs`, `benches/tableau-msd-stim.rs`:

- Warmup: 3s -> 1s
- Measurement time: 5s -> 3s
- Sample size: 100 -> 50

This preserves measurement accuracy (all per-iteration times are sub-300us) while cutting total suite time roughly in half.

### 2. New file: `benches/micro.rs`

Single benchmark binary with 6 criterion benchmark groups. All benchmarks use `Byte8F64<2>` config with `usize` index type on a 32-qubit tableau.

Same Criterion settings as above (1s warmup, 3s measurement, 50 samples).

#### Group 1: `gates/single-qubit`

Each gate applied once to qubit 0 on a fresh 32-qubit `GeneralizedTableau`, using `iter_batched_ref` with `fork` for setup isolation.

Benchmarks: `h`, `s`, `s_adj`, `x`, `y`, `z`, `sqrt_x`, `sqrt_x_adj`, `sqrt_y`, `sqrt_y_adj`

#### Group 2: `gates/two-qubit`

Each gate applied to qubits (0, 1) on a fresh tableau.

Benchmarks: `cnot`, `cz`, `cy`

#### Group 3: `gates/non-clifford`

Each gate applied to a tableau where qubit 0 is in |+> state (via prior H gate) to trigger the branching code path. For `rxx`, both qubits 0 and 1 are in |+>.

Benchmarks:
- `t` -- T gate on qubit 0
- `t_adj` -- T-adjoint on qubit 0
- `rx` -- rx(pi/4) on qubit 0 (calls `rotate_1(Pauli::X, 0, pi/4)`)
- `rxx` -- rxx(pi/4) on qubits (0, 1) (calls `rotate_2` with XX axis)
- `u3` -- u3(pi/4, pi/4, pi/4) on qubit 0

#### Group 4: `internals`

Sparse vector operations benchmarked on a coefficient vector with ~16 entries (from a 32-qubit tableau after 4 T gates on qubits in |+>).

Note: `compute_decomposition` and `compute_phase` are `pub(crate)` and cannot be benchmarked directly from `benches/`. They are exercised indirectly via the generalized measurement benchmark in Group 5.

Benchmarks:
- `sparse_vec_normalize` -- normalize a cloned coefficient vector
- `sparse_vec_trim` -- trim coefficients below threshold on a cloned vector

#### Group 5: `measurement`

Three distinct measurement code paths:

- `deterministic` -- measure qubit in |0> (no anticommuting stabilizer, fast path)
- `random` -- measure qubit in |+> (anticommuting stabilizer found, tableau update required)
- `generalized` -- measure on GeneralizedTableau after 4 T gates (coefficient-aware path with overlap computation)

#### Group 6: `noise`

Noise channels on a fresh 32-qubit `GeneralizedTableau`. Probabilities set to 1.0 (or summing to 1.0) to ensure the noise always fires, avoiding RNG-dependent early returns that would skew measurements.

Benchmarks:
- `depolarize` -- depolarize(0, 1.0)
- `pauli_error` -- pauli_error(0, [1/3, 1/3, 1/3])
- `two_qubit_pauli_error` -- two_qubit_pauli_error(0, 1, [1/15; 15])
- `depolarize2` -- depolarize2(0, 1, 1.0)
- `loss_channel` -- loss_channel(0, 1.0)
- `correlated_loss_channel` -- correlated_loss_channel(0, 1, [0.5, 0.3, 0.2])

### 3. Cargo.toml addition

Add a `[[bench]]` entry for the new benchmark:

```toml
[[bench]]
name = "micro"
harness = false
```

## Benchmark count estimate

- Group 1: 10 benchmarks
- Group 2: 3 benchmarks
- Group 3: 5 benchmarks
- Group 4: 2 benchmarks
- Group 5: 3 benchmarks
- Group 6: 6 benchmarks
- **Total: 29 benchmarks**

At ~4s per benchmark (1s warmup + 3s measurement), the new file takes ~2 minutes. Combined with the tuned existing benchmarks (~1 minute), the full suite runs in ~3 minutes.

## Usage

```bash
# Run only microbenchmarks
cargo bench --bench micro

# Run a specific group
cargo bench --bench micro -- "gates/single-qubit"

# Run a specific benchmark
cargo bench --bench micro -- "gates/single-qubit/h"

# Run all benchmarks
cargo bench --package ppvm-tableau
```
