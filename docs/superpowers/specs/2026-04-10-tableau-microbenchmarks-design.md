# Tableau Microbenchmark Suite

## Goal

Add a focused microbenchmark suite for the generalized tableau implementation (`ppvm-tableau`) to track performance improvements and regressions. The benchmarks should provide statistically meaningful measurements while keeping total wall-clock time under 5 minutes.

## Changes

### 1. Tune Criterion settings on all existing benchmark files

Apply to `benches/tableau.rs`, `benches/tableau-msd.rs`, `benches/tableau-msd-stim.rs`:

- Warmup: 3s -> 1s
- Measurement time: 5s -> 3s
- Sample size: 100 -> 50

This preserves measurement accuracy (all per-iteration times are sub-300us) while cutting total suite time roughly in half.

### 2. New file: `benches/micro.rs`

Single benchmark binary with 6 criterion benchmark groups. Tableau benchmarks (Groups 1-5) use `Byte8F64<2>` config with `usize` index type on a 32-qubit tableau. Sparse vector benchmarks (Group 6) test across multiple index types.

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

#### Group 4: `measurement`

Three distinct measurement code paths:

- `deterministic` -- measure qubit in |0> (no anticommuting stabilizer, fast path)
- `random` -- measure qubit in |+> (anticommuting stabilizer found, tableau update required)
- `generalized` -- measure on GeneralizedTableau after 4 T gates (coefficient-aware path with overlap computation)

Note: `compute_decomposition` and `compute_phase` are `pub(crate)` and cannot be benchmarked directly from `benches/`. They are exercised indirectly via the `generalized` measurement benchmark above.

#### Group 5: `noise`

Noise channels on a fresh 32-qubit `GeneralizedTableau`. Probabilities set to 1.0 (or summing to 1.0) to ensure the noise always fires, avoiding RNG-dependent early returns that would skew measurements.

Benchmarks:
- `depolarize` -- depolarize(0, 1.0)
- `pauli_error` -- pauli_error(0, [1/3, 1/3, 1/3])
- `two_qubit_pauli_error` -- two_qubit_pauli_error(0, 1, [1/15; 15])
- `depolarize2` -- depolarize2(0, 1, 1.0)
- `loss_channel` -- loss_channel(0, 1.0)
- `correlated_loss_channel` -- correlated_loss_channel(0, 1, [0.5, 0.3, 0.2])

#### Group 6: `sparse-vec`

Benchmarks for `SparseVector` trait operations on `Vec<(Complex64, I)>`, focusing on how index type size affects performance. Each operation is benchmarked with 16 pre-populated entries (realistic post-4-T-gate coefficient count).

**Index types** (3 levels):
- `usize` (64-bit) -- baseline, up to 64 qubits
- `u128` (128-bit) -- medium, up to 128 qubits
- `U256` (`bnum::types::U256`, 256-bit) -- large, common big-int type from the Python interface

**Operations** (8 per index type):
- `unsafe_insert` -- append a new entry (O(1))
- `add_or_insert/existing` -- update an existing entry by index scan (O(n))
- `add_or_insert/new` -- insert a new entry after full scan (O(n))
- `get` -- look up an entry by index (O(n))
- `mul_by` -- scale all entries (O(n))
- `mul_element_by` -- find and scale one entry (O(n))
- `trim` -- filter entries below threshold (O(n))
- `normalize` -- compute norm and rescale (O(n), two passes)

**Setup**: Each benchmark clones a pre-built vector with 16 entries using `iter_batched_ref`. Indices are spread across the index space (not sequential) to reflect realistic tableau usage where indices are bitstrings with gaps.

### 3. Cargo.toml addition

Add a `[[bench]]` entry for the new benchmark:

```toml
[[bench]]
name = "micro"
harness = false
```

## Benchmark count estimate

- Group 1 (gates/single-qubit): 10 benchmarks
- Group 2 (gates/two-qubit): 3 benchmarks
- Group 3 (gates/non-clifford): 5 benchmarks
- Group 4 (measurement): 3 benchmarks
- Group 5 (noise): 6 benchmarks
- Group 6 (sparse-vec): 24 benchmarks (8 ops x 3 index types)
- **Total: 51 benchmarks**

At ~4s per benchmark (1s warmup + 3s measurement), the new file takes ~3.5 minutes. Combined with the tuned existing benchmarks (~1 minute), the full suite runs in ~4.5 minutes.

## Usage

```bash
# Run only microbenchmarks
cargo bench --bench micro

# Run a specific group
cargo bench --bench micro -- "gates/single-qubit"

# Run only sparse vector benchmarks
cargo bench --bench micro -- "sparse-vec"

# Run sparse vector benchmarks for a specific index type
cargo bench --bench micro -- "sparse-vec/U256"

# Run a specific benchmark
cargo bench --bench micro -- "gates/single-qubit/h"

# Run all benchmarks
cargo bench --package ppvm-tableau
```
