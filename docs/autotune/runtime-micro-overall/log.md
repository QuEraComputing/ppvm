# runtime-micro-overall

## Target

Improve the overall `ppvm-runtime` `micro` benchmark suite, using the first fresh full-suite run in this session as the baseline.

## Architecture Notes

- The slowest representative microbenchmarks in the baseline run clustered around the shared `PauliSum` transform paths: `clifford/single/h` at `1.053 µs`, `clifford/two-qubit/cnot` at `1.050 µs`, and `scaling/cnot/1000` at `1.418 µs`.
- A quick ad-hoc profile showed cloning was not the dominant cost for representative operations on the benchmark state shape (`clone_only_ns≈72-83ns` versus total operation timings in the `170-295ns` range in the profiler), so the first optimization hypothesis focused on transform mechanics rather than clone elimination.
- Single-word profiling did not show `rehash()` dominating gate cost, so rehash-specific work was deprioritized.

## Iterations

### owned-bijective-clifford-transform

- Status: discard
- Hypothesis: Clifford conjugation is bijective, so consuming owned map entries should avoid clone-heavy `map_add` work and speed up `h`/`cnot`/related scaling benchmarks.
- Result: strong regression.
- Evidence:
  - `clifford/single/h`: `1.053 µs -> 1.362 µs`
  - `clifford/two-qubit/cnot`: `1.050 µs -> 1.441 µs`
  - `scaling/cnot/1000`: `1.418 µs -> 2.533 µs`
- Takeaway: the current `map_add` path plus existing map behavior is materially better than the owned-entry transform on this workload.

### noise-factor-specialization

- Status: keep
- Hypothesis: the noise group spends avoidable time in repeated branchy factor selection and, for `depolarize2`, unnecessary construction of a 15-entry probability array followed by the generic two-qubit channel.
- Result: keep.
- Implemented:
  - Simplified `pauli_error` and `depolarize` to use bit-based Pauli classification and precomputed single-qubit factors.
  - Replaced `depolarize2` with the closed-form two-qubit depolarizing factor `1 - 16p/15` for any non-identity two-qubit Pauli term.
  - Reverted the attempted generic `two_qubit_pauli_error` factor table after it added too much fixed overhead.
- Evidence from the final full-suite run:
  - `noise/pauli_error`: `279.81 ns -> 213.57 ns`
  - `noise/depolarize`: `272.73 ns -> 217.39 ns`
  - `noise/depolarize2`: `509.21 ns -> 269.71 ns`
  - `noise/amplitude_damping`: `524.70 ns -> 458.56 ns`
  - `noise/two_qubit_pauli_error`: `523.07 ns -> 457.41 ns`

## Current Best

- Keep the noise specialization in `crates/ppvm-runtime/src/sum/noise.rs`.
- Discard the owned-entry Clifford transform idea unless new profiling reveals a different map implementation bottleneck.
