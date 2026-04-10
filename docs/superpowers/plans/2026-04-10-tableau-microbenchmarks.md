# Tableau Microbenchmark Suite Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a microbenchmark suite (`benches/micro.rs`) for `ppvm-tableau` and tune existing benchmark Criterion settings for faster CI runs.

**Architecture:** Single new benchmark binary with 6 criterion groups covering gates, measurement, noise, and sparse vector operations across index types. Existing benchmark files get Criterion config tuned (1s warmup, 3s measurement, 50 samples).

**Tech Stack:** Rust, Criterion 0.7.0, bnum (U256), ppvm-tableau, ppvm-runtime

---

### Task 1: Tune Criterion settings on existing benchmark files

**Files:**
- Modify: `crates/ppvm-tableau/benches/tableau.rs`
- Modify: `crates/ppvm-tableau/benches/tableau-msd.rs`
- Modify: `crates/ppvm-tableau/benches/tableau-msd-stim.rs`

- [ ] **Step 1: Add Criterion config to `benches/tableau.rs`**

Add `use std::time::Duration;` at the top of the file and replace the `criterion_group!` macro at the bottom:

```rust
// old:
criterion_group!(benches, tableau_scaling_benchmarks);

// new:
criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50);
    targets = tableau_scaling_benchmarks
}
```

- [ ] **Step 2: Add Criterion config to `benches/tableau-msd.rs`**

Add `use std::time::Duration;` at the top and replace the `criterion_group!` macro:

```rust
// old:
criterion_group!(benches, msd_benchmarks);

// new:
criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50);
    targets = msd_benchmarks
}
```

- [ ] **Step 3: Add Criterion config to `benches/tableau-msd-stim.rs`**

Add `use std::time::Duration;` at the top and replace the `criterion_group!` macro:

```rust
// old:
criterion_group!(benches, msd_stim_benchmarks);

// new:
criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50);
    targets = msd_stim_benchmarks
}
```

- [ ] **Step 4: Verify existing benchmarks still compile and run**

Run: `cargo bench --package ppvm-tableau --bench tableau -- --test`
Run: `cargo bench --package ppvm-tableau --bench tableau-msd -- --test`
Run: `cargo bench --package ppvm-tableau --bench tableau-msd-stim -- --test`

Expected: All three compile and run a quick test pass (the `--test` flag runs each benchmark exactly once to verify it works without doing full measurement).

- [ ] **Step 5: Commit**

```bash
git add crates/ppvm-tableau/benches/tableau.rs crates/ppvm-tableau/benches/tableau-msd.rs crates/ppvm-tableau/benches/tableau-msd-stim.rs
git commit -m "Tune Criterion settings on existing benchmarks (1s warmup, 3s measurement, 50 samples)"
```

---

### Task 2: Add `[[bench]]` entry and scaffold `benches/micro.rs`

**Files:**
- Modify: `crates/ppvm-tableau/Cargo.toml`
- Create: `crates/ppvm-tableau/benches/micro.rs`

- [ ] **Step 1: Add bench entry to Cargo.toml**

Append to `crates/ppvm-tableau/Cargo.toml`:

```toml
[[bench]]
name = "micro"
harness = false
```

- [ ] **Step 2: Create minimal `benches/micro.rs` with Criterion config and empty group**

```rust
use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<Byte8F64<2>, usize>;

const N_QUBITS: usize = 32;

fn configure() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50)
}

fn bench_single_qubit_gates(c: &mut Criterion) {
    let tab = Tab::new(N_QUBITS, 1e-10);
    let mut group = c.benchmark_group("gates/single-qubit");

    group.bench_function("h", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.h(0), criterion::BatchSize::SmallInput);
    });

    group.finish();
}

criterion_group! {
    name = benches;
    config = configure();
    targets = bench_single_qubit_gates
}
criterion_main!(benches);
```

- [ ] **Step 3: Verify scaffold compiles and runs**

Run: `cargo bench --package ppvm-tableau --bench micro -- --test`

Expected: Compiles and runs a single quick test pass for the `h` benchmark.

- [ ] **Step 4: Commit**

```bash
git add crates/ppvm-tableau/Cargo.toml crates/ppvm-tableau/benches/micro.rs
git commit -m "Scaffold micro benchmark with Criterion config and single test gate"
```

---

### Task 3: Implement Group 1 — single-qubit gate benchmarks

**Files:**
- Modify: `crates/ppvm-tableau/benches/micro.rs`

- [ ] **Step 1: Complete `bench_single_qubit_gates` with all 10 gates**

Replace the `bench_single_qubit_gates` function:

```rust
fn bench_single_qubit_gates(c: &mut Criterion) {
    let tab = Tab::new(N_QUBITS, 1e-10);
    let mut group = c.benchmark_group("gates/single-qubit");

    group.bench_function("h", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.h(0), criterion::BatchSize::SmallInput);
    });
    group.bench_function("s", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.s(0), criterion::BatchSize::SmallInput);
    });
    group.bench_function("s_adj", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.s_adj(0), criterion::BatchSize::SmallInput);
    });
    group.bench_function("x", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.x(0), criterion::BatchSize::SmallInput);
    });
    group.bench_function("y", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.y(0), criterion::BatchSize::SmallInput);
    });
    group.bench_function("z", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.z(0), criterion::BatchSize::SmallInput);
    });
    group.bench_function("sqrt_x", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.sqrt_x(0), criterion::BatchSize::SmallInput);
    });
    group.bench_function("sqrt_x_adj", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.sqrt_x_adj(0), criterion::BatchSize::SmallInput);
    });
    group.bench_function("sqrt_y", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.sqrt_y(0), criterion::BatchSize::SmallInput);
    });
    group.bench_function("sqrt_y_adj", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.sqrt_y_adj(0), criterion::BatchSize::SmallInput);
    });

    group.finish();
}
```

- [ ] **Step 2: Verify it compiles and runs**

Run: `cargo bench --package ppvm-tableau --bench micro -- "gates/single-qubit" --test`

Expected: All 10 single-qubit gate benchmarks pass the test run.

- [ ] **Step 3: Commit**

```bash
git add crates/ppvm-tableau/benches/micro.rs
git commit -m "Add single-qubit gate microbenchmarks (h, s, s_adj, x, y, z, sqrt_x/y and adjoints)"
```

---

### Task 4: Implement Group 2 — two-qubit gate benchmarks

**Files:**
- Modify: `crates/ppvm-tableau/benches/micro.rs`

- [ ] **Step 1: Add `bench_two_qubit_gates` function**

Add this function before the `criterion_group!` macro:

```rust
fn bench_two_qubit_gates(c: &mut Criterion) {
    let tab = Tab::new(N_QUBITS, 1e-10);
    let mut group = c.benchmark_group("gates/two-qubit");

    group.bench_function("cnot", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.cnot(0, 1), criterion::BatchSize::SmallInput);
    });
    group.bench_function("cz", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.cz(0, 1), criterion::BatchSize::SmallInput);
    });
    group.bench_function("cy", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.cy(0, 1), criterion::BatchSize::SmallInput);
    });

    group.finish();
}
```

- [ ] **Step 2: Register in `criterion_group!`**

Update the targets:

```rust
criterion_group! {
    name = benches;
    config = configure();
    targets = bench_single_qubit_gates, bench_two_qubit_gates
}
```

- [ ] **Step 3: Verify it compiles and runs**

Run: `cargo bench --package ppvm-tableau --bench micro -- "gates/two-qubit" --test`

Expected: All 3 two-qubit gate benchmarks pass.

- [ ] **Step 4: Commit**

```bash
git add crates/ppvm-tableau/benches/micro.rs
git commit -m "Add two-qubit gate microbenchmarks (cnot, cz, cy)"
```

---

### Task 5: Implement Group 3 — non-Clifford gate benchmarks

**Files:**
- Modify: `crates/ppvm-tableau/benches/micro.rs`

- [ ] **Step 1: Add `bench_non_clifford_gates` function**

Add this function. The setup puts qubit 0 (and qubit 1 for `rxx`) into |+> via H, which triggers the branching path for non-Clifford gates:

```rust
fn bench_non_clifford_gates(c: &mut Criterion) {
    // Qubit 0 in |+> to trigger branching
    let mut tab = Tab::new(N_QUBITS, 1e-10);
    tab.h(0);

    // Both qubits 0 and 1 in |+> for rxx
    let mut tab_2q = Tab::new(N_QUBITS, 1e-10);
    tab_2q.h(0);
    tab_2q.h(1);

    let mut group = c.benchmark_group("gates/non-clifford");

    let pi_4 = std::f64::consts::FRAC_PI_4;

    group.bench_function("t", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.t(0), criterion::BatchSize::SmallInput);
    });
    group.bench_function("t_adj", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.t_adj(0), criterion::BatchSize::SmallInput);
    });
    group.bench_function("rx", |b| {
        b.iter_batched_ref(|| tab.fork(None), |t| t.rx(0, pi_4), criterion::BatchSize::SmallInput);
    });
    group.bench_function("rxx", |b| {
        b.iter_batched_ref(
            || tab_2q.fork(None),
            |t| t.rxx(0, 1, pi_4),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("u3", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.u3(0, pi_4, pi_4, pi_4),
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}
```

- [ ] **Step 2: Register in `criterion_group!`**

```rust
criterion_group! {
    name = benches;
    config = configure();
    targets = bench_single_qubit_gates, bench_two_qubit_gates, bench_non_clifford_gates
}
```

- [ ] **Step 3: Verify it compiles and runs**

Run: `cargo bench --package ppvm-tableau --bench micro -- "gates/non-clifford" --test`

Expected: All 5 non-Clifford benchmarks pass.

- [ ] **Step 4: Commit**

```bash
git add crates/ppvm-tableau/benches/micro.rs
git commit -m "Add non-Clifford gate microbenchmarks (t, t_adj, rx, rxx, u3)"
```

---

### Task 6: Implement Group 4 — measurement benchmarks

**Files:**
- Modify: `crates/ppvm-tableau/benches/micro.rs`

- [ ] **Step 1: Add `bench_measurement` function**

Three measurement paths: deterministic (qubit in |0>), random (qubit in |+>), and generalized (after T gates with coefficients):

```rust
fn bench_measurement(c: &mut Criterion) {
    // Deterministic: qubit 0 is |0>, measurement is deterministic
    let tab_det = Tab::new(N_QUBITS, 1e-10);

    // Random: qubit 0 is |+>, measurement triggers tableau update
    let mut tab_rand = Tab::new(N_QUBITS, 1e-10);
    tab_rand.h(0);

    // Generalized: 4 T gates create ~16 coefficient branches
    let mut tab_gen = Tab::new(N_QUBITS, 1e-10);
    for i in 0..4 {
        tab_gen.h(i);
        tab_gen.t(i);
    }

    let mut group = c.benchmark_group("measurement");

    group.bench_function("deterministic", |b| {
        b.iter_batched_ref(
            || tab_det.fork(None),
            |t| t.measure(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("random", |b| {
        b.iter_batched_ref(
            || tab_rand.fork(None),
            |t| t.measure(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("generalized", |b| {
        b.iter_batched_ref(
            || tab_gen.fork(None),
            |t| t.measure(0),
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}
```

- [ ] **Step 2: Register in `criterion_group!`**

```rust
criterion_group! {
    name = benches;
    config = configure();
    targets = bench_single_qubit_gates, bench_two_qubit_gates, bench_non_clifford_gates, bench_measurement
}
```

- [ ] **Step 3: Verify it compiles and runs**

Run: `cargo bench --package ppvm-tableau --bench micro -- "measurement" --test`

Expected: All 3 measurement benchmarks pass.

- [ ] **Step 4: Commit**

```bash
git add crates/ppvm-tableau/benches/micro.rs
git commit -m "Add measurement microbenchmarks (deterministic, random, generalized)"
```

---

### Task 7: Implement Group 5 — noise benchmarks

**Files:**
- Modify: `crates/ppvm-tableau/benches/micro.rs`

- [ ] **Step 1: Add `bench_noise` function**

Probabilities set to 1.0 (or summing to 1.0) to ensure noise always fires:

```rust
fn bench_noise(c: &mut Criterion) {
    let tab = Tab::new(N_QUBITS, 1e-10);
    let mut group = c.benchmark_group("noise");

    group.bench_function("depolarize", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.depolarize(0, 1.0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("pauli_error", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.pauli_error(0, [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("two_qubit_pauli_error", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.two_qubit_pauli_error(0, 1, [1.0 / 15.0; 15]),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("depolarize2", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.depolarize2(0, 1, 1.0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("loss_channel", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.loss_channel(0, 1.0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("correlated_loss_channel", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.correlated_loss_channel(0, 1, [0.5, 0.3, 0.2]),
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}
```

- [ ] **Step 2: Register in `criterion_group!`**

```rust
criterion_group! {
    name = benches;
    config = configure();
    targets = bench_single_qubit_gates, bench_two_qubit_gates, bench_non_clifford_gates,
              bench_measurement, bench_noise
}
```

- [ ] **Step 3: Verify it compiles and runs**

Run: `cargo bench --package ppvm-tableau --bench micro -- "noise" --test`

Expected: All 6 noise benchmarks pass.

- [ ] **Step 4: Commit**

```bash
git add crates/ppvm-tableau/benches/micro.rs
git commit -m "Add noise channel microbenchmarks (depolarize, pauli_error, loss, correlated loss)"
```

---

### Task 8: Implement Group 6 — sparse vector benchmarks

**Files:**
- Modify: `crates/ppvm-tableau/benches/micro.rs`

- [ ] **Step 1: Add `bnum` import and sparse vector helper**

Add at the top of the file:

```rust
use bnum::types::U256;
use num::complex::Complex64;
```

Add a generic helper function that builds a pre-populated sparse vector with 16 entries whose indices are spread across the index space (bit-shifted to simulate realistic bitstring gaps):

```rust
fn make_sparse_vec<I: TableauIndex>(n: usize) -> Vec<(Complex64, I)> {
    let mut vec: Vec<(Complex64, I)> = SparseVector::new();
    for k in 0..n {
        // Spread indices across the space by shifting: index = k << 2
        // This avoids sequential indices and better simulates real tableau usage
        let index = I::from(k as u8) << 2;
        let value = Complex64::new(1.0 / (k as f64 + 1.0), 0.1 * k as f64);
        vec.unsafe_insert(index, value);
    }
    vec
}
```

- [ ] **Step 2: Add `bench_sparse_vec` function with generic inner function**

Add a generic inner function that benchmarks all 8 operations for a given index type, then call it for each type:

```rust
fn bench_sparse_vec_for_type<I: TableauIndex + std::fmt::Debug>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
    type_name: &str,
) {
    let vec16 = make_sparse_vec::<I>(16);
    // An index known to exist (index = 4 << 2 = 16 -> I::from(4) << 2, i.e. the 5th element)
    let existing_index = I::from(4u8) << 2;
    // An index known not to exist
    let new_index = I::from(99u8);

    group.bench_function(format!("{type_name}/unsafe_insert"), |b| {
        b.iter_batched_ref(
            || vec16.clone(),
            |v| v.unsafe_insert(new_index, Complex64::new(1.0, 0.0)),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function(format!("{type_name}/add_or_insert_existing"), |b| {
        b.iter_batched_ref(
            || vec16.clone(),
            |v| v.add_or_insert(existing_index, Complex64::new(0.5, 0.0)),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function(format!("{type_name}/add_or_insert_new"), |b| {
        b.iter_batched_ref(
            || vec16.clone(),
            |v| v.add_or_insert(new_index, Complex64::new(0.5, 0.0)),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function(format!("{type_name}/get"), |b| {
        b.iter_batched_ref(
            || vec16.clone(),
            |v| v.get(&existing_index),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function(format!("{type_name}/mul_by"), |b| {
        b.iter_batched_ref(
            || vec16.clone(),
            |v| v.mul_by(Complex64::new(2.0, 0.0)),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function(format!("{type_name}/mul_element_by"), |b| {
        b.iter_batched_ref(
            || vec16.clone(),
            |v| v.mul_element_by(existing_index, Complex64::new(2.0, 0.0)),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function(format!("{type_name}/trim"), |b| {
        b.iter_batched_ref(
            || vec16.clone(),
            |v| v.trim(Complex64::new(0.05, 0.0)),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function(format!("{type_name}/normalize"), |b| {
        b.iter_batched_ref(
            || vec16.clone(),
            |v| v.normalize(),
            criterion::BatchSize::SmallInput,
        );
    });
}

fn bench_sparse_vec(c: &mut Criterion) {
    let mut group = c.benchmark_group("sparse-vec");

    bench_sparse_vec_for_type::<usize>(&mut group, "usize");
    bench_sparse_vec_for_type::<u128>(&mut group, "u128");
    bench_sparse_vec_for_type::<U256>(&mut group, "U256");

    group.finish();
}
```

- [ ] **Step 3: Register in `criterion_group!`**

```rust
criterion_group! {
    name = benches;
    config = configure();
    targets = bench_single_qubit_gates, bench_two_qubit_gates, bench_non_clifford_gates,
              bench_measurement, bench_noise, bench_sparse_vec
}
```

- [ ] **Step 4: Verify it compiles and runs**

Run: `cargo bench --package ppvm-tableau --bench micro -- "sparse-vec" --test`

Expected: All 24 sparse vector benchmarks (8 ops x 3 types) pass the test run.

- [ ] **Step 5: Commit**

```bash
git add crates/ppvm-tableau/benches/micro.rs
git commit -m "Add sparse vector microbenchmarks across usize, u128, and U256 index types"
```

---

### Task 9: Full suite validation

**Files:** None (verification only)

- [ ] **Step 1: Run the full micro benchmark suite end-to-end**

Run: `cargo bench --package ppvm-tableau --bench micro -- --test`

Expected: All 51 benchmarks (10 + 3 + 5 + 3 + 6 + 24) pass the test run without errors.

- [ ] **Step 2: Run all benchmarks to verify total wall-clock time**

Run: `cargo bench --package ppvm-tableau 2>&1 | tail -5`

Expected: All benchmarks complete. Total time should be under 5 minutes.

- [ ] **Step 3: Spot-check a few benchmark results for sanity**

Run: `cargo bench --package ppvm-tableau --bench micro -- "gates/single-qubit/h"`

Expected: Result shows a time in the low microsecond range (< 10 µs), consistent with a single gate on 32 qubits.

- [ ] **Step 4: Commit (no-op if nothing changed)**

No commit needed unless issues were found and fixed.
