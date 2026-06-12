// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Baseline + batched-vs-loop benches for every Clifford / CliffordExtensions
//! gate on `Tableau<T>`, at qubit counts spanning the per-word storage boundary.
//! Re-run after promoting the batched methods to a trait to measure the speedup
//! on gates that newly gain a batched implementation.

use std::time::Duration;

use criterion::{BatchSize, BenchmarkId, Criterion, criterion_group, criterion_main};
use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_runtime::prelude::*;
use ppvm_tableau::prelude::*;

// Two u64 words = up to 128 qubits, with the boundary at qubit 64.
type Tab = Tableau<Byte8F64<2>>;

const SIZES: &[usize] = &[128];

fn configure() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50)
}

/// Non-trivial starting tableau so bit patterns aren't all identity.
fn setup(n: usize) -> Tab {
    let mut tab = Tab::new(n);
    for q in 0..n.min(8) {
        tab.h(q);
    }
    for q in (0..n.min(8)).step_by(2) {
        tab.s(q);
    }
    tab
}

fn indices_all(n: usize) -> Vec<usize> {
    (0..n).collect()
}

fn indices_every_other(n: usize) -> Vec<usize> {
    (0..n).step_by(2).collect()
}

/// Disjoint adjacent pairs covering all qubits: (0,1), (2,3), ...
fn pairs_all(n: usize) -> Vec<(usize, usize)> {
    (0..n)
        .step_by(2)
        .filter(|&i| i + 1 < n)
        .map(|i| (i, i + 1))
        .collect()
}

/// Every other disjoint pair: (0,1), (4,5), (8,9), ...
fn pairs_every_other(n: usize) -> Vec<(usize, usize)> {
    (0..n)
        .step_by(4)
        .filter(|&i| i + 1 < n)
        .map(|i| (i, i + 1))
        .collect()
}

macro_rules! bench_single_loop {
    ($group:expr, $tab:expr, $n:expr, $gran:expr, $method:ident, $idx:expr) => {{
        $group.bench_with_input(
            BenchmarkId::new(
                concat!("loop/", stringify!($method)),
                format!("n={}/{}", $n, $gran),
            ),
            $idx,
            |b, idx| {
                b.iter_batched_ref(
                    || $tab.clone(),
                    |t| {
                        for &q in idx.iter() {
                            t.$method(q);
                        }
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }};
}

macro_rules! bench_single_batch {
    ($group:expr, $tab:expr, $n:expr, $gran:expr, $method:ident, $idx:expr) => {{
        $group.bench_with_input(
            BenchmarkId::new(
                concat!("batch/", stringify!($method)),
                format!("n={}/{}", $n, $gran),
            ),
            $idx,
            |b, idx| {
                b.iter_batched_ref(
                    || $tab.clone(),
                    |t| t.$method(idx),
                    BatchSize::SmallInput,
                );
            },
        );
    }};
}

macro_rules! bench_pair_loop {
    ($group:expr, $tab:expr, $n:expr, $gran:expr, $method:ident, $pairs:expr) => {{
        $group.bench_with_input(
            BenchmarkId::new(
                concat!("loop/", stringify!($method)),
                format!("n={}/{}", $n, $gran),
            ),
            $pairs,
            |b, pairs| {
                b.iter_batched_ref(
                    || $tab.clone(),
                    |t| {
                        for &(c, x) in pairs.iter() {
                            t.$method(c, x);
                        }
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }};
}

macro_rules! bench_pair_batch {
    ($group:expr, $tab:expr, $n:expr, $gran:expr, $method:ident, $pairs:expr) => {{
        $group.bench_with_input(
            BenchmarkId::new(
                concat!("batch/", stringify!($method)),
                format!("n={}/{}", $n, $gran),
            ),
            $pairs,
            |b, pairs| {
                b.iter_batched_ref(
                    || $tab.clone(),
                    |t| t.$method(pairs),
                    BatchSize::SmallInput,
                );
            },
        );
    }};
}

fn bench_clifford_batch(c: &mut Criterion) {
    let mut group = c.benchmark_group("clifford-batch");
    for &n in SIZES {
        let tab = setup(n);
        let all = indices_all(n);
        let half = indices_every_other(n);
        let pairs_full = pairs_all(n);
        let pairs_half = pairs_every_other(n);

        // --- Clifford trait ---
        // x, y, z, s, cnot have no batched impl today.
        // h and cz do — bench both paths.
        bench_single_loop!(group, tab, n, "all", x, &all);
        bench_single_loop!(group, tab, n, "every_other", x, &half);
        bench_single_loop!(group, tab, n, "all", y, &all);
        bench_single_loop!(group, tab, n, "every_other", y, &half);
        bench_single_loop!(group, tab, n, "all", z, &all);
        bench_single_loop!(group, tab, n, "every_other", z, &half);
        bench_single_loop!(group, tab, n, "all", h, &all);
        bench_single_loop!(group, tab, n, "every_other", h, &half);
        bench_single_batch!(group, tab, n, "all", h_batch, &all);
        bench_single_batch!(group, tab, n, "every_other", h_batch, &half);
        bench_single_loop!(group, tab, n, "all", s, &all);
        bench_single_loop!(group, tab, n, "every_other", s, &half);

        bench_pair_loop!(group, tab, n, "all", cnot, &pairs_full);
        bench_pair_loop!(group, tab, n, "every_other", cnot, &pairs_half);
        bench_pair_loop!(group, tab, n, "all", cz, &pairs_full);
        bench_pair_loop!(group, tab, n, "every_other", cz, &pairs_half);
        bench_pair_batch!(group, tab, n, "all", cz_batch, &pairs_full);
        bench_pair_batch!(group, tab, n, "every_other", cz_batch, &pairs_half);

        // --- CliffordExtensions trait ---
        // s_adj, cy have no batched impl today.
        // sqrt_x, sqrt_x_adj, sqrt_y, sqrt_y_adj do — bench both paths.
        bench_single_loop!(group, tab, n, "all", s_adj, &all);
        bench_single_loop!(group, tab, n, "every_other", s_adj, &half);
        bench_single_loop!(group, tab, n, "all", sqrt_x, &all);
        bench_single_loop!(group, tab, n, "every_other", sqrt_x, &half);
        bench_single_batch!(group, tab, n, "all", sqrt_x_batch, &all);
        bench_single_batch!(group, tab, n, "every_other", sqrt_x_batch, &half);
        bench_single_loop!(group, tab, n, "all", sqrt_x_adj, &all);
        bench_single_loop!(group, tab, n, "every_other", sqrt_x_adj, &half);
        bench_single_batch!(group, tab, n, "all", sqrt_x_adj_batch, &all);
        bench_single_batch!(group, tab, n, "every_other", sqrt_x_adj_batch, &half);
        bench_single_loop!(group, tab, n, "all", sqrt_y, &all);
        bench_single_loop!(group, tab, n, "every_other", sqrt_y, &half);
        bench_single_batch!(group, tab, n, "all", sqrt_y_batch, &all);
        bench_single_batch!(group, tab, n, "every_other", sqrt_y_batch, &half);
        bench_single_loop!(group, tab, n, "all", sqrt_y_adj, &all);
        bench_single_loop!(group, tab, n, "every_other", sqrt_y_adj, &half);
        bench_single_batch!(group, tab, n, "all", sqrt_y_adj_batch, &all);
        bench_single_batch!(group, tab, n, "every_other", sqrt_y_adj_batch, &half);

        bench_pair_loop!(group, tab, n, "all", cy, &pairs_full);
        bench_pair_loop!(group, tab, n, "every_other", cy, &pairs_half);
    }
    group.finish();
}

criterion_group! {
    name = benches;
    config = configure();
    targets = bench_clifford_batch
}
criterion_main!(benches);
