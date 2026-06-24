// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Two-qubit-rotation (`rotate_2`) throughput benchmark.
//!
//! `rotate_2` (RXX/RYY/RZZ) is the only caller of the "apply" coefficient
//! accumulation (`compute_coefficients_after_pauli_apply`). The headline
//! `stim-circuits` bench is T-gate heavy and never hits this path, so this bench
//! exercises it directly: a branchy brickwork of non-Clifford two-qubit
//! rotations whose coefficient vector grows into the thousands, making the
//! per-`rotate_2` apply cost dominate.
//!
//! Matches the cultivation config (`ByteFxHashF64<8>, usize`) so the workload is
//! representative of a real branchy run.

use std::f64::consts::PI;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<ByteFxHashF64<8>, usize>;

/// Branchy brickwork of two-qubit rotations on `n` qubits over `layers` layers.
fn rot2_brickwork(n: usize, layers: usize) -> Tab {
    let mut tab: Tab = GeneralizedTableau::new_with_seed(n, 1e-10, 1);
    for q in (0..n).step_by(2) {
        tab.h(q);
    }
    for layer in 0..layers {
        for a in (0..n.saturating_sub(1)).step_by(2) {
            tab.rxx(a, a + 1, 0.3 * PI);
            tab.ryy(a, a + 1, 0.4 * PI);
        }
        for a in (1..n.saturating_sub(1)).step_by(2) {
            tab.rzz(a, a + 1, 0.25 * PI);
            tab.rxx(a, a + 1, 0.15 * PI);
        }
        if layer % 2 == 0 {
            for q in (1..n).step_by(2) {
                tab.h(q);
            }
        }
    }
    tab
}

pub fn rot2_apply_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("rot2-apply");
    // (qubits, layers): each grows the coefficient vector to a different scale.
    for &(n, layers) in &[(8usize, 4usize), (10, 4), (12, 3)] {
        let m = rot2_brickwork(n, layers).coefficients.len();
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("n{n}_l{layers}_m{m}")),
            &(n, layers),
            |b, &(n, layers)| b.iter(|| rot2_brickwork(n, layers)),
        );
    }
    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(2))
        .sample_size(30);
    targets = rot2_apply_benchmarks
}
criterion_main!(benches);
