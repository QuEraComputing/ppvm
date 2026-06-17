// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;

pub fn benchmark_suite_tableau(c: &mut Criterion, name: impl AsRef<str>) {
    let mut group = c.benchmark_group(name.as_ref());

    for n_qubits in (32..129).step_by(32) {
        let tableau = GeneralizedTableau::<Byte8F64<2>, usize>::new(n_qubits, 1e-10);

        group.bench_function(format!("tableau-scaling-{}", n_qubits), |b| {
            b.iter_batched_ref(
                || tableau.fork(None),
                |tableau| {
                    tableau.h(0);
                    tableau.t(0);
                    for i in 0..n_qubits - 1 {
                        tableau.cnot([i, i + 1]);
                    }

                    // some more T gates
                    tableau.t(n_qubits - 1);
                    tableau.t(n_qubits - 2);

                    tableau.measure(0)
                },
                criterion::BatchSize::SmallInput,
            );
        });

        let mut tab = GeneralizedTableau::<Byte8F64<2>, usize>::new(n_qubits, 1e-10);
        tab.h(0);
        tab.t(0);
        for i in 0..n_qubits - 1 {
            tab.cnot([i, i + 1]);
        }

        // some more T gates
        tab.t(n_qubits - 1);
        tab.t(n_qubits - 2);
        group.bench_function(format!("tableau-measure-scaling-{}", n_qubits), |b| {
            b.iter_batched_ref(
                || tab.fork(None),
                |tab| {
                    for i in 0..n_qubits {
                        tab.measure(i);
                    }
                },
                criterion::BatchSize::SmallInput,
            );
        });
    }

    let mut tableau = GeneralizedTableau::<Byte8F64<2>, usize>::new(85, 1e-10);
    for i in 0..10 {
        // make sure it branches with t gates
        tableau.h(i);
    }
    for tgate_num in [1, 4, 8, 12] {
        group.bench_function(format!("tableau-t-gate-{}", tgate_num), |b| {
            b.iter_batched_ref(
                || tableau.fork(None),
                |tableau| {
                    for i in 0..tgate_num {
                        tableau.t(i);
                    }
                },
                criterion::BatchSize::SmallInput,
            );
        });

        group.bench_function(format!("tableau-measure-t-gate-{}", tgate_num), |b| {
            b.iter_batched_ref(
                || {
                    let mut tab = tableau.fork(None);
                    for i in 0..tgate_num {
                        tab.t(i);
                    }
                    tab
                },
                |tab| {
                    for i in 0..10 {
                        tab.measure(i);
                    }
                },
                criterion::BatchSize::SmallInput,
            );
        });
    }

    group.finish();

    // Large-coefficient T-gate benchmarks: measures the cost of a single T gate
    // applied to a pre-built state with a known number of coefficients.
    // This isolates the coefficient branching hot path for rayon testing.
    let mut large_group = c.benchmark_group("large-tgate");
    large_group
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(5))
        .sample_size(20);

    // Use 128 qubits (u128 index) so we have enough room for many distinct T targets
    type LargeTab = GeneralizedTableau<Byte8F64<2>, u128>;

    for n_tgates in [12, 14, 16, 18] {
        // Pre-build the state with (n_tgates - 1) T gates already applied,
        // then benchmark applying the last one.
        let mut setup: LargeTab = GeneralizedTableau::new(128, 1e-10);
        for i in 0..n_tgates {
            setup.h(i);
        }
        for i in 0..n_tgates - 1 {
            setup.t(i);
        }
        let last_qubit = n_tgates - 1;
        let n_coeffs = setup.coefficients.len();

        large_group.bench_function(format!("single-t-on-{n_coeffs}-coeffs"), |b| {
            b.iter_batched_ref(
                || setup.fork(Some(0)),
                |tab| tab.t(last_qubit),
                criterion::BatchSize::LargeInput,
            );
        });
    }

    large_group.finish();

    // Combined benchmark: fused Clifford batches + variable T gates + measurement.
    // Shows how fusion and rayon compose in a realistic circuit.
    let mut combined_group = c.benchmark_group("fused-tgate-circuit");
    combined_group
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(5))
        .sample_size(20);

    for n_tgates in [8, 12, 16] {
        // Circuit: batched Cliffords → H (ensure branching) → T gates → batched Cliffords → measure.
        // H is applied AFTER Cliffords to ensure T gates see non-stabilizer states and branch.
        let n_qubits: usize = 85;

        // Pre-build the Clifford portion so the setup cost is excluded
        let mut setup = GeneralizedTableau::<Byte8F64<2>, u128>::new(n_qubits, 1e-10);
        let block1: Vec<usize> = (0..17).collect();
        let block2: Vec<usize> = (17..34).collect();
        setup.sqrt_y_batch(&block1);
        setup.sqrt_x_batch(&block2);
        setup.cz_block_pairs(0, 17, 17);
        // H on qubits that will get T gates — applied after Cliffords to guarantee branching
        for i in 0..n_tgates {
            setup.h(i);
        }

        combined_group.bench_function(format!("fused-{n_tgates}t-{n_qubits}q"), |b| {
            b.iter_batched_ref(
                || setup.fork(Some(0)),
                |tab| {
                    // T gates (coefficient branching — rayon target)
                    for i in 0..n_tgates {
                        tab.t(i);
                    }

                    // More batched Cliffords after T gates
                    tab.sqrt_x_adj_batch(&block1);
                    tab.sqrt_y_adj_batch(&block2);

                    // Measure all
                    for i in 0..n_qubits {
                        tab.measure(i);
                    }
                },
                criterion::BatchSize::LargeInput,
            );
        });
    }

    combined_group.finish();
}

pub fn tableau_scaling_benchmarks(c: &mut Criterion) {
    benchmark_suite_tableau(c, "ByteF64FxIndexMap<8, CoefficientThreshold>");
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50);
    targets = tableau_scaling_benchmarks
}
criterion_main!(benches);
