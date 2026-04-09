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
                        tableau.cnot(i, i + 1);
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
            tab.cnot(i, i + 1);
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
}

pub fn tableau_scaling_benchmarks(c: &mut Criterion) {
    benchmark_suite_tableau(c, "ByteF64FxIndexMap<8, CoefficientThreshold>");
}

criterion_group!(benches, tableau_scaling_benchmarks);
criterion_main!(benches);
