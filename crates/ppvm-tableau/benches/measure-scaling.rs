use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<Byte8F64<2>, u128>;

fn build_circuit(n_qubits: usize, n_t_gates: usize) -> Tab {
    let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
    for i in 0..n_t_gates {
        tab.h(i);
        tab.t(i);
    }
    for i in 0..n_qubits - 1 {
        tab.cz(i, i + 1);
    }
    tab
}

fn bench_measure_scaling(c: &mut Criterion) {
    let n_qubits = 85;
    let mut group = c.benchmark_group("measure-scaling");

    for n_t in [8, 10, 12, 14] {
        let tab = build_circuit(n_qubits, n_t);
        group.bench_function(format!("t{n_t}-{n_qubits}q"), |b| {
            b.iter_batched_ref(
                || tab.fork(Some(42)),
                |t| {
                    for q in 0..n_qubits {
                        let _ = t.measure(q);
                    }
                },
                criterion::BatchSize::SmallInput,
            );
        });
    }

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_secs(2))
        .measurement_time(Duration::from_secs(5))
        .sample_size(50);
    targets = bench_measure_scaling
}
criterion_main!(benches);
