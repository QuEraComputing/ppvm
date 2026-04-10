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

criterion_group! {
    name = benches;
    config = configure();
    targets = bench_single_qubit_gates
}
criterion_main!(benches);
