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

fn bench_non_clifford_gates(c: &mut Criterion) {
    let mut tab = Tab::new(N_QUBITS, 1e-10);
    tab.h(0);

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

fn bench_measurement(c: &mut Criterion) {
    let tab_det = Tab::new(N_QUBITS, 1e-10);

    let mut tab_rand = Tab::new(N_QUBITS, 1e-10);
    tab_rand.h(0);

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

criterion_group! {
    name = benches;
    config = configure();
    targets = bench_single_qubit_gates, bench_two_qubit_gates, bench_non_clifford_gates,
              bench_measurement, bench_noise
}
criterion_main!(benches);
