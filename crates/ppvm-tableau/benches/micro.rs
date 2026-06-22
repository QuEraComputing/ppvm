// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::time::Duration;

use bnum::types::U256;
use criterion::{Criterion, criterion_group, criterion_main};
use num::complex::Complex64;
use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<Byte8F64<2>, usize>;

const N_QUBITS: usize = 32;

fn configure() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50)
}

fn make_sparse_vec<I: TableauIndex>(n: usize) -> Vec<(Complex64, I)> {
    let mut vec: Vec<(Complex64, I)> = SparseVector::new();
    for k in 0..n {
        let index = <I as From<u8>>::from(k as u8) << 2;
        // Mix large and small entries so trim benchmarks exercise actual removal
        let value = if k % 4 == 0 {
            Complex64::new(0.001, 0.0) // small — will be trimmed at cutoff 0.05
        } else {
            Complex64::new(1.0 / (k as f64 + 1.0), 0.1 * k as f64)
        };
        vec.unsafe_insert(index, value);
    }
    vec
}

fn bench_single_qubit_gates(c: &mut Criterion) {
    let tab = Tab::new(N_QUBITS, 1e-10);
    let mut group = c.benchmark_group("gates/single-qubit");

    group.bench_function("h", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.h(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("s", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.s(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("s_dag", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.s_dag(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("x", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.x(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("y", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.y(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("z", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.z(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("sqrt_x", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.sqrt_x(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("sqrt_x_dag", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.sqrt_x_dag(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("sqrt_y", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.sqrt_y(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("sqrt_y_dag", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.sqrt_y_dag(0),
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}

fn bench_two_qubit_gates(c: &mut Criterion) {
    let tab = Tab::new(N_QUBITS, 1e-10);
    let mut group = c.benchmark_group("gates/two-qubit");

    group.bench_function("cnot", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.cnot([0, 1]),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("cz", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.cz([0, 1]),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("cy", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.cy([0, 1]),
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}

fn bench_non_clifford_gates(c: &mut Criterion) {
    // |+> state: Z-axis rotations (t, t_dag, u3) anticommute with X stabilizer → branching
    let mut tab_plus = Tab::new(N_QUBITS, 1e-10);
    tab_plus.h(0);

    // |0> state: X-axis rotations (rx) anticommute with Z stabilizer → branching
    let tab_zero = Tab::new(N_QUBITS, 1e-10);

    // |0,0> state: XX anticommutes with Z⊗I stabilizer → branching for rxx
    let tab_zero_2q = Tab::new(N_QUBITS, 1e-10);

    let mut group = c.benchmark_group("gates/non-clifford");

    let pi_4 = std::f64::consts::FRAC_PI_4;

    group.bench_function("t", |b| {
        b.iter_batched_ref(
            || tab_plus.fork(None),
            |t| t.t(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("t_dag", |b| {
        b.iter_batched_ref(
            || tab_plus.fork(None),
            |t| t.t_dag(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("rx", |b| {
        b.iter_batched_ref(
            || tab_zero.fork(None),
            |t| t.rx(0, pi_4),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("rxx", |b| {
        b.iter_batched_ref(
            || tab_zero_2q.fork(None),
            |t| t.rxx([0, 1], pi_4),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("u3", |b| {
        b.iter_batched_ref(
            || tab_plus.fork(None),
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
            |t| t.depolarize1(0, 1.0),
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
            |t| t.two_qubit_pauli_error([0, 1], [1.0 / 15.0; 15]),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("depolarize2", |b| {
        b.iter_batched_ref(
            || tab.fork(None),
            |t| t.depolarize2([0, 1], 1.0),
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

fn bench_sparse_vec_for_type<I: TableauIndex + std::fmt::Debug>(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
    type_name: &str,
) {
    let vec16 = make_sparse_vec::<I>(16);
    let existing_index = <I as From<u8>>::from(4u8) << 2;
    let new_index = <I as From<u8>>::from(99u8);

    group.bench_function(format!("{type_name}/unsafe_insert"), |b| {
        b.iter_batched_ref(
            || vec16.clone(),
            |v| v.unsafe_insert(new_index, Complex64::new(1.0, 0.0)),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function(format!("{type_name}/add_or_insert/existing"), |b| {
        b.iter_batched_ref(
            || vec16.clone(),
            |v| v.add_or_insert(existing_index, Complex64::new(0.5, 0.0)),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function(format!("{type_name}/add_or_insert/new"), |b| {
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

criterion_group! {
    name = benches;
    config = configure();
    targets = bench_single_qubit_gates, bench_two_qubit_gates, bench_non_clifford_gates,
              bench_measurement, bench_noise, bench_sparse_vec
}
criterion_main!(benches);
