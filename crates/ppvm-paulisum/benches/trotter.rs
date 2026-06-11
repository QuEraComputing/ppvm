// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_runtime::prelude::*;
use ppvm_runtime::strategy::CoefficientThreshold;
use rayon::current_num_threads;

fn trotter_func<T: Config<Coeff = f64, Strategy = CoefficientThreshold>>(
    state: &mut PauliSum<T>,
    n: usize,
    total_time: &f64,
    dt: &f64,
    interaction_strength: &f64,
    external_field: &f64,
    noise_params: [f64; 3],
) {
    let steps = (total_time / dt) as usize;

    let theta_zz = dt * interaction_strength;
    let theta_x = dt * external_field;
    for _ in 0..steps {
        // perform trotter step

        // truncate after each gate application to be consistent with PP.jl
        for i in 0..n {
            state.rx(i, theta_x);
            state.truncate();

            state.pauli_error(i, noise_params);
            state.truncate();
        }
        for i in 0..n - 1 {
            state.rzz(i, i + 1, theta_zz);
            state.truncate();

            state.pauli_error(i, noise_params);
            state.truncate();

            state.pauli_error(i + 1, noise_params);
            state.truncate();
        }
    }
}

pub fn benchmark_suite_trotter<T: Config<Coeff = f64, Strategy = CoefficientThreshold>>(
    c: &mut Criterion,
    name: impl AsRef<str>,
) {
    let mut group = c.benchmark_group(name.as_ref());

    // parameters
    let n_qubits = 12;
    let h = 1.0;
    let dt = 0.1 / h;
    let time = 1.0 / h;
    let j = 1.0 / 8.0 * h;

    let strat = CoefficientThreshold(1e-6);
    let mut state: PauliSum<T> = PauliSum::builder()
        .n_qubits(n_qubits)
        .strategy(strat)
        .capacity(n_qubits.pow(2))
        .build();

    // initial state: let's calculate the expectation value of Sum(Z(i))
    let tmp = vec!['I'; n_qubits];
    for i in 0..n_qubits {
        let mut zi = tmp.clone();
        zi[i] = 'Z';
        let zi_string: String = zi.iter().collect();
        state += (zi_string, 1.0);
    }

    println!("Initial state has {} terms", state.len());

    let noise_params = [1e-4 / 4.0; 3];

    group.bench_function("trotter", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                trotter_func(state, n_qubits, &time, &dt, &j, &h, noise_params);
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}

pub fn trotter_benchmarks(c: &mut Criterion) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build_global()
        .unwrap();
    println!("Using {} threads", current_num_threads());
    benchmark_suite_trotter::<config::gxhash::ByteF64<2, CoefficientThreshold>>(
        c,
        "ByteF64GxHashMap<2, CoefficientThreshold>",
    );
    benchmark_suite_trotter::<config::fxhash::ByteF64<2, CoefficientThreshold>>(
        c,
        "ByteF64FxHashMap<2, CoefficientThreshold>",
    );
    benchmark_suite_trotter::<config::dashmap::ByteFxHashF64<2, CoefficientThreshold>>(
        c,
        "ByteF64FxDashMap<2, CoefficientThreshold>",
    );
    benchmark_suite_trotter::<config::dashmap::ByteGxHashF64<2, CoefficientThreshold>>(
        c,
        "ByteF64GxDashMap<2, CoefficientThreshold>",
    );
    benchmark_suite_trotter::<config::indexmap::ByteFxHashF64<2, CoefficientThreshold>>(
        c,
        "ByteF64FxIndexMap<2, CoefficientThreshold>",
    );
    benchmark_suite_trotter::<config::indexmap::ByteGxHashF64<2, CoefficientThreshold>>(
        c,
        "ByteF64GxIndexMap<2, CoefficientThreshold>",
    );
}

criterion_group!(benches, trotter_benchmarks);
criterion_main!(benches);
