// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_pauli_sum::prelude::*;
use rayon::current_num_threads;

fn layer<T: Config>(state: &mut PauliSum<T>, n: usize) {
    for i in 0..n {
        state.rz(i, 1.1);
        state.ry(i, 2.1);
        state.rz(i, 1.1);
    }
}

fn entangle<T: Config>(state: &mut PauliSum<T>, n: usize) {
    for i in 0..n {
        state.cnot(i, (i + 1) % n)
    }
}

fn random_circuit<T: Config>(state: &mut PauliSum<T>, n_qubits: usize, circuit_depth: usize) {
    for _ in 0..circuit_depth {
        layer(state, n_qubits);
        entangle(state, n_qubits);
    }

    layer(state, n_qubits);

    // let zero_state_pattern: PauliPattern = "Z?*".into();
    // state.trace(&zero_state_pattern);
}

pub fn benchmark_suite_random_circuit<T: Config>(c: &mut Criterion, name: impl AsRef<str>) {
    let mut group = c.benchmark_group(name.as_ref());

    // parameters
    let n_qubits = 4;
    let circuit_depth = 2;

    let mut state: PauliSum<T> = PauliSum::builder().n_qubits(n_qubits).build();

    // initial state: let's calculate the expectation value of Sum(Z(i))
    let mut term = T::PauliWordType::new(n_qubits);
    term.set(0, Pauli::Z);
    term.set(1, Pauli::Z);
    state += (term.clone(), T::Coeff::from(1.0));

    group.bench_function("random-circuit", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                random_circuit(state, n_qubits, circuit_depth);
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}

pub fn random_circuit_benchmarks(c: &mut Criterion) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build_global()
        .unwrap();
    println!("Using {} threads", current_num_threads());
    benchmark_suite_random_circuit::<config::gxhash::ByteF64<2>>(c, "ByteF64GxHashMap<2>");
    benchmark_suite_random_circuit::<config::fxhash::ByteF64<2>>(c, "ByteF64FxHashMap<2>");
    benchmark_suite_random_circuit::<config::dashmap::ByteFxHashF64<2>>(c, "ByteF64FxDashMap<2>");
    benchmark_suite_random_circuit::<config::dashmap::ByteGxHashF64<2>>(c, "ByteF64GxDashMap<2>");
    benchmark_suite_random_circuit::<config::indexmap::ByteFxHashF64<2>>(c, "ByteF64FxIndexMap<2>");
    benchmark_suite_random_circuit::<config::indexmap::ByteGxHashF64<2>>(c, "ByteF64GxIndexMap<2>");
}

criterion_group!(benches, random_circuit_benchmarks);
criterion_main!(benches);
