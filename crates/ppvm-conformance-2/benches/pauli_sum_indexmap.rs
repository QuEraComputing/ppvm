// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Same-build legacy/new insertion-ordered backend comparison.

use std::hint::black_box;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use ppvm_conformance_2::{random_terms, seeded_rng};
use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;
use ppvm_pauli_sum::sum::PauliSum as OldSumT;
use ppvm_pauli_sum_2::{IndexPauliSum, PauliWord};
use ppvm_traits::traits::{Clifford as OldClifford, RotationOne as OldRotation};
use ppvm_traits_2::{Clifford as NewClifford, RotationOne as NewRotation};

type Old = OldSumT<ByteFxHashF64<8>>;
type New = IndexPauliSum<8>;

const N: usize = 16;
const TERMS: usize = 1_000;

fn terms() -> Vec<(String, f64)> {
    random_terms(&mut seeded_rng(0x1d_ea), N, TERMS)
}

fn old_sum(terms: &[(String, f64)]) -> Old {
    let mut sum = Old::builder().n_qubits(N).capacity(TERMS * 2).build();
    for (word, coeff) in terms {
        sum += (word.as_str(), *coeff);
    }
    sum
}

fn new_sum(terms: &[(String, f64)]) -> New {
    let mut sum = New::with_capacity(N, ppvm_pauli_sum_2::NoPolicy, TERMS * 2);
    for (word, coeff) in terms {
        sum += (PauliWord::from(word.as_str()), *coeff);
    }
    sum
}

fn bench_build_and_terms(c: &mut Criterion) {
    let terms = terms();
    let mut group = c.benchmark_group("pauli_sum_indexmap/build");
    group.bench_function("old", |b| b.iter(|| black_box(old_sum(black_box(&terms)))));
    group.bench_function("new", |b| b.iter(|| black_box(new_sum(black_box(&terms)))));
    group.finish();

    let (old, new) = (old_sum(&terms), new_sum(&terms));
    let mut group = c.benchmark_group("pauli_sum_indexmap/ordered_terms");
    group.bench_function("old", |b| {
        b.iter(|| {
            black_box(
                old.iter()
                    .map(|(k, c)| (k.to_string(), *c))
                    .collect::<Vec<_>>(),
            )
        })
    });
    group.bench_function("new", |b| {
        b.iter(|| {
            black_box(
                new.iter()
                    .map(|(k, c)| (k.to_string(), c))
                    .collect::<Vec<_>>(),
            )
        })
    });
    group.finish();
}

fn bench_gates(c: &mut Criterion) {
    let terms = terms();
    let (old, new) = (old_sum(&terms), new_sum(&terms));
    let mut group = c.benchmark_group("pauli_sum_indexmap/gates");
    group.bench_function("old/cnot", |b| {
        let mut sum = old.clone();
        b.iter(|| sum.cnot(black_box(0), black_box(1)));
    });
    group.bench_function("new/cnot", |b| {
        let mut sum = new.clone();
        b.iter(|| sum.cnot(black_box(0), black_box(1)));
    });
    group.bench_function("old/rx", |b| {
        b.iter_batched_ref(
            || old.clone(),
            |sum| sum.rx(black_box(0), black_box(0.37)),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("new/rx", |b| {
        b.iter_batched_ref(
            || new.clone(),
            |sum| sum.rx(black_box(0), black_box(0.37)),
            BatchSize::SmallInput,
        )
    });
    group.finish();
}

criterion_group!(benches, bench_build_and_terms, bench_gates);
criterion_main!(benches);
