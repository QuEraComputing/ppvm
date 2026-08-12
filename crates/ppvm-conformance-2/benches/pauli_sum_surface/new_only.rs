// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hint::black_box;

use criterion::{BatchSize, Criterion};
use num::Complex;
use ppvm_pauli_sum_2::{HashMapStore, NoPolicy, Sum};

use super::{CAPACITY, N, NewKey, NewSum, build_new, terms};

type ComplexSum = Sum<HashMapStore<NewKey, Complex<f64>>, NoPolicy>;

pub fn bench(c: &mut Criterion) {
    reduce(c);
    hermitian(c);
    singleton_products(c);
    products(c);
}

fn reduce(c: &mut Criterion) {
    let mut data = terms(0, 192);
    let cancelled = "ZZZZZZZZ".to_string();
    data.push((cancelled.clone(), 1.0));
    data.push((cancelled, -1.0));
    let seed: NewSum = build_new(&data);
    let mut probe = seed.clone();
    let before = probe.len();
    probe.reduce();
    assert!(probe.len() < before);
    let mut group = c.benchmark_group("pauli_sum_surface/new_only/reduce");
    group.bench_function("new", |b| {
        b.iter_batched_ref(|| seed.clone(), |s| s.reduce(), BatchSize::LargeInput)
    });
    group.finish();
}

fn complex(offset: usize, count: usize) -> ComplexSum {
    let mut sum = ComplexSum::with_capacity(N, NoPolicy, CAPACITY);
    for (i, (word, coeff)) in terms(offset, count).into_iter().enumerate() {
        sum += (
            NewKey::from(word.as_str()),
            Complex::new(coeff, (i % 11) as f64 / 13.0 - 0.4),
        );
    }
    sum
}

fn hermitian(c: &mut Criterion) {
    let (a, b) = (complex(0, 192), complex(96, 192));
    let value = a.hermitian_overlap(&b);
    assert!((value.conj() - b.hermitian_overlap(&a)).norm() <= 1e-9);
    let mut group = c.benchmark_group("pauli_sum_surface/new_only/hermitian_overlap");
    group.bench_function("new", |bencher| {
        bencher.iter(|| black_box(a.hermitian_overlap(black_box(&b))))
    });
    group.finish();
}

fn singleton_products(c: &mut Criterion) {
    let lhs = complex(0, 48);
    let mut rhs = ComplexSum::with_capacity(N, NoPolicy, CAPACITY);
    rhs += (NewKey::from("XYZXYZXY"), Complex::new(0.75, -0.125));
    assert_eq!(rhs.len(), 1);
    let allocated = lhs.multiply(&rhs);
    let mut in_place = lhs.clone();
    in_place *= &rhs;
    assert_eq!(allocated, in_place);

    let mut group = c.benchmark_group("pauli_sum_surface/new_only/sum_product_single_rhs_allocate");
    group.bench_function("new", |bencher| {
        bencher.iter(|| black_box(lhs.multiply(black_box(&rhs))))
    });
    group.finish();

    let mut group =
        c.benchmark_group("pauli_sum_surface/new_only/sum_product_single_rhs_mul_assign");
    group.bench_function("new", |bencher| {
        bencher.iter_batched_ref(
            || lhs.clone(),
            |sum| *sum *= black_box(&rhs),
            BatchSize::LargeInput,
        )
    });
    group.finish();
}

fn products(c: &mut Criterion) {
    let (a, b) = (complex(0, 48), complex(128, 24));
    let expected = a.multiply(&b);
    assert!(!expected.is_empty());

    let mut group = c.benchmark_group("pauli_sum_surface/new_only/sum_product_allocate");
    group.bench_function("new", |bencher| {
        bencher.iter(|| black_box(a.multiply(black_box(&b))))
    });
    group.finish();

    let mut group = c.benchmark_group("pauli_sum_surface/new_only/sum_product_into");
    group.bench_function("new", |bencher| {
        bencher.iter_batched(
            || ComplexSum::with_capacity(N, NoPolicy, CAPACITY),
            |mut acc| {
                a.multiply_into(black_box(&b), &mut acc);
                black_box(acc)
            },
            BatchSize::LargeInput,
        )
    });
    group.finish();

    let mut group = c.benchmark_group("pauli_sum_surface/new_only/sum_product_in_place");
    group.bench_function("new", |bencher| {
        bencher.iter_batched_ref(
            || a.clone(),
            |sum| sum.multiply_in_place(black_box(&b)),
            BatchSize::LargeInput,
        )
    });
    group.finish();
}
