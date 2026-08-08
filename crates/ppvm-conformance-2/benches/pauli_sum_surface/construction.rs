// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hint::black_box;

use criterion::{BatchSize, Criterion};
use ppvm_traits::traits::NoStrategy;

use super::{
    CAPACITY, N, NewKey, NewSum, OldKey, OldSum, assert_pair, build_new, build_old, keyed_new,
    keyed_old, terms,
};

pub fn bench(c: &mut Criterion) {
    constructor(c);
    parse(c);
    build(c);
    clone_sum(c);
    add_term(c);
    add_sum(c);
    extend(c);
}

fn constructor(c: &mut Criterion) {
    let mut group = c.benchmark_group("pauli_sum_surface/construct/empty");
    group.bench_function("old", |b| {
        b.iter(|| {
            black_box(
                OldSum::builder()
                    .n_qubits(N)
                    .strategy(NoStrategy)
                    .capacity(CAPACITY)
                    .build(),
            )
        })
    });
    group.bench_function("new", |b| {
        b.iter(|| {
            black_box(NewSum::with_capacity(
                N,
                ppvm_pauli_sum_2::NoPolicy,
                CAPACITY,
            ))
        })
    });
    group.finish();
}

fn parse(c: &mut Criterion) {
    let text = "XYZIYZXI";
    let old = OldKey::from(text);
    let new = NewKey::from(text);
    assert_eq!(old.to_string(), new.to_string());
    let mut group = c.benchmark_group("pauli_sum_surface/construct/parse_word");
    group.bench_function("old", |b| {
        b.iter(|| black_box(OldKey::from(black_box(text))))
    });
    group.bench_function("new", |b| {
        b.iter(|| black_box(NewKey::from(black_box(text))))
    });
    group.finish();
}

fn build(c: &mut Criterion) {
    let data = terms(0, 192);
    let old_terms = keyed_old(&data);
    let new_terms = keyed_new(&data);
    assert_pair(&build_old(&data), &build_new(&data));
    let mut group = c.benchmark_group("pauli_sum_surface/construct/build_support");
    group.bench_function("old", |b| {
        b.iter(|| {
            let mut sum = OldSum::builder()
                .n_qubits(N)
                .strategy(NoStrategy)
                .capacity(CAPACITY)
                .build();
            for (key, coeff) in &old_terms {
                sum += (*key, *coeff);
            }
            black_box(sum)
        })
    });
    group.bench_function("new", |b| {
        b.iter(|| {
            let mut sum = NewSum::with_capacity(N, ppvm_pauli_sum_2::NoPolicy, CAPACITY);
            for (key, coeff) in &new_terms {
                sum += (key.clone(), *coeff);
            }
            black_box(sum)
        })
    });
    group.finish();
}

fn clone_sum(c: &mut Criterion) {
    let data = terms(0, 192);
    let (old, new) = (build_old(&data), build_new(&data));
    assert_pair(&old, &new);
    let mut group = c.benchmark_group("pauli_sum_surface/construct/clone");
    group.bench_function("old", |b| b.iter(|| black_box(old.clone())));
    group.bench_function("new", |b| b.iter(|| black_box(new.clone())));
    group.finish();
}

fn add_term(c: &mut Criterion) {
    let data = terms(0, 192);
    let (old, new) = (build_old(&data), build_new(&data));
    let old_term = (OldKey::from("XYZIYZXI"), 0.75);
    let new_term = (NewKey::from("XYZIYZXI"), 0.75);
    let mut op = old.clone();
    op += old_term;
    let mut np = new.clone();
    np.add_term(new_term.0.clone(), new_term.1);
    assert_pair(&op, &np);
    let mut group = c.benchmark_group("pauli_sum_surface/add/term");
    group.bench_function("old", |b| {
        b.iter_batched_ref(|| old.clone(), |s| *s += old_term, BatchSize::LargeInput)
    });
    group.bench_function("new", |b| {
        b.iter_batched_ref(
            || new.clone(),
            |s| s.add_term(new_term.0.clone(), new_term.1),
            BatchSize::LargeInput,
        )
    });
    group.finish();
}

fn add_sum(c: &mut Criterion) {
    let (lhs, rhs) = (terms(0, 96), terms(256, 96));
    let (old, old_rhs) = (build_old(&lhs), build_old(&rhs));
    let (new, new_rhs) = (build_new(&lhs), build_new(&rhs));
    let mut op = old.clone();
    op += old_rhs.clone();
    let mut np = new.clone();
    np.add_sum(&new_rhs);
    assert_pair(&op, &np);
    let mut group = c.benchmark_group("pauli_sum_surface/add/sum_disjoint");
    group.bench_function("old", |b| {
        b.iter_batched(
            || (old.clone(), old_rhs.clone()),
            |(mut s, rhs)| {
                s += rhs;
                black_box(s)
            },
            BatchSize::LargeInput,
        )
    });
    group.bench_function("new", |b| {
        b.iter_batched_ref(
            || new.clone(),
            |s| s.add_sum(black_box(&new_rhs)),
            BatchSize::LargeInput,
        )
    });
    group.finish();
}

fn extend(c: &mut Criterion) {
    let data = terms(0, 192);
    let (old, new) = (build_old(&data), build_new(&data));
    let extra = terms(512, 32);
    let (old_extra, new_extra) = (keyed_old(&extra), keyed_new(&extra));
    let mut op = old.clone();
    op.extend(old_extra.clone());
    let mut np = new.clone();
    np.extend(new_extra.clone());
    assert_pair(&op, &np);
    let mut group = c.benchmark_group("pauli_sum_surface/add/extend");
    group.bench_function("old", |b| {
        b.iter_batched_ref(
            || old.clone(),
            |s| s.extend(old_extra.clone()),
            BatchSize::LargeInput,
        )
    });
    group.bench_function("new", |b| {
        b.iter_batched_ref(
            || new.clone(),
            |s| s.extend(new_extra.clone()),
            BatchSize::LargeInput,
        )
    });
    group.finish();
}
