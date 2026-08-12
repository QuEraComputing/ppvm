// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hint::black_box;

use criterion::Criterion;

use super::{NewKey, OldKey, assert_pair, build_new, build_old, terms};

pub fn bench(c: &mut Criterion) {
    let data = terms(0, 192);
    let (old, new) = (build_old(&data), build_new(&data));
    assert_pair(&old, &new);
    let word = &data[64].0;
    let (old_key, new_key) = (OldKey::from(word.as_str()), NewKey::from(word.as_str()));

    let mut group = c.benchmark_group("pauli_sum_surface/inspect/get");
    group.bench_function("old", |b| {
        b.iter(|| black_box(old.data().get(black_box(&old_key))))
    });
    group.bench_function("new", |b| {
        b.iter(|| black_box(new.get(black_box(&new_key))))
    });
    group.finish();

    let mut group = c.benchmark_group("pauli_sum_surface/inspect/contains_key");
    group.bench_function("old", |b| {
        b.iter(|| black_box(old.data().contains_key(black_box(&old_key))))
    });
    group.bench_function("new", |b| {
        b.iter(|| black_box(new.contains_key(black_box(&new_key))))
    });
    group.finish();

    let old_coeff = *old.data().get(&old_key).expect("prepared key is present");
    let new_coeff = new.get(&new_key).expect("prepared key is present");
    assert_eq!(
        old.contains(&old_key, &old_coeff),
        new.contains(&new_key, &new_coeff)
    );
    assert!(old.contains(&old_key, &old_coeff));
    let mut group = c.benchmark_group("pauli_sum_surface/inspect/contains_key_value");
    group.bench_function("old", |b| {
        b.iter(|| black_box(old.contains(black_box(&old_key), black_box(&old_coeff))))
    });
    group.bench_function("new", |b| {
        b.iter(|| black_box(new.contains(black_box(&new_key), black_box(&new_coeff))))
    });
    group.finish();

    let (old_equal, new_equal) = (old.clone(), new.clone());
    assert!(old == old_equal);
    assert!(new == new_equal);
    let mut group = c.benchmark_group("pauli_sum_surface/inspect/equality_equal_support");
    group.bench_function("old", |b| {
        b.iter(|| black_box(black_box(&old) == black_box(&old_equal)))
    });
    group.bench_function("new", |b| {
        b.iter(|| black_box(black_box(&new) == black_box(&new_equal)))
    });
    group.finish();

    let mut group = c.benchmark_group("pauli_sum_surface/inspect/metadata");
    group.bench_function("old", |b| {
        b.iter(|| black_box((old.n_qubits(), old.capacity(), old.len(), old.is_empty())))
    });
    group.bench_function("new", |b| {
        b.iter(|| black_box((new.n_sites(), new.capacity(), new.len(), new.is_empty())))
    });
    group.finish();

    let expected: f64 = old.iter().map(|(_, c)| *c).sum();
    let mut got = 0.0;
    new.for_each_ref(|_, c| got += c);
    assert!((expected - got).abs() <= 1e-9);
    let mut group = c.benchmark_group("pauli_sum_surface/inspect/borrowed_traversal");
    group.bench_function("old", |b| {
        b.iter(|| black_box(old.iter().fold(0.0, |acc, (_, c)| acc + *c)))
    });
    group.bench_function("new", |b| {
        b.iter(|| {
            let mut total = 0.0;
            new.for_each_ref(|_, c| total += c);
            black_box(total)
        })
    });
    group.finish();
}
