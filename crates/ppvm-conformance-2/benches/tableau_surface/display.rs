// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;

use super::*;

pub fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("tableau-surface/display");
    let (old_bare, new_bare) = prepared_bare(16);
    let (old_gen, new_gen) = prepared_gen(16);

    assert_eq!(old_bare.to_string(), new_bare.to_string());
    group.bench_function("bare/old", |b| {
        b.iter(|| std::hint::black_box(format!("{old_bare}")))
    });
    group.bench_function("bare/new", |b| {
        b.iter(|| std::hint::black_box(format!("{new_bare}")))
    });

    assert_eq!(old_gen.to_string(), new_gen.to_string());
    group.bench_function("generalized/old", |b| {
        b.iter(|| std::hint::black_box(format!("{old_gen}")))
    });
    group.bench_function("generalized/new", |b| {
        b.iter(|| std::hint::black_box(format!("{new_gen}")))
    });
    group.finish();
}
