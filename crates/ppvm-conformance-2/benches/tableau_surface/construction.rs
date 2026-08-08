// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BenchmarkId, Criterion};
use ppvm_conformance_2::tableau::Driver;

use super::*;

pub fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("tableau-surface/construction");

    for n in WIDTHS {
        let (old_bare, new_bare) = prepared_bare(n);
        let (old_gen, new_gen) = prepared_gen(n);
        assert_bare_eq(&old_bare.clone(), &new_bare.clone());
        assert_gen_eq(&old_gen.clone(), &new_gen.clone());
        assert_gen_eq(&old_gen.fork(Some(SEED + 1)), &new_gen.fork(Some(SEED + 1)));
        let (old_entropy_bare, new_entropy_bare) = (OldBare::new(n), NewBare::new(n));
        assert_bare_eq(&old_entropy_bare, &new_entropy_bare);
        let (old_entropy_gen, new_entropy_gen) =
            (OldGen::new(n, THRESHOLD), NewGen::new(n, THRESHOLD));
        assert_gen_eq(&old_entropy_gen, &new_entropy_gen);

        group.bench_with_input(BenchmarkId::new("bare/new/old", n), &n, |b, &n| {
            b.iter(|| std::hint::black_box(OldBare::new_with_seed(n, SEED)))
        });
        group.bench_with_input(BenchmarkId::new("bare/new/new", n), &n, |b, &n| {
            b.iter(|| std::hint::black_box(NewBare::new_with_seed(n, SEED)))
        });
        group.bench_with_input(BenchmarkId::new("generalized/new/old", n), &n, |b, &n| {
            b.iter(|| std::hint::black_box(<OldGen as Driver>::new_seeded(n, THRESHOLD, SEED)))
        });
        group.bench_with_input(BenchmarkId::new("generalized/new/new", n), &n, |b, &n| {
            b.iter(|| std::hint::black_box(<NewGen as Driver>::new_seeded(n, THRESHOLD, SEED)))
        });
        group.bench_with_input(BenchmarkId::new("bare/new-entropy/old", n), &n, |b, &n| {
            b.iter(|| std::hint::black_box(OldBare::new(n)))
        });
        group.bench_with_input(BenchmarkId::new("bare/new-entropy/new", n), &n, |b, &n| {
            b.iter(|| std::hint::black_box(NewBare::new(n)))
        });
        group.bench_with_input(
            BenchmarkId::new("generalized/new-entropy/old", n),
            &n,
            |b, &n| b.iter(|| std::hint::black_box(OldGen::new(n, THRESHOLD))),
        );
        group.bench_with_input(
            BenchmarkId::new("generalized/new-entropy/new", n),
            &n,
            |b, &n| b.iter(|| std::hint::black_box(NewGen::new(n, THRESHOLD))),
        );

        group.bench_function(format!("bare/clone/{n}/old"), |b| {
            b.iter(|| std::hint::black_box(old_bare.clone()))
        });
        group.bench_function(format!("bare/clone/{n}/new"), |b| {
            b.iter(|| std::hint::black_box(new_bare.clone()))
        });
        group.bench_function(format!("generalized/clone/{n}/old"), |b| {
            b.iter(|| std::hint::black_box(old_gen.clone()))
        });
        group.bench_function(format!("generalized/clone/{n}/new"), |b| {
            b.iter(|| std::hint::black_box(new_gen.clone()))
        });
        group.bench_function(format!("generalized/fork/{n}/old"), |b| {
            b.iter(|| std::hint::black_box(old_gen.fork(Some(SEED + 1))))
        });
        group.bench_function(format!("generalized/fork/{n}/new"), |b| {
            b.iter(|| std::hint::black_box(new_gen.fork(Some(SEED + 1))))
        });

        let (mut old_check, mut new_check) = (old_gen.clone(), new_gen.clone());
        old_check.reset_all();
        new_check.reset_all();
        assert_gen_eq(&old_check, &new_check);
        bench_mut_pair!(
            group,
            format!("generalized/reset_all/{n}"),
            old_gen,
            new_gen,
            |t: &mut OldGen| t.reset_all(),
            |t: &mut NewGen| t.reset_all()
        );
        bench_mut_pair!(
            group,
            format!("bare/reset_all/{n}"),
            old_bare,
            new_bare,
            |t: &mut OldBare| t.reset_all(),
            |t: &mut NewBare| t.reset_all()
        );
    }
    group.finish();
}
