// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::BatchSize;

use super::*;

pub fn bench(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    old: &OldGen,
    new: &NewGen,
    targets: &[usize],
) {
    let mut oc = old.fork(Some(SEED + 2));
    let mut nc = new.fork(Some(SEED + 2));
    let mut os = ppvm_tableau::measure::MeasureScratch::<u128, f64>::new();
    let mut ns = ppvm_tableau_2::MeasureScratch::<u128>::new();
    assert_eq!(
        oc.measure_many_with_scratch(targets, &mut os),
        nc.measure_many_with_scratch(targets, &mut ns)
    );
    assert_gen_eq(&oc, &nc);

    group.bench_function("generalized/measure_many_with_scratch/old", |b| {
        b.iter_batched(
            || {
                (
                    old.fork(Some(SEED + 2)),
                    ppvm_tableau::measure::MeasureScratch::new(),
                )
            },
            |(mut t, mut s)| std::hint::black_box(t.measure_many_with_scratch(targets, &mut s)),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("generalized/measure_many_with_scratch/new", |b| {
        b.iter_batched(
            || {
                (
                    new.fork(Some(SEED + 2)),
                    ppvm_tableau_2::MeasureScratch::new(),
                )
            },
            |(mut t, mut s)| std::hint::black_box(t.measure_many_with_scratch(targets, &mut s)),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("generalized/measure_all_with_scratch/old", |b| {
        b.iter_batched(
            || {
                (
                    old.fork(Some(SEED + 2)),
                    ppvm_tableau::measure::MeasureScratch::new(),
                )
            },
            |(mut t, mut s)| std::hint::black_box(t.measure_all_with_scratch(&mut s)),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("generalized/measure_all_with_scratch/new", |b| {
        b.iter_batched(
            || {
                (
                    new.fork(Some(SEED + 2)),
                    ppvm_tableau_2::MeasureScratch::new(),
                )
            },
            |(mut t, mut s)| std::hint::black_box(t.measure_all_with_scratch(&mut s)),
            BatchSize::SmallInput,
        )
    });
}
