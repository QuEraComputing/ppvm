// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hint::black_box;

use criterion::{BatchSize, BenchmarkId, Criterion};

use super::support::{assert_same, branch_pair};

#[path = "sampling_parallel.rs"]
mod parallel;

const SHOT_COUNTS: [usize; 4] = [1, 16, 128, 1024];

pub fn register(c: &mut Criterion) {
    let (old, new) = branch_pair(8);

    let (mut old_check, mut new_check) = (old.clone(), new.clone());
    let mut old_sampler = old_check.sampler();
    let mut new_sampler = new_check.sampler();
    assert_eq!(old_sampler.entries.len(), new_sampler.entries.len());
    assert_eq!(old_sampler.sample(), new_sampler.sample());
    assert_same(&old_check, &new_check);

    let mut construction = c.benchmark_group("mixture/sampler/construction");
    construction.bench_function("old", |b| {
        b.iter_batched(
            || old.clone(),
            |mut state| black_box(state.sampler()),
            BatchSize::SmallInput,
        )
    });
    construction.bench_function("new", |b| {
        b.iter_batched(
            || new.clone(),
            |mut state| black_box(state.sampler()),
            BatchSize::SmallInput,
        )
    });
    construction.finish();

    let mut sample = c.benchmark_group("mixture/sampler/single_sample");
    sample.bench_function("old", |b| {
        b.iter_batched(
            || {
                let mut state = old.clone();
                state.sampler()
            },
            |mut sampler| black_box(sampler.sample()),
            BatchSize::SmallInput,
        )
    });
    sample.bench_function("new", |b| {
        b.iter_batched(
            || {
                let mut state = new.clone();
                state.sampler()
            },
            |mut sampler| black_box(sampler.sample()),
            BatchSize::SmallInput,
        )
    });
    sample.finish();

    adaptive_dispatch(c, &old, &new);
    shot_scaling(c, &old, &new);
    parallel::register(c);
}

fn adaptive_dispatch(
    c: &mut Criterion,
    old: &ppvm_conformance_2::mixture::Old,
    new: &ppvm_conformance_2::mixture::New,
) {
    let (mut old_check, mut new_check) = (old.clone(), new.clone());
    assert_eq!(
        old_check.sampler().sample_shots(128),
        new_check.sampler().sample_shots(128),
    );
    let mut group = c.benchmark_group("mixture/sampler/adaptive_sample_shots_128");
    group.bench_function("old", |b| {
        b.iter_batched(
            || {
                let mut state = old.clone();
                state.sampler()
            },
            |mut sampler| black_box(sampler.sample_shots(128)),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("new", |b| {
        b.iter_batched(
            || {
                let mut state = new.clone();
                state.sampler()
            },
            |mut sampler| black_box(sampler.sample_shots(128)),
            BatchSize::SmallInput,
        )
    });
    group.finish();
}

fn shot_scaling(
    c: &mut Criterion,
    old: &ppvm_conformance_2::mixture::Old,
    new: &ppvm_conformance_2::mixture::New,
) {
    let mut group = c.benchmark_group("mixture/sampler/serial_shot_scaling");
    for shots in SHOT_COUNTS {
        let (mut old_check, mut new_check) = (old.clone(), new.clone());
        assert_eq!(
            old_check.sampler().sample_shots_serial(shots),
            new_check.sampler().sample_shots_serial(shots),
        );
        group.bench_with_input(BenchmarkId::new("old", shots), &shots, |b, &shots| {
            b.iter_batched(
                || {
                    let mut state = old.clone();
                    state.sampler()
                },
                |mut sampler| black_box(sampler.sample_shots_serial(shots)),
                BatchSize::SmallInput,
            )
        });
        group.bench_with_input(BenchmarkId::new("new", shots), &shots, |b, &shots| {
            b.iter_batched(
                || {
                    let mut state = new.clone();
                    state.sampler()
                },
                |mut sampler| black_box(sampler.sample_shots_serial(shots)),
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}
