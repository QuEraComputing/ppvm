// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::{collections::BTreeMap, hint::black_box};

use criterion::{BatchSize, BenchmarkId, Criterion};

use super::super::support::branch_pair;

const SHOT_COUNTS: [usize; 4] = [1, 16, 128, 1024];
const CROSSOVER_SHOTS: [usize; 11] = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];

pub fn register(c: &mut Criterion) {
    parallel_scaling(c);
    serial_parallel_crossover(c);
}

fn parallel_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("mixture/sampler/parallel_branch_shot_scaling");
    for branches in [1usize, 8, 64] {
        let (old, new) = branch_pair(branches);
        for shots in SHOT_COUNTS {
            assert_parallel_parity(&old, &new, shots);
            let parameter = format!("{branches}_branches/{shots}_shots");
            group.bench_with_input(BenchmarkId::new("old", &parameter), &shots, |b, &shots| {
                b.iter_batched(
                    || {
                        let mut state = old.clone();
                        state.sampler()
                    },
                    |mut sampler| black_box(sampler.sample_shots_parallel(shots)),
                    BatchSize::SmallInput,
                )
            });
            group.bench_with_input(BenchmarkId::new("new", &parameter), &shots, |b, &shots| {
                b.iter_batched(
                    || {
                        let mut state = new.clone();
                        state.sampler()
                    },
                    |mut sampler| black_box(sampler.sample_shots_parallel(shots)),
                    BatchSize::SmallInput,
                )
            });
        }
    }
    group.finish();
}

fn serial_parallel_crossover(c: &mut Criterion) {
    let (old, new) = branch_pair(8);
    let mut group = c.benchmark_group("mixture/sampler/serial_parallel_crossover");
    for shots in CROSSOVER_SHOTS {
        assert_parallel_parity(&old, &new, shots);
        group.bench_with_input(
            BenchmarkId::new("old_serial", shots),
            &shots,
            |b, &shots| {
                b.iter_batched(
                    || {
                        let mut state = old.clone();
                        state.sampler()
                    },
                    |mut sampler| black_box(sampler.sample_shots_serial(shots)),
                    BatchSize::SmallInput,
                )
            },
        );
        group.bench_with_input(
            BenchmarkId::new("old_parallel", shots),
            &shots,
            |b, &shots| {
                b.iter_batched(
                    || {
                        let mut state = old.clone();
                        state.sampler()
                    },
                    |mut sampler| black_box(sampler.sample_shots_parallel(shots)),
                    BatchSize::SmallInput,
                )
            },
        );
        group.bench_with_input(
            BenchmarkId::new("new_serial", shots),
            &shots,
            |b, &shots| {
                b.iter_batched(
                    || {
                        let mut state = new.clone();
                        state.sampler()
                    },
                    |mut sampler| black_box(sampler.sample_shots_serial(shots)),
                    BatchSize::SmallInput,
                )
            },
        );
        group.bench_with_input(
            BenchmarkId::new("new_parallel", shots),
            &shots,
            |b, &shots| {
                b.iter_batched(
                    || {
                        let mut state = new.clone();
                        state.sampler()
                    },
                    |mut sampler| black_box(sampler.sample_shots_parallel(shots)),
                    BatchSize::SmallInput,
                )
            },
        );
    }
    group.finish();
}

fn assert_parallel_parity(
    old: &ppvm_conformance_2::mixture::Old,
    new: &ppvm_conformance_2::mixture::New,
    shots: usize,
) {
    let (mut old_serial, mut old_parallel) = (old.clone(), old.clone());
    let (mut new_serial, mut new_parallel) = (new.clone(), new.clone());
    let old_serial = old_serial.sampler().sample_shots_serial(shots);
    let old_parallel = old_parallel.sampler().sample_shots_parallel(shots);
    let new_serial = new_serial.sampler().sample_shots_serial(shots);
    let new_parallel = new_parallel.sampler().sample_shots_parallel(shots);

    assert_eq!(old_serial, old_parallel);
    assert_eq!(new_serial, new_parallel);
    assert_eq!(distribution(old_parallel), distribution(new_parallel));
}

fn distribution(samples: Vec<Vec<Option<bool>>>) -> BTreeMap<Vec<Option<bool>>, usize> {
    let mut counts = BTreeMap::new();
    for sample in samples {
        *counts.entry(sample).or_default() += 1;
    }
    counts
}
