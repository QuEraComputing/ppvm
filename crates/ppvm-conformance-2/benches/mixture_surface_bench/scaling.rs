// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hint::black_box;

use criterion::{BatchSize, BenchmarkId, Criterion};
use ppvm_traits::traits::{Clifford as OldClifford, LossChannel as OldLoss};
use ppvm_traits_2::{Clifford as NewClifford, LossChannel as NewLoss};

use super::support::{assert_same, branch_pair};

pub fn register(c: &mut Criterion) {
    clifford_branch_scaling(c);
    branching_channel_scaling(c);
    sampler_branch_scaling(c);
}

fn clifford_branch_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("mixture/scaling/branches_h");
    for branches in [1usize, 2, 4, 8, 16] {
        let (old, new) = branch_pair(branches);
        let (mut old_check, mut new_check) = (old.clone(), new.clone());
        old_check.h(7);
        new_check.h(7);
        assert_same(&old_check, &new_check);
        group.bench_with_input(BenchmarkId::new("old", branches), &branches, |b, _| {
            b.iter_batched(
                || old.clone(),
                |mut state| {
                    state.h(7);
                    black_box(state)
                },
                BatchSize::SmallInput,
            )
        });
        group.bench_with_input(BenchmarkId::new("new", branches), &branches, |b, _| {
            b.iter_batched(
                || new.clone(),
                |mut state| {
                    state.h(7);
                    black_box(state)
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}

fn branching_channel_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("mixture/scaling/branches_loss");
    for branches in [1usize, 2, 4, 8, 16] {
        let (old, new) = branch_pair(branches);
        let qubit = branches.ilog2() as usize + 1;
        let (mut old_check, mut new_check) = (old.clone(), new.clone());
        old_check.loss_channel(qubit, 0.2);
        new_check.loss_channel(qubit, 0.2);
        assert_eq!(old_check.len(), 2 * branches);
        assert_eq!(new_check.len(), 2 * branches);
        assert_same(&old_check, &new_check);
        group.bench_with_input(BenchmarkId::new("old", branches), &branches, |b, _| {
            b.iter_batched(
                || old.clone(),
                |mut state| {
                    state.loss_channel(qubit, 0.2);
                    black_box(state)
                },
                BatchSize::SmallInput,
            )
        });
        group.bench_with_input(BenchmarkId::new("new", branches), &branches, |b, _| {
            b.iter_batched(
                || new.clone(),
                |mut state| {
                    state.loss_channel(qubit, 0.2);
                    black_box(state)
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}

fn sampler_branch_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("mixture/scaling/branches_sampler_construction");
    for branches in [1usize, 2, 4, 8, 16] {
        let (old, new) = branch_pair(branches);
        let (mut old_check, mut new_check) = (old.clone(), new.clone());
        let mut old_sampler = old_check.sampler();
        let mut new_sampler = new_check.sampler();
        assert_eq!(old_sampler.entries.len(), new_sampler.entries.len());
        assert_eq!(old_sampler.sample(), new_sampler.sample());
        group.bench_with_input(BenchmarkId::new("old", branches), &branches, |b, _| {
            b.iter_batched(
                || old.clone(),
                |mut state| black_box(state.sampler()),
                BatchSize::SmallInput,
            )
        });
        group.bench_with_input(BenchmarkId::new("new", branches), &branches, |b, _| {
            b.iter_batched(
                || new.clone(),
                |mut state| black_box(state.sampler()),
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}
