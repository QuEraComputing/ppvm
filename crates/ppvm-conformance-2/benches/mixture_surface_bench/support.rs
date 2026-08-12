// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hint::black_box;

use criterion::{BatchSize, Criterion};
use ppvm_conformance_2::mixture::{
    New, Old, assert_snapshots_close, new, new_snapshot, old, old_snapshot,
};
use ppvm_traits::traits::LossChannel as OldLossChannel;
use ppvm_traits_2::LossChannel as NewLossChannel;

pub const SEED: u64 = 0x5eed_cafe;
pub const CUTOFF: f64 = 1e-14;

pub fn pair() -> (Old, New) {
    let pair = (old(SEED, CUTOFF), new(SEED, CUTOFF));
    assert_same(&pair.0, &pair.1);
    pair
}

pub fn branch_pair(branches: usize) -> (Old, New) {
    assert!(branches.is_power_of_two());
    let (mut old, mut new) = pair();
    for qubit in 0..branches.ilog2() as usize {
        old.loss_channel(qubit, 0.2);
        new.loss_channel(qubit, 0.2, &mut ppvm_conformance_2::analytic_rng());
    }
    assert_eq!(old.len(), branches);
    assert_eq!(new.len(), branches);
    assert_same(&old, &new);
    (old, new)
}

pub fn assert_same(old: &Old, new: &New) {
    assert_snapshots_close(old_snapshot(old), new_snapshot(new));
}

pub fn bench_mut<OldOp, NewOp>(
    c: &mut Criterion,
    name: &str,
    old_base: &Old,
    new_base: &New,
    old_op: OldOp,
    new_op: NewOp,
) where
    OldOp: Fn(&mut Old) + Copy,
    NewOp: Fn(&mut New) + Copy,
{
    let (mut old_check, mut new_check) = (old_base.clone(), new_base.clone());
    old_op(&mut old_check);
    new_op(&mut new_check);
    assert_same(&old_check, &new_check);

    let mut group = c.benchmark_group(name);
    group.bench_function("old", |b| {
        b.iter_batched(
            || old_base.clone(),
            |mut state| {
                old_op(&mut state);
                black_box(state)
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("new", |b| {
        b.iter_batched(
            || new_base.clone(),
            |mut state| {
                new_op(&mut state);
                black_box(state)
            },
            BatchSize::SmallInput,
        )
    });
    group.finish();
}

pub fn bench_output<OldOp, NewOp, Output>(
    c: &mut Criterion,
    name: &str,
    old_base: &Old,
    new_base: &New,
    old_op: OldOp,
    new_op: NewOp,
) where
    OldOp: Fn(&mut Old) -> Output + Copy,
    NewOp: Fn(&mut New) -> Output + Copy,
    Output: PartialEq + std::fmt::Debug,
{
    let (mut old_check, mut new_check) = (old_base.clone(), new_base.clone());
    assert_eq!(old_op(&mut old_check), new_op(&mut new_check));
    assert_same(&old_check, &new_check);

    let mut group = c.benchmark_group(name);
    group.bench_function("old", |b| {
        b.iter_batched(
            || old_base.clone(),
            |mut state| black_box(old_op(&mut state)),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("new", |b| {
        b.iter_batched(
            || new_base.clone(),
            |mut state| black_box(new_op(&mut state)),
            BatchSize::SmallInput,
        )
    });
    group.finish();
}
