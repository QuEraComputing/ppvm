// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hint::black_box;

use criterion::Criterion;
use ppvm_conformance_2::mixture::{New, Old, new, old};
use ppvm_traits::traits::LossChannel as OldLoss;
use ppvm_traits_2::LossChannel as NewLoss;

use super::support::{CUTOFF, SEED, assert_same, bench_mut, branch_pair, pair};

pub fn register(c: &mut Criterion) {
    constructors(c);

    let (old, new) = branch_pair(4);
    assert_same(&old.clone(), &new.clone());
    let mut clone = c.benchmark_group("mixture/lifecycle/clone");
    clone.bench_function("old", |b| b.iter(|| black_box(old.clone())));
    clone.bench_function("new", |b| b.iter(|| black_box(new.clone())));
    clone.finish();

    let mut len = c.benchmark_group("mixture/lifecycle/len");
    assert_eq!(old.len(), new.len());
    len.bench_function("old", |b| b.iter(|| black_box(old.len())));
    len.bench_function("new", |b| b.iter(|| black_box(new.len())));
    len.finish();

    let mut empty = c.benchmark_group("mixture/lifecycle/is_empty");
    assert_eq!(old.is_empty(), new.is_empty());
    empty.bench_function("old", |b| b.iter(|| black_box(old.is_empty())));
    empty.bench_function("new", |b| b.iter(|| black_box(new.is_empty())));
    empty.finish();

    normalize(c);
    truncate(c);
}

fn constructors(c: &mut Criterion) {
    let old_check = old(SEED, CUTOFF);
    let new_check = new(SEED, CUTOFF);
    assert_same(&old_check, &new_check);
    let mut seeded = c.benchmark_group("mixture/lifecycle/new_with_seed");
    seeded.bench_function("old", |b| b.iter(|| black_box(old(SEED, CUTOFF))));
    seeded.bench_function("new", |b| b.iter(|| black_box(new(SEED, CUTOFF))));
    seeded.finish();

    let old_check = Old::new(12, 1e-12, CUTOFF);
    let new_check = New::new(12, 1e-12, CUTOFF);
    assert_same(&old_check, &new_check);
    let mut entropy = c.benchmark_group("mixture/lifecycle/new");
    entropy.bench_function("old", |b| b.iter(|| black_box(Old::new(12, 1e-12, CUTOFF))));
    entropy.bench_function("new", |b| b.iter(|| black_box(New::new(12, 1e-12, CUTOFF))));
    entropy.finish();
}

fn normalize(c: &mut Criterion) {
    let (mut old, mut new) = branch_pair(2);
    old.entries.entries[0].1 = 3.0;
    old.entries.entries[1].1 = 1.0;
    new.entries[0].1 = 3.0;
    new.entries[1].1 = 1.0;
    assert_same(&old, &new);
    bench_mut(
        c,
        "mixture/lifecycle/normalize_probabilities",
        &old,
        &new,
        Old::normalize_probabilities,
        New::normalize_probabilities,
    );
}

fn truncate(c: &mut Criterion) {
    let (mut old_state, mut new_state) = (old(SEED, 0.1), new(SEED, 0.1));
    old_state.loss_channel(0, 0.2);
    new_state.loss_channel(0, 0.2);
    old_state.entries.entries[0].1 = 0.95;
    old_state.entries.entries[1].1 = 0.05;
    new_state.entries[0].1 = 0.95;
    new_state.entries[1].1 = 0.05;
    assert_same(&old_state, &new_state);
    bench_mut(
        c,
        "mixture/lifecycle/truncate_cutoff",
        &old_state,
        &new_state,
        Old::truncate,
        New::truncate,
    );

    // Keep the empty-cutoff boundary represented outside the timed path.
    let (old_empty, new_empty) = (old(SEED, 1.0), new(SEED, 1.0));
    assert!(old_empty.is_empty() && new_empty.is_empty());
    let _ = pair();
}
