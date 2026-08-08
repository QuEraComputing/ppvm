// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hint::black_box;

use criterion::{BatchSize, Criterion};
use ppvm_conformance_2::mixture::{New, new};
use ppvm_traits_2::RotXY;

use super::support::{CUTOFF, SEED};

pub fn register(c: &mut Criterion) {
    // NEW-ONLY, NO OLD TWIN: old exposes its store field but has no mixture
    // iterator method.
    let iter_state = new(SEED, CUTOFF);
    assert_eq!(iter_state.iter().count(), iter_state.len());
    let mut iter = c.benchmark_group("mixture/new_only_no_old_twin/iter");
    iter.bench_function("new", |b| b.iter(|| black_box(iter_state.iter().count())));
    iter.finish();

    // NEW-ONLY, NO OLD TWIN: the old mixture does not implement `RotXY`.
    let state = new(SEED, CUTOFF);
    let mut r = c.benchmark_group("mixture/new_only_no_old_twin/rot_xy_r");
    r.bench_function("new", |b| {
        b.iter_batched(
            || state.clone(),
            |mut state: New| {
                state.r(5, 0.3, -0.4);
                black_box(state)
            },
            BatchSize::SmallInput,
        )
    });
    r.finish();
}
