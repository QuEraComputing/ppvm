// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_sym::Term as OldTerm;
use ppvm_sym_2::Term as NewTerm;

use super::super::assert_real;

pub(super) fn bench(c: &mut Criterion) {
    let old = OldTerm::from(1.0) + OldTerm::var(1).cos();
    let new = NewTerm::from(1.0) + NewTerm::var(1).cos();

    let mut old_max = old.clone();
    let mut new_max = new.clone();
    old_max.set_max_sin(0);
    new_max.set_max_sin(0);
    old_max *= OldTerm::var(0).sin();
    new_max *= NewTerm::var(0).sin();
    assert_real(
        old_max.eval(&[0.3, -0.7]).unwrap(),
        new_max.eval(&[0.3, -0.7]).unwrap(),
    );
    assert_eq!(new_max.max_sin(), 0);

    let mut old_eps = old.clone();
    let mut new_eps = new.clone();
    old_eps.set_min_eps(1e-3);
    new_eps.set_min_eps(1e-3);
    old_eps += OldTerm::var(0).sin() * 1e-9;
    new_eps += NewTerm::var(0).sin() * 1e-9;
    assert_real(
        old_eps.eval(&[0.3, -0.7]).unwrap(),
        new_eps.eval(&[0.3, -0.7]).unwrap(),
    );
    assert_eq!(new_eps.min_eps(), 1e-3);

    let mut group = c.benchmark_group("sym/surface/term_setters");
    group.bench_function("new/set_max_sin", |b| {
        b.iter_batched(
            || new.clone(),
            |mut term| {
                term.set_max_sin(3);
                term
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/set_max_sin", |b| {
        b.iter_batched(
            || old.clone(),
            |mut term| {
                term.set_max_sin(3);
                term
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("new/set_min_eps", |b| {
        b.iter_batched(
            || new.clone(),
            |mut term| {
                term.set_min_eps(1e-9);
                term
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/set_min_eps", |b| {
        b.iter_batched(
            || old.clone(),
            |mut term| {
                term.set_min_eps(1e-9);
                term
            },
            BatchSize::SmallInput,
        )
    });
    group.finish();
}
