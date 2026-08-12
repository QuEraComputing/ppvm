// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_sym::{Prod as OldProd, Sum as OldSum};
use ppvm_sym_2::{Prod as NewProd, Sum as NewSum};

use super::super::assert_real;

fn old_fixture() -> OldSum {
    let mut sum = OldSum::new();
    sum.add_const(2.0, f64::EPSILON);
    for var in 0..16 {
        sum.add_term(OldProd::sin(var), 1.0, usize::MAX, f64::EPSILON);
    }
    sum
}

fn new_fixture() -> NewSum {
    let mut sum = NewSum::new();
    sum.add_const(2.0, f64::EPSILON);
    for var in 0..16 {
        sum.add_term(NewProd::sin(var), 1.0, usize::MAX, f64::EPSILON);
    }
    sum
}

pub(super) fn bench(c: &mut Criterion) {
    let old = old_fixture();
    let new = new_fixture();
    let vals = [0.3; 16];

    let mut old_added = old.clone();
    let mut new_added = new.clone();
    old_added.add_const(1.5, f64::EPSILON);
    new_added.add_const(1.5, f64::EPSILON);
    assert_real(
        old_added.eval(&vals).unwrap(),
        new_added.eval(&vals).unwrap(),
    );

    let mut group = c.benchmark_group("sym/surface/sum");
    group.bench_function("new/add_const", |b| {
        b.iter_batched(
            || new.clone(),
            |mut sum| {
                sum.add_const(1.5, f64::EPSILON);
                sum
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/add_const", |b| {
        b.iter_batched(
            || old.clone(),
            |mut sum| {
                sum.add_const(1.5, f64::EPSILON);
                sum
            },
            BatchSize::SmallInput,
        )
    });

    let mut old_added = old.clone();
    let mut new_added = new.clone();
    old_added.add_term(OldProd::cos(15), 1.5, usize::MAX, f64::EPSILON);
    new_added.add_term(NewProd::cos(15), 1.5, usize::MAX, f64::EPSILON);
    assert_real(
        old_added.eval(&vals).unwrap(),
        new_added.eval(&vals).unwrap(),
    );
    group.bench_function("new/add_term", |b| {
        b.iter_batched(
            || new.clone(),
            |mut sum| {
                sum.add_term(NewProd::cos(15), 1.5, usize::MAX, f64::EPSILON);
                sum
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/add_term", |b| {
        b.iter_batched(
            || old.clone(),
            |mut sum| {
                sum.add_term(OldProd::cos(15), 1.5, usize::MAX, f64::EPSILON);
                sum
            },
            BatchSize::SmallInput,
        )
    });

    let mut old_multiplied = old.clone();
    let mut new_multiplied = new.clone();
    old_multiplied.mul_term(OldProd::cos(15), 1.5, usize::MAX, f64::EPSILON);
    new_multiplied.mul_term(NewProd::cos(15), 1.5, usize::MAX, f64::EPSILON);
    assert_real(
        old_multiplied.eval(&vals).unwrap(),
        new_multiplied.eval(&vals).unwrap(),
    );
    group.bench_function("new/mul_term", |b| {
        b.iter_batched(
            || new.clone(),
            |mut sum| {
                sum.mul_term(NewProd::cos(15), 1.5, usize::MAX, f64::EPSILON);
                sum
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/mul_term", |b| {
        b.iter_batched(
            || old.clone(),
            |mut sum| {
                sum.mul_term(OldProd::cos(15), 1.5, usize::MAX, f64::EPSILON);
                sum
            },
            BatchSize::SmallInput,
        )
    });

    let mut old_scaled = old.clone();
    let mut new_scaled = new.clone();
    old_scaled *= 1.5;
    new_scaled *= 1.5;
    assert_real(
        old_scaled.eval(&vals).unwrap(),
        new_scaled.eval(&vals).unwrap(),
    );
    group.bench_function("new/mul_scalar", |b| {
        b.iter_batched(
            || new.clone(),
            |mut sum| {
                sum *= 1.5;
                sum
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/mul_scalar", |b| {
        b.iter_batched(
            || old.clone(),
            |mut sum| {
                sum *= 1.5;
                sum
            },
            BatchSize::SmallInput,
        )
    });
    group.finish();
}
