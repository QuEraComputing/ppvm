// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_sym::Prod as OldProd;
use ppvm_sym_2::Prod as NewProd;

use super::super::assert_real;

pub(super) fn bench(c: &mut Criterion) {
    let mut old = OldProd::sin(0);
    let mut new = NewProd::sin(0);
    for var in 1..8 {
        old.mul_cos(var);
        new.mul_cos(var);
    }

    let mut old_sin = old.clone();
    let mut new_sin = new.clone();
    old_sin.mul_sin(3);
    new_sin.mul_sin(3);
    assert_real(
        old_sin.eval(&[0.3; 8]).unwrap(),
        new_sin.eval(&[0.3; 8]).unwrap(),
    );

    let mut group = c.benchmark_group("sym/surface/product");
    let mut old_phase = old.clone();
    let mut new_phase = new.clone();
    old_phase.add_phase(1);
    new_phase.add_phase(1);
    old_phase.add_phase(3);
    new_phase.add_phase(3);
    assert_eq!(old_phase, old);
    assert_eq!(new_phase, new);
    group.bench_function("new/add_phase", |b| {
        b.iter_batched(
            || new.clone(),
            |mut product| {
                product.add_phase(1);
                product
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/add_phase", |b| {
        b.iter_batched(
            || old.clone(),
            |mut product| {
                product.add_phase(1);
                product
            },
            BatchSize::SmallInput,
        )
    });

    group.bench_function("new/mul_sin", |b| {
        b.iter_batched(
            || new.clone(),
            |mut product| {
                product.mul_sin(3);
                product
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/mul_sin", |b| {
        b.iter_batched(
            || old.clone(),
            |mut product| {
                product.mul_sin(3);
                product
            },
            BatchSize::SmallInput,
        )
    });

    let mut old_cos = old.clone();
    let mut new_cos = new.clone();
    old_cos.mul_cos(3);
    new_cos.mul_cos(3);
    assert_real(
        old_cos.eval(&[0.3; 8]).unwrap(),
        new_cos.eval(&[0.3; 8]).unwrap(),
    );
    group.bench_function("new/mul_cos", |b| {
        b.iter_batched(
            || new.clone(),
            |mut product| {
                product.mul_cos(3);
                product
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/mul_cos", |b| {
        b.iter_batched(
            || old.clone(),
            |mut product| {
                product.mul_cos(3);
                product
            },
            BatchSize::SmallInput,
        )
    });

    let old_rhs = OldProd::sin(7);
    let new_rhs = NewProd::sin(7);
    assert_real(
        (old.clone() * old_rhs.clone()).eval(&[0.3; 8]).unwrap(),
        (new.clone() * new_rhs.clone()).eval(&[0.3; 8]).unwrap(),
    );
    group.bench_function("new/multiply", |b| {
        b.iter_batched(
            || (new.clone(), new_rhs.clone()),
            |(lhs, rhs)| lhs * rhs,
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/multiply", |b| {
        b.iter_batched(
            || (old.clone(), old_rhs.clone()),
            |(lhs, rhs)| lhs * rhs,
            BatchSize::SmallInput,
        )
    });
    group.finish();
}
