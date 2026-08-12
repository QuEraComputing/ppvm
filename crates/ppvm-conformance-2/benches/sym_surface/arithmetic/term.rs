// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_conformance_2::sym::NewSymKey;
use ppvm_sym::Term as OldTerm;
use ppvm_sym_2::Term as NewTerm;

use super::super::{assert_real, new_sum, old_sum};

fn old_fixture() -> OldTerm {
    (0..8).fold(OldTerm::from(2.0), |term, var| {
        term + OldTerm::var(var).sin() + OldTerm::var(var).cos()
    })
}

fn new_fixture() -> NewTerm {
    (0..8).fold(NewTerm::from(2.0), |term, var| {
        term + NewTerm::var(var).sin() + NewTerm::var(var).cos()
    })
}

pub(super) fn bench(c: &mut Criterion) {
    bench_term(c);
    bench_pauli_sum_coefficient(c);
}

fn bench_term(c: &mut Criterion) {
    let old = old_fixture();
    let new = new_fixture();
    let vals = [0.3; 8];
    let old_rhs = OldTerm::var(7).sin() + OldTerm::from(0.5);
    let new_rhs = NewTerm::var(7).sin() + NewTerm::from(0.5);

    assert_real(
        (old.clone() + old_rhs.clone()).eval(&vals).unwrap(),
        (new.clone() + new_rhs.clone()).eval(&vals).unwrap(),
    );
    let mut group = c.benchmark_group("sym/surface/term");
    group.bench_function("new/add", |b| {
        b.iter_batched(
            || (new.clone(), new_rhs.clone()),
            |(lhs, rhs)| lhs + rhs,
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/add", |b| {
        b.iter_batched(
            || (old.clone(), old_rhs.clone()),
            |(lhs, rhs)| lhs + rhs,
            BatchSize::SmallInput,
        )
    });

    assert_real(
        (old.clone() * old_rhs.clone()).eval(&vals).unwrap(),
        (new.clone() * new_rhs.clone()).eval(&vals).unwrap(),
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

    assert_real(
        (old.clone() * 1.5).eval(&vals).unwrap(),
        (new.clone() * 1.5).eval(&vals).unwrap(),
    );
    group.bench_function("new/mul_scalar", |b| {
        b.iter_batched(|| new.clone(), |term| term * 1.5, BatchSize::SmallInput)
    });
    group.bench_function("old/mul_scalar", |b| {
        b.iter_batched(|| old.clone(), |term| term * 1.5, BatchSize::SmallInput)
    });
    group.finish();
}

fn bench_pauli_sum_coefficient(c: &mut Criterion) {
    let mut old = old_sum(4);
    old += ("ZIII", OldTerm::var(0).sin());
    let mut new = new_sum(4);
    new += (NewSymKey::from("ZIII"), NewTerm::var(0).sin());
    let mut old_expected = old.clone();
    old_expected *= OldTerm::from(1.5);
    let mut new_expected = new.clone();
    new_expected *= NewTerm::from(1.5);
    assert_eq!(old_expected.data().len(), new_expected.len());

    let mut group = c.benchmark_group("sym/surface/pauli_sum");
    group.bench_function("new/mul_coefficient", |b| {
        b.iter_batched(
            || new.clone(),
            |mut sum| {
                sum *= NewTerm::from(1.5);
                sum
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/mul_coefficient", |b| {
        b.iter_batched(
            || old.clone(),
            |mut sum| {
                sum *= OldTerm::from(1.5);
                sum
            },
            BatchSize::SmallInput,
        )
    });
    group.finish();
}
