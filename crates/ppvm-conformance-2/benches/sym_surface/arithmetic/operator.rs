// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Public operator spellings kept separate from `add_const` / `add_term`.

use criterion::{BatchSize, Criterion};
use ppvm_sym::{Prod as OldProd, Sum as OldSum, Term as OldTerm};
use ppvm_sym_2::{Prod as NewProd, Sum as NewSum, Term as NewTerm};

use super::super::assert_real;

pub(super) fn bench(c: &mut Criterion) {
    let vals = [0.3; 8];
    let old_sum = OldSum::new();
    let new_sum = NewSum::new();

    let mut old_expected = old_sum.clone();
    let mut new_expected = new_sum.clone();
    old_expected += 1.5;
    new_expected += 1.5;
    assert_real(
        old_expected.eval(&vals).unwrap(),
        new_expected.eval(&vals).unwrap(),
    );

    let mut group = c.benchmark_group("sym/surface/operator_add");
    group.bench_function("new/sum_add_coefficient", |b| {
        b.iter_batched(
            || new_sum.clone(),
            |mut sum| {
                sum += 1.5;
                sum
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/sum_add_coefficient", |b| {
        b.iter_batched(
            || old_sum.clone(),
            |mut sum| {
                sum += 1.5;
                sum
            },
            BatchSize::SmallInput,
        )
    });

    let mut old_expected = old_sum.clone();
    let mut new_expected = new_sum.clone();
    old_expected += OldProd::sin(7);
    new_expected += NewProd::sin(7);
    assert_real(
        old_expected.eval(&vals).unwrap(),
        new_expected.eval(&vals).unwrap(),
    );
    group.bench_function("new/sum_add_term", |b| {
        b.iter_batched(
            || new_sum.clone(),
            |mut sum| {
                sum += NewProd::sin(7);
                sum
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/sum_add_term", |b| {
        b.iter_batched(
            || old_sum.clone(),
            |mut sum| {
                sum += OldProd::sin(7);
                sum
            },
            BatchSize::SmallInput,
        )
    });

    // A map-backed receiver avoids old's documented single-monomial `+= f64`
    // bug, keeping this an honest comparable path.
    let old_term = OldTerm::from(1.0) + OldTerm::var(0).sin();
    let new_term = NewTerm::from(1.0) + NewTerm::var(0).sin();
    let mut old_expected = old_term.clone();
    let mut new_expected = new_term.clone();
    old_expected += 1.5;
    new_expected += 1.5;
    assert_real(
        old_expected.eval(&vals).unwrap(),
        new_expected.eval(&vals).unwrap(),
    );
    group.bench_function("new/term_add_coefficient", |b| {
        b.iter_batched(
            || new_term.clone(),
            |mut term| {
                term += 1.5;
                term
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/term_add_coefficient", |b| {
        b.iter_batched(
            || old_term.clone(),
            |mut term| {
                term += 1.5;
                term
            },
            BatchSize::SmallInput,
        )
    });

    let old_rhs = OldTerm::var(1).cos();
    let new_rhs = NewTerm::var(1).cos();
    assert_real(
        (old_term.clone() - old_rhs.clone()).eval(&vals).unwrap(),
        (new_term.clone() - new_rhs.clone()).eval(&vals).unwrap(),
    );
    group.bench_function("new/term_subtract_term", |b| {
        b.iter_batched(
            || (new_term.clone(), new_rhs.clone()),
            |(lhs, rhs)| lhs - rhs,
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/term_subtract_term", |b| {
        b.iter_batched(
            || (old_term.clone(), old_rhs.clone()),
            |(lhs, rhs)| lhs - rhs,
            BatchSize::SmallInput,
        )
    });

    assert_real(
        (old_term.clone() - 1.5).eval(&vals).unwrap(),
        (new_term.clone() - 1.5).eval(&vals).unwrap(),
    );
    group.bench_function("new/term_subtract_coefficient", |b| {
        b.iter_batched(
            || new_term.clone(),
            |term| term - 1.5,
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/term_subtract_coefficient", |b| {
        b.iter_batched(
            || old_term.clone(),
            |term| term - 1.5,
            BatchSize::SmallInput,
        )
    });

    assert_real(
        (-old_term.clone()).eval(&vals).unwrap(),
        (-new_term.clone()).eval(&vals).unwrap(),
    );
    group.bench_function("new/term_negate", |b| {
        b.iter_batched(|| new_term.clone(), |term| -term, BatchSize::SmallInput)
    });
    group.bench_function("old/term_negate", |b| {
        b.iter_batched(|| old_term.clone(), |term| -term, BatchSize::SmallInput)
    });
    group.finish();
}
