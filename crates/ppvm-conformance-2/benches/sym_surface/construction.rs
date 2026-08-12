// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_sym::{Prod as OldProd, Sum as OldSum, Term as OldTerm};
use ppvm_sym_2::{Prod as NewProd, Sum as NewSum, Term as NewTerm};

use super::assert_real;

pub(super) fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("sym/surface/construct");

    assert_eq!(
        OldProd::new().eval(&[]).unwrap(),
        NewProd::new().eval(&[]).unwrap()
    );
    group.bench_function("new/prod_new", |b| b.iter(NewProd::new));
    group.bench_function("old/prod_new", |b| b.iter(OldProd::new));

    assert_real(
        OldProd::sin(7).eval(&[0.0; 8]).unwrap(),
        NewProd::sin(7).eval(&[0.0; 8]).unwrap(),
    );
    group.bench_function("new/prod_sin", |b| b.iter(|| NewProd::sin(7)));
    group.bench_function("old/prod_sin", |b| b.iter(|| OldProd::sin(7)));

    assert_real(
        OldProd::cos(7).eval(&[0.0; 8]).unwrap(),
        NewProd::cos(7).eval(&[0.0; 8]).unwrap(),
    );
    group.bench_function("new/prod_cos", |b| b.iter(|| NewProd::cos(7)));
    group.bench_function("old/prod_cos", |b| b.iter(|| OldProd::cos(7)));

    assert_eq!(
        OldSum::new().eval(&[]).unwrap(),
        NewSum::new().eval(&[]).unwrap()
    );
    group.bench_function("new/sum_new", |b| b.iter(NewSum::new));
    group.bench_function("old/sum_new", |b| b.iter(OldSum::new));

    assert_eq!(
        OldTerm::from_f64(2.5).eval(&[]).unwrap(),
        NewTerm::from_f64(2.5).eval(&[]).unwrap()
    );
    group.bench_function("new/term_constant", |b| b.iter(|| NewTerm::from_f64(2.5)));
    group.bench_function("old/term_constant", |b| b.iter(|| OldTerm::from_f64(2.5)));

    assert_eq!(
        OldTerm::var(7).eval(&[0.0; 8]).unwrap(),
        NewTerm::var(7).eval(&[0.0; 8]).unwrap()
    );
    group.bench_function("new/term_variable", |b| b.iter(|| NewTerm::var(7)));
    group.bench_function("old/term_variable", |b| b.iter(|| OldTerm::var(7)));

    assert_real(
        OldTerm::var(7).sin().eval(&[0.3; 8]).unwrap(),
        NewTerm::var(7).sin().eval(&[0.3; 8]).unwrap(),
    );
    group.bench_function("new/promote_sin", |b| b.iter(|| NewTerm::var(7).sin()));
    group.bench_function("old/promote_sin", |b| b.iter(|| OldTerm::var(7).sin()));

    assert_real(
        OldTerm::var(7).cos().eval(&[0.3; 8]).unwrap(),
        NewTerm::var(7).cos().eval(&[0.3; 8]).unwrap(),
    );
    group.bench_function("new/promote_cos", |b| b.iter(|| NewTerm::var(7).cos()));
    group.bench_function("old/promote_cos", |b| b.iter(|| OldTerm::var(7).cos()));

    assert_real(
        OldTerm::from_f64(0.3).sin().eval(&[]).unwrap(),
        NewTerm::from_f64(0.3).sin().eval(&[]).unwrap(),
    );
    group.bench_function("new/fold_sin_constant", |b| {
        b.iter(|| NewTerm::from_f64(0.3).sin())
    });
    group.bench_function("old/fold_sin_constant", |b| {
        b.iter(|| OldTerm::from_f64(0.3).sin())
    });

    assert_real(
        OldTerm::from_f64(0.3).cos().eval(&[]).unwrap(),
        NewTerm::from_f64(0.3).cos().eval(&[]).unwrap(),
    );
    group.bench_function("new/fold_cos_constant", |b| {
        b.iter(|| NewTerm::from_f64(0.3).cos())
    });
    group.bench_function("old/fold_cos_constant", |b| {
        b.iter(|| OldTerm::from_f64(0.3).cos())
    });

    group.finish();
}
