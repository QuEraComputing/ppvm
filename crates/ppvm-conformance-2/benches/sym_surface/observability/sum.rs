// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_sym::{Prod as OldProd, Sum as OldSum};
use ppvm_sym_2::{Prod as NewProd, Sum as NewSum};

use super::super::assert_real;

fn fixtures() -> (OldSum, NewSum) {
    let mut old = OldSum::new();
    old.add_const(2.0, f64::EPSILON);
    old.add_term(OldProd::sin(0), 1.5, usize::MAX, f64::EPSILON);
    old.add_term(OldProd::cos(1), -0.5, usize::MAX, f64::EPSILON);
    old.add_term(
        OldProd::sin(2) * OldProd::cos(2),
        0.25,
        usize::MAX,
        f64::EPSILON,
    );

    let mut new = NewSum::new();
    new.add_const(2.0, f64::EPSILON);
    new.add_term(NewProd::sin(0), 1.5, usize::MAX, f64::EPSILON);
    new.add_term(NewProd::cos(1), -0.5, usize::MAX, f64::EPSILON);
    new.add_term(
        NewProd::sin(2) * NewProd::cos(2),
        0.25,
        usize::MAX,
        f64::EPSILON,
    );
    (old, new)
}

pub(super) fn bench(c: &mut Criterion) {
    let (old, new) = fixtures();
    let old_equal = old.clone();
    let new_equal = new.clone();
    assert_eq!(old_equal, old);
    assert_eq!(new_equal, new);
    assert_real(
        old.eval(&[0.3, -0.7, 1.1]).unwrap(),
        new.eval(&[0.3, -0.7, 1.1]).unwrap(),
    );
    assert_eq!(old.to_string(), new.to_string());

    let mut group = c.benchmark_group("sym/surface/observable/sum");
    group.bench_function("new/clone", |b| b.iter(|| new.clone()));
    group.bench_function("old/clone", |b| b.iter(|| old.clone()));
    group.bench_function("new/equality", |b| b.iter(|| new == new_equal));
    group.bench_function("old/equality", |b| b.iter(|| old == old_equal));
    group.bench_function("new/display", |b| b.iter(|| new.to_string()));
    group.bench_function("old/display", |b| b.iter(|| old.to_string()));
    group.finish();
}
