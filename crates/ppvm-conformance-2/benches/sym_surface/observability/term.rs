// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_conformance_2::sym::{NewTerm, OldTerm};

use super::super::assert_real;

fn fixtures() -> (OldTerm, NewTerm) {
    let old = OldTerm::from(2.0)
        + OldTerm::var(0).sin() * 1.5
        + OldTerm::var(1).cos() * -0.5
        + OldTerm::var(2).sin() * OldTerm::var(2).cos() * 0.25;
    let new = NewTerm::from(2.0)
        + NewTerm::var(0).sin() * 1.5
        + NewTerm::var(1).cos() * -0.5
        + NewTerm::var(2).sin() * NewTerm::var(2).cos() * 0.25;
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

    let mut group = c.benchmark_group("sym/surface/observable/term");
    group.bench_function("new/clone", |b| b.iter(|| new.clone()));
    group.bench_function("old/clone", |b| b.iter(|| old.clone()));
    group.bench_function("new/equality", |b| b.iter(|| new == new_equal));
    group.bench_function("old/equality", |b| b.iter(|| old == old_equal));
    group.bench_function("new/display", |b| b.iter(|| new.to_string()));
    group.bench_function("old/display", |b| b.iter(|| old.to_string()));
    group.finish();
}
