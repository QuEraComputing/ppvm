// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use num::complex::Complex64;

use super::bench_mut_pair;

pub(super) type OldAmplitudes = Vec<(Complex64, u128)>;
pub(super) type NewAmplitudes = ppvm_tableau_2::Amplitudes<u128>;

pub(super) fn pair() -> (OldAmplitudes, NewAmplitudes) {
    let mut old = Vec::with_capacity(256);
    let mut new = NewAmplitudes::new();
    new.reserve(256);
    for i in 0..128u128 {
        let value = Complex64::new(1.0 + i as f64 / 128.0, (i % 7) as f64 / 16.0);
        ppvm_tableau::sparsevec::SparseVector::unsafe_insert(&mut old, i, value);
        new.unsafe_insert(i, value);
    }
    assert_eq!(old.as_slice(), new.entries());
    (old, new)
}

pub(super) fn assert_equal(old: &OldAmplitudes, new: &NewAmplitudes) {
    assert_eq!(old.as_slice(), new.entries());
}

macro_rules! mutation {
    ($group:expr, $name:expr, $old:expr, $new:expr, $old_op:expr, $new_op:expr) => {{
        let (mut old_check, mut new_check) = (($old).clone(), ($new).clone());
        ($old_op)(&mut old_check);
        ($new_op)(&mut new_check);
        assert_equal(&old_check, &new_check);
        bench_mut_pair!($group, $name, $old, $new, $old_op, $new_op);
    }};
}

pub fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("tableau-surface/sparse-amplitudes");
    let (old, new) = pair();
    let factor = Complex64::new(0.75, -0.25);

    mutation!(
        group,
        "unsafe_insert",
        old,
        new,
        |v: &mut OldAmplitudes| {
            ppvm_tableau::sparsevec::SparseVector::unsafe_insert(v, 200, Complex64::new(0.5, 0.25))
        },
        |v: &mut NewAmplitudes| v.unsafe_insert(200, Complex64::new(0.5, 0.25))
    );
    mutation!(
        group,
        "add_or_insert/hit",
        old,
        new,
        |v: &mut OldAmplitudes| {
            ppvm_tableau::sparsevec::SparseVector::add_or_insert(v, 64, Complex64::new(0.5, 0.25))
        },
        |v: &mut NewAmplitudes| v.add_or_insert(64, Complex64::new(0.5, 0.25))
    );
    mutation!(
        group,
        "add_or_insert/miss",
        old,
        new,
        |v: &mut OldAmplitudes| {
            ppvm_tableau::sparsevec::SparseVector::add_or_insert(v, 200, Complex64::new(0.5, 0.25))
        },
        |v: &mut NewAmplitudes| v.add_or_insert(200, Complex64::new(0.5, 0.25))
    );
    mutation!(
        group,
        "mul_by",
        old,
        new,
        |v: &mut OldAmplitudes| { ppvm_tableau::sparsevec::SparseVector::mul_by(v, factor) },
        |v: &mut NewAmplitudes| v.mul_by(factor)
    );
    mutation!(
        group,
        "mul_element_by/hit",
        old,
        new,
        |v: &mut OldAmplitudes| {
            ppvm_tableau::sparsevec::SparseVector::mul_element_by(v, 64, factor)
        },
        |v: &mut NewAmplitudes| v.mul_element_by(64, factor)
    );
    mutation!(
        group,
        "mul_element_by/miss",
        old,
        new,
        |v: &mut OldAmplitudes| {
            ppvm_tableau::sparsevec::SparseVector::mul_element_by(v, 200, factor)
        },
        |v: &mut NewAmplitudes| v.mul_element_by(200, factor)
    );
    mutation!(
        group,
        "trim",
        old,
        new,
        |v: &mut OldAmplitudes| {
            ppvm_tableau::sparsevec::SparseVector::trim(v, Complex64::new(1.5, 0.0))
        },
        |v: &mut NewAmplitudes| v.trim(Complex64::new(1.5, 0.0))
    );
    mutation!(
        group,
        "retain",
        old,
        new,
        |v: &mut OldAmplitudes| {
            ppvm_tableau::sparsevec::SparseVector::retain(v, |(_, i)| *i % 3 != 0)
        },
        |v: &mut NewAmplitudes| v.retain_entries(|(_, i)| *i % 3 != 0)
    );
    mutation!(
        group,
        "reserve",
        old,
        new,
        |v: &mut OldAmplitudes| { ppvm_tableau::sparsevec::SparseVector::reserve(v, 512) },
        |v: &mut NewAmplitudes| v.reserve(512)
    );
    mutation!(
        group,
        "normalize",
        old,
        new,
        |v: &mut OldAmplitudes| ppvm_tableau::sparsevec::SparseVector::normalize(v),
        |v: &mut NewAmplitudes| v.normalize()
    );

    let old_hit = ppvm_tableau::sparsevec::SparseVector::get(&old, &64);
    let new_hit = new.get(&64);
    assert_eq!(old_hit, new_hit);
    group.bench_function("get/hit/old", |b| {
        b.iter(|| std::hint::black_box(ppvm_tableau::sparsevec::SparseVector::get(&old, &64)))
    });
    group.bench_function("get/hit/new", |b| {
        b.iter(|| std::hint::black_box(new.get(&64)))
    });
    let old_miss = ppvm_tableau::sparsevec::SparseVector::get(&old, &200);
    let new_miss = new.get(&200);
    assert_eq!(old_miss, new_miss);
    group.bench_function("get/miss/old", |b| {
        b.iter(|| std::hint::black_box(ppvm_tableau::sparsevec::SparseVector::get(&old, &200)))
    });
    group.bench_function("get/miss/new", |b| {
        b.iter(|| std::hint::black_box(new.get(&200)))
    });
    super::sparse_amplitudes_access::bench(&mut group, &old, &new);
    group.finish();
}
