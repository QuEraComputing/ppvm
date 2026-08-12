// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hint::black_box;

use approx::{AbsDiffEq, RelativeEq};
use criterion::Criterion;

use super::{NewKey, OldKey, assert_pair, build_new, build_old, terms};

pub fn bench(c: &mut Criterion) {
    approximate(c);
    formatting(c);
}

fn approximate(c: &mut Criterion) {
    const ABS_EPSILON: f64 = 1e-9;
    const REL_EPSILON: f64 = 1e-12;
    const MAX_RELATIVE: f64 = 1e-8;

    let data = terms(0, 192);
    let (old, new) = (build_old(&data), build_new(&data));
    let (old_equal, new_equal) = (old.clone(), new.clone());
    let (mut old_near, mut new_near) = (old.clone(), new.clone());
    let word = data[64].0.as_str();
    old_near += (OldKey::from(word), 1e-10);
    new_near += (NewKey::from(word), 1e-10);
    assert_pair(&old, &new);
    assert_pair(&old_near, &new_near);

    bench_abs(c, "equal", &old, &old_equal, &new, &new_equal, ABS_EPSILON);
    bench_abs(c, "near", &old, &old_near, &new, &new_near, ABS_EPSILON);
    bench_relative(
        c,
        "equal",
        &old,
        &old_equal,
        &new,
        &new_equal,
        REL_EPSILON,
        MAX_RELATIVE,
    );
    bench_relative(
        c,
        "near",
        &old,
        &old_near,
        &new,
        &new_near,
        REL_EPSILON,
        MAX_RELATIVE,
    );
}

fn bench_abs(
    c: &mut Criterion,
    fixture: &str,
    old: &super::OldSum,
    old_rhs: &super::OldSum,
    new: &super::NewSum,
    new_rhs: &super::NewSum,
    epsilon: f64,
) {
    let (old_result, new_result) = (
        old.abs_diff_eq(old_rhs, epsilon),
        new.abs_diff_eq(new_rhs, epsilon),
    );
    assert_eq!(old_result, new_result);
    assert!(old_result);
    let mut group = c.benchmark_group(format!("pauli_sum_surface/compare/abs_diff_eq_{fixture}"));
    group.bench_function("old", |b| {
        b.iter(|| black_box(old.abs_diff_eq(black_box(old_rhs), epsilon)))
    });
    group.bench_function("new", |b| {
        b.iter(|| black_box(new.abs_diff_eq(black_box(new_rhs), epsilon)))
    });
    group.finish();
}

#[allow(clippy::too_many_arguments)]
fn bench_relative(
    c: &mut Criterion,
    fixture: &str,
    old: &super::OldSum,
    old_rhs: &super::OldSum,
    new: &super::NewSum,
    new_rhs: &super::NewSum,
    epsilon: f64,
    max_relative: f64,
) {
    let (old_result, new_result) = (
        old.relative_eq(old_rhs, epsilon, max_relative),
        new.relative_eq(new_rhs, epsilon, max_relative),
    );
    assert_eq!(old_result, new_result);
    assert!(old_result);
    let mut group = c.benchmark_group(format!("pauli_sum_surface/compare/relative_eq_{fixture}"));
    group.bench_function("old", |b| {
        b.iter(|| black_box(old.relative_eq(black_box(old_rhs), epsilon, max_relative)))
    });
    group.bench_function("new", |b| {
        b.iter(|| black_box(new.relative_eq(black_box(new_rhs), epsilon, max_relative)))
    });
    group.finish();
}

fn formatting(c: &mut Criterion) {
    // One term per weight avoids unspecified hash iteration order within equal
    // weight classes while exercising the complete weight-sort path.
    let data: Vec<(String, f64)> = (0..=8)
        .map(|weight| {
            let word: String = (0..8)
                .map(|site| if site < weight { 'X' } else { 'I' })
                .collect();
            (word, weight as f64 + 0.125)
        })
        .collect();
    let (old, new) = (build_old(&data), build_new(&data));
    assert_pair(&old, &new);

    let old_display = old.to_string();
    let new_display = new.to_string();
    assert_eq!(old_display, new_display);
    let mut group = c.benchmark_group("pauli_sum_surface/format/display");
    group.bench_function("old", |b| b.iter(|| black_box(old.to_string())));
    group.bench_function("new", |b| b.iter(|| black_box(new.to_string())));
    group.finish();

    let old_debug = format!("{old:?}");
    let new_debug = format!("{new:?}");
    assert_eq!(old_debug, new_debug);
    let mut group = c.benchmark_group("pauli_sum_surface/format/debug");
    group.bench_function("old", |b| b.iter(|| black_box(format!("{old:?}"))));
    group.bench_function("new", |b| b.iter(|| black_box(format!("{new:?}"))));
    group.finish();
}
