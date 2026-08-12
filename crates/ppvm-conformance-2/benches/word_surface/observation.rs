// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Equality and `Display` costs over prebuilt, matched values.

use criterion::Criterion;
use ppvm_pauli_sum_2::PauliPattern as NewPattern;
use ppvm_pauli_word::pattern::PauliPattern as OldPattern;
use std::hint::black_box;

use super::common::*;

pub fn bench(c: &mut Criterion) {
    ordinary(c);
    lossy(c);
    phased(c);
    pattern(c);
}

fn ordinary(c: &mut Criterion) {
    let text = ordinary_string(WIDTH);
    let (new_a, new_b) = (NewWord::from(text.as_str()), NewWord::from(text.as_str()));
    let (old_a, old_b) = (OldWord::from(text.as_str()), OldWord::from(text.as_str()));
    assert_eq!(new_a == new_b, old_a == old_b);
    assert_eq!(new_a.to_string(), old_a.to_string());
    let mut g = c.benchmark_group("word_surface/ordinary/observation/256");
    g.bench_function("new/equality", |b| {
        b.iter(|| black_box(black_box(&new_a) == black_box(&new_b)))
    });
    g.bench_function("old/equality", |b| {
        b.iter(|| black_box(black_box(&old_a) == black_box(&old_b)))
    });
    g.bench_function("new/display", |b| {
        b.iter(|| black_box(black_box(&new_a).to_string()))
    });
    g.bench_function("old/display", |b| {
        b.iter(|| black_box(black_box(&old_a).to_string()))
    });
    g.finish();
}

fn lossy(c: &mut Criterion) {
    let text = lossy_string(WIDTH);
    let (new_a, new_b) = (NewLossy::from(text.as_str()), NewLossy::from(text.as_str()));
    let (old_a, old_b) = (OldLossy::from(text.as_str()), OldLossy::from(text.as_str()));
    assert_eq!(new_a == new_b, old_a == old_b);
    assert_eq!(new_a.to_string(), old_a.to_string());
    let mut g = c.benchmark_group("word_surface/lossy/observation/256");
    g.bench_function("new/equality", |b| {
        b.iter(|| black_box(black_box(&new_a) == black_box(&new_b)))
    });
    g.bench_function("old/equality", |b| {
        b.iter(|| black_box(black_box(&old_a) == black_box(&old_b)))
    });
    g.bench_function("new/display", |b| {
        b.iter(|| black_box(black_box(&new_a).to_string()))
    });
    g.bench_function("old/display", |b| {
        b.iter(|| black_box(black_box(&old_a).to_string()))
    });
    g.finish();
}

fn phased(c: &mut Criterion) {
    let text = phased_string(WIDTH);
    let (new_a, new_b) = (
        NewPhased::from(text.as_str()),
        NewPhased::from(text.as_str()),
    );
    let (old_a, old_b) = (
        OldPhased::from(text.as_str()),
        OldPhased::from(text.as_str()),
    );
    assert_eq!(new_a == new_b, old_a == old_b);
    assert_eq!(new_a.to_string(), old_a.to_string());
    let mut g = c.benchmark_group("word_surface/phased/observation/256");
    g.bench_function("new/equality", |b| {
        b.iter(|| black_box(black_box(&new_a) == black_box(&new_b)))
    });
    g.bench_function("old/equality", |b| {
        b.iter(|| black_box(black_box(&old_a) == black_box(&old_b)))
    });
    g.bench_function("new/display", |b| {
        b.iter(|| black_box(black_box(&new_a).to_string()))
    });
    g.bench_function("old/display", |b| {
        b.iter(|| black_box(black_box(&old_a).to_string()))
    });
    g.finish();
}

fn pattern(c: &mut Criterion) {
    const SOURCE: &str = "[XYZ]?{64}X127[XY]*";
    let (new_a, new_b) = (
        NewPattern::parse(SOURCE).unwrap(),
        NewPattern::parse(SOURCE).unwrap(),
    );
    let (old_a, old_b) = (
        OldPattern::parse(SOURCE).unwrap(),
        OldPattern::parse(SOURCE).unwrap(),
    );
    assert_eq!(new_a == new_b, old_a == old_b);
    assert_eq!(new_a.to_string(), old_a.to_string());
    let mut g = c.benchmark_group("word_surface/pattern/observation");
    g.bench_function("new/equality", |b| {
        b.iter(|| black_box(black_box(&new_a) == black_box(&new_b)))
    });
    g.bench_function("old/equality", |b| {
        b.iter(|| black_box(black_box(&old_a) == black_box(&old_b)))
    });
    g.bench_function("new/display", |b| {
        b.iter(|| black_box(black_box(&new_a).to_string()))
    });
    g.bench_function("old/display", |b| {
        b.iter(|| black_box(black_box(&old_a).to_string()))
    });
    g.finish();
}
