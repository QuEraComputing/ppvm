// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Comparable old/new pattern operations. Pattern and word setup is untimed;
//! enumeration deliberately includes iterator construction and consumption.

use criterion::{BenchmarkId, Criterion};
use ppvm_pauli_sum_2::PauliPattern as NewPattern;
use ppvm_pauli_word::pattern::{Contains, PauliPattern as OldPattern};
use std::hint::black_box;

use super::common::*;

pub fn bench(c: &mut Criterion) {
    parse(c);
    matching(c);
    bounded_enumeration(c);
}

fn parse(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/pattern/parse");
    for (label, source) in [
        ("indexed", "X1Y17Z255"),
        ("optional_repeat", "[XYZ]?{64}"),
        ("star", "[XY]*Z255"),
    ] {
        let old = OldPattern::parse(source).unwrap();
        let new = NewPattern::parse(source).unwrap();
        assert_eq!(old.to_string(), new.to_string(), "{label}");
        g.bench_with_input(BenchmarkId::new("new", label), &source, |b, source| {
            b.iter(|| black_box(NewPattern::parse(black_box(*source)).unwrap()))
        });
        g.bench_with_input(BenchmarkId::new("old", label), &source, |b, source| {
            b.iter(|| black_box(OldPattern::parse(black_box(*source)).unwrap()))
        });
    }
    g.finish();
}

fn matching(c: &mut Criterion) {
    const SOURCE: &str = "[XYZ]?{256}";
    let old_pattern = OldPattern::parse(SOURCE).unwrap();
    let new_pattern = NewPattern::parse(SOURCE).unwrap();
    let text = ordinary_string(WIDTH);
    let old_word = OldWord::from(text.as_str());
    let new_word = NewWord::from(text.as_str());
    let old_lossy = OldLossy::from(text.as_str());
    let new_lossy = NewLossy::from(text.as_str());
    assert_eq!(
        old_pattern.contains(&old_word),
        new_pattern.matches(&new_word)
    );
    assert_eq!(
        old_pattern.contains(&old_lossy),
        new_pattern.matches(&new_lossy)
    );

    let mut g = c.benchmark_group("word_surface/pattern/match_contains/256");
    g.bench_function("new/ordinary", |b| {
        b.iter(|| black_box(new_pattern.matches(black_box(&new_word))))
    });
    g.bench_function("old/ordinary", |b| {
        b.iter(|| black_box(old_pattern.contains(black_box(&old_word))))
    });
    g.bench_function("new/lossy_present", |b| {
        b.iter(|| black_box(new_pattern.matches(black_box(&new_lossy))))
    });
    g.bench_function("old/lossy_present", |b| {
        b.iter(|| black_box(old_pattern.contains(black_box(&old_lossy))))
    });
    g.finish();
}

fn bounded_enumeration(c: &mut Criterion) {
    const SOURCE: &str = "[XYZ]{4}";
    const WIDTH: usize = 8;
    let old_pattern = OldPattern::parse(SOURCE).unwrap();
    let new_pattern = NewPattern::parse(SOURCE).unwrap();
    let old_output: Vec<_> = old_pattern
        .enumerate_matches::<Storage>(WIDTH)
        .map(|word| word.to_string())
        .collect();
    let new_output: Vec<_> = new_pattern
        .enumerate_matches::<Storage>(WIDTH)
        .map(|word| word.to_string())
        .collect();
    assert_eq!(old_output, new_output);

    let mut g = c.benchmark_group("word_surface/pattern/bounded_enumeration/8");
    g.bench_function("new/enumerate_all", |b| {
        b.iter(|| {
            black_box(&new_pattern)
                .enumerate_matches::<Storage>(black_box(WIDTH))
                .fold(0usize, |count, word| {
                    black_box(word);
                    count + 1
                })
        })
    });
    g.bench_function("old/enumerate_all", |b| {
        b.iter(|| {
            black_box(&old_pattern)
                .enumerate_matches::<Storage>(black_box(WIDTH))
                .fold(0usize, |count, word| {
                    black_box(word);
                    count + 1
                })
        })
    });
    g.finish();
}
