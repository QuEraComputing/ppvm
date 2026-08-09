// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Explicit phased construction and wrapper decomposition/access.

use criterion::{BatchSize, Criterion};
use std::hint::black_box;

use super::common::*;

pub fn bench(c: &mut Criterion) {
    explicit_phase(c);
    access_and_decompose(c);
}

fn explicit_phase(c: &mut Criterion) {
    let text = ordinary_string(WIDTH);
    let new_word = NewWord::from(text.as_str());
    let old_word = OldWord::from(text.as_str());
    let new = NewPhased::with_phase(new_word, ppvm_traits_2::Phase::NegI);
    let old = OldPhased::build_from_word(old_word, 3);
    assert_eq!(new.to_string(), old.to_string());

    let mut g = c.benchmark_group("word_surface/construct/phased_explicit_phase/256");
    g.bench_function("new/from_word_and_phase", |b| {
        b.iter_batched(
            || new_word,
            |word| {
                black_box(NewPhased::with_phase(
                    black_box(word),
                    black_box(ppvm_traits_2::Phase::NegI),
                ))
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/from_word_and_phase", |b| {
        b.iter_batched(
            || old_word,
            |word| black_box(OldPhased::build_from_word(black_box(word), black_box(3))),
            BatchSize::SmallInput,
        )
    });
    g.finish();
}

fn access_and_decompose(c: &mut Criterion) {
    let text = phased_string(WIDTH);
    let new = NewPhased::from(text.as_str());
    let old = OldPhased::from(text.as_str());
    assert_eq!(new.word().to_string(), old.word.to_string());
    let (new_word, new_phase) = new.clone().into_parts();
    let (old_word, old_phase) = (old.word, old.phase);
    assert_eq!(new_word.to_string(), old_word.to_string());
    assert_eq!(new_phase.exponent(), old_phase);

    let mut g = c.benchmark_group("word_surface/phased/wrapper/256");
    g.bench_function("new/inner_word", |b| {
        b.iter(|| black_box(black_box(&new).word()))
    });
    g.bench_function("old/inner_word", |b| {
        b.iter(|| black_box(&black_box(&old).word))
    });
    g.bench_function("new/into_parts", |b| {
        b.iter_batched(
            || new.clone(),
            |word| black_box(black_box(word).into_parts()),
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/into_parts", |b| {
        b.iter_batched(
            || old,
            |word| black_box((black_box(word).word, black_box(word).phase)),
            BatchSize::SmallInput,
        )
    });
    g.finish();
}
