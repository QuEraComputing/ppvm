// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, BenchmarkId, Criterion};
use ppvm_traits::traits::PauliWordTrait;
use ppvm_traits_2::{Indexable, Word};
use std::hint::black_box;

use super::common::*;

pub fn bench(c: &mut Criterion) {
    parse_and_new(c);
    phased_from_word(c);
    clone_copy(c);
}

fn phased_from_word(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/construct/phased_from_word/256");
    let s = ordinary_string(WIDTH);
    let new_word = NewWord::from(s.as_str());
    let old_word = OldWord::from(s.as_str());
    let new = NewPhased::new(new_word);
    let old = OldPhased::build_from_word(old_word, 0);
    assert_eq!(new.to_string(), old.to_string());
    g.bench_function("new/from_existing_word", |b| {
        b.iter_batched(
            || new_word,
            |word| black_box(NewPhased::new(black_box(word))),
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/from_existing_word", |b| {
        b.iter_batched(
            || old_word,
            |word| black_box(OldPhased::build_from_word(black_box(word), 0)),
            BatchSize::SmallInput,
        )
    });
    g.finish();
}

fn parse_and_new(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/construct");
    for width in [8, 64, WIDTH] {
        let plain = ordinary_string(width);
        let lossy = lossy_string(width);
        let phased = phased_string(width);

        assert_eq!(NewWord::from(plain.as_str()).n_sites(), width);
        assert_eq!(OldWord::from(plain.as_str()).n_qubits(), width);
        g.bench_with_input(
            BenchmarkId::new("ordinary/new/parse", width),
            &plain,
            |b, s| b.iter(|| black_box(NewWord::from(black_box(s.as_str())))),
        );
        g.bench_with_input(
            BenchmarkId::new("ordinary/old/parse", width),
            &plain,
            |b, s| b.iter(|| black_box(OldWord::from(black_box(s.as_str())))),
        );

        assert_eq!(NewLossy::from(lossy.as_str()).n_sites(), width);
        assert_eq!(OldLossy::from(lossy.as_str()).n_qubits(), width);
        g.bench_with_input(
            BenchmarkId::new("lossy/new/parse", width),
            &lossy,
            |b, s| b.iter(|| black_box(NewLossy::from(black_box(s.as_str())))),
        );
        g.bench_with_input(
            BenchmarkId::new("lossy/old/parse", width),
            &lossy,
            |b, s| b.iter(|| black_box(OldLossy::from(black_box(s.as_str())))),
        );

        assert_eq!(NewPhased::from(phased.as_str()).n_sites(), width);
        assert_eq!(OldPhased::from(phased.as_str()).n_qubits(), width);
        g.bench_with_input(
            BenchmarkId::new("phased/new/parse", width),
            &phased,
            |b, s| b.iter(|| black_box(NewPhased::from(black_box(s.as_str())))),
        );
        g.bench_with_input(
            BenchmarkId::new("phased/old/parse", width),
            &phased,
            |b, s| b.iter(|| black_box(OldPhased::from(black_box(s.as_str())))),
        );

        assert_eq!(NewWord::new(width).n_sites(), width);
        assert_eq!(<OldWord as PauliWordTrait>::new(width).n_qubits(), width);
        assert_eq!(NewLossy::new(width).n_sites(), width);
        assert_eq!(<OldLossy as PauliWordTrait>::new(width).n_qubits(), width);
        assert_eq!(NewPhased::new(NewWord::new(width)).n_sites(), width);
        assert_eq!(OldPhased::new(width).n_qubits(), width);
        g.bench_function(BenchmarkId::new("ordinary/new/new_identity", width), |b| {
            b.iter(|| black_box(NewWord::new(black_box(width))))
        });
        g.bench_function(BenchmarkId::new("ordinary/old/new_identity", width), |b| {
            b.iter(|| black_box(<OldWord as PauliWordTrait>::new(black_box(width))))
        });
        g.bench_function(BenchmarkId::new("lossy/new/new_identity", width), |b| {
            b.iter(|| black_box(NewLossy::new(black_box(width))))
        });
        g.bench_function(BenchmarkId::new("lossy/old/new_identity", width), |b| {
            b.iter(|| black_box(<OldLossy as PauliWordTrait>::new(black_box(width))))
        });
        g.bench_function(BenchmarkId::new("phased/new/new_identity", width), |b| {
            b.iter(|| black_box(NewPhased::new(NewWord::new(black_box(width)))))
        });
        g.bench_function(BenchmarkId::new("phased/old/new_identity", width), |b| {
            b.iter(|| black_box(OldPhased::new(black_box(width))))
        });
    }
    g.finish();
}

fn clone_copy(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/clone_copy/256");
    let plain = ordinary_string(WIDTH);
    let lossy = lossy_string(WIDTH);
    let phased = phased_string(WIDTH);

    let ordinary_cold = NewWord::from(plain.as_str());
    let ordinary_warm = NewWord::from(plain.as_str());
    black_box(ordinary_warm.key_hash());
    assert_eq!(ordinary_cold.clone(), ordinary_cold);
    assert_eq!(ordinary_warm.clone(), ordinary_warm);
    g.bench_function("ordinary/new/clone_cold", |b| {
        b.iter(|| black_box(*black_box(&ordinary_cold)))
    });
    g.bench_function("ordinary/new/clone_warm", |b| {
        b.iter(|| black_box(*black_box(&ordinary_warm)))
    });
    let ordinary_old = OldWord::from(plain.as_str());
    assert_eq!(*black_box(&ordinary_old), ordinary_old);
    // Old is eagerly hashed, so its copy-equivalent cost is identical in both
    // paired cache-state rows.
    g.bench_function("ordinary/old/clone_cold", |b| {
        b.iter(|| black_box(*black_box(&ordinary_old)))
    });
    g.bench_function("ordinary/old/clone_warm", |b| {
        b.iter(|| black_box(*black_box(&ordinary_old)))
    });

    let lossy_cold = NewLossy::from(lossy.as_str());
    let lossy_warm = NewLossy::from(lossy.as_str());
    black_box(lossy_warm.key_hash());
    assert_eq!(lossy_cold.clone(), lossy_cold);
    assert_eq!(lossy_warm.clone(), lossy_warm);
    g.bench_function("lossy/new/clone_cold", |b| {
        b.iter(|| black_box(black_box(&lossy_cold).clone()))
    });
    g.bench_function("lossy/new/clone_warm", |b| {
        b.iter(|| black_box(black_box(&lossy_warm).clone()))
    });
    let lossy_old = OldLossy::from(lossy.as_str());
    assert_eq!(lossy_old.clone(), lossy_old);
    g.bench_function("lossy/old/clone_cold", |b| {
        b.iter(|| black_box(black_box(&lossy_old).clone()))
    });
    g.bench_function("lossy/old/clone_warm", |b| {
        b.iter(|| black_box(black_box(&lossy_old).clone()))
    });

    let phased_cold = NewPhased::from(phased.as_str());
    let phased_warm = NewPhased::from(phased.as_str());
    black_box(phased_warm.word().key_hash());
    assert_eq!(phased_cold.clone(), phased_cold);
    assert_eq!(phased_warm.clone(), phased_warm);
    g.bench_function("phased/new/clone_cold", |b| {
        b.iter(|| black_box(black_box(&phased_cold).clone()))
    });
    g.bench_function("phased/new/clone_warm", |b| {
        b.iter(|| black_box(black_box(&phased_warm).clone()))
    });
    let phased_old = OldPhased::from(phased.as_str());
    assert_eq!(*black_box(&phased_old), phased_old);
    g.bench_function("phased/old/clone_cold", |b| {
        b.iter(|| black_box(*black_box(&phased_old)))
    });
    g.bench_function("phased/old/clone_warm", |b| {
        b.iter(|| black_box(*black_box(&phased_old)))
    });
    g.finish();
}
