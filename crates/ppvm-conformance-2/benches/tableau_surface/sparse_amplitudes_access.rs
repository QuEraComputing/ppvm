// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::BatchSize;
use num::complex::Complex64;

use super::sparse_amplitudes::{NewAmplitudes, OldAmplitudes, assert_equal};

fn digest_ref<'a>(iter: impl Iterator<Item = &'a (Complex64, u128)>) -> (usize, u128, f64) {
    iter.fold((0, 0, 0.0), |(n, index, norm), (value, i)| {
        (n + 1, index ^ *i, norm + value.norm_sqr())
    })
}

fn digest_owned(iter: impl IntoIterator<Item = (Complex64, u128)>) -> (usize, u128, f64) {
    iter.into_iter()
        .fold((0, 0, 0.0), |(n, index, norm), (value, i)| {
            (n + 1, index ^ i, norm + value.norm_sqr())
        })
}

pub fn bench(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    old: &OldAmplitudes,
    new: &NewAmplitudes,
) {
    let old_new: OldAmplitudes = ppvm_tableau::sparsevec::SparseVector::new();
    let new_new = NewAmplitudes::new();
    assert_equal(&old_new, &new_new);
    group.bench_function("new/old", |b| {
        b.iter(|| {
            std::hint::black_box(<OldAmplitudes as ppvm_tableau::sparsevec::SparseVector<
                _,
                _,
            >>::new())
        })
    });
    group.bench_function("new/new", |b| {
        b.iter(|| std::hint::black_box(NewAmplitudes::new()))
    });

    let old_default = OldAmplitudes::default();
    let new_default = NewAmplitudes::default();
    assert_equal(&old_default, &new_default);
    group.bench_function("default/old", |b| {
        b.iter(|| std::hint::black_box(OldAmplitudes::default()))
    });
    group.bench_function("default/new", |b| {
        b.iter(|| std::hint::black_box(NewAmplitudes::default()))
    });

    assert_equal(&old.clone(), &new.clone());
    group.bench_function("clone/old", |b| {
        b.iter(|| std::hint::black_box(old.clone()))
    });
    group.bench_function("clone/new", |b| {
        b.iter(|| std::hint::black_box(new.clone()))
    });

    let old_clone = old.clone();
    let new_clone = new.clone();
    assert!(old == &old_clone);
    assert!(new == &new_clone);
    group.bench_function("equality/old", |b| {
        b.iter(|| std::hint::black_box(old == &old_clone))
    });
    group.bench_function("equality/new", |b| {
        b.iter(|| std::hint::black_box(new == &new_clone))
    });

    assert_eq!(ppvm_tableau::sparsevec::SparseVector::len(old), new.len());
    group.bench_function("len/old", |b| {
        b.iter(|| std::hint::black_box(ppvm_tableau::sparsevec::SparseVector::len(old)))
    });
    group.bench_function("len/new", |b| b.iter(|| std::hint::black_box(new.len())));
    assert_eq!(
        ppvm_tableau::sparsevec::SparseVector::is_empty(old),
        new.is_empty()
    );
    group.bench_function("is_empty/old", |b| {
        b.iter(|| std::hint::black_box(ppvm_tableau::sparsevec::SparseVector::is_empty(old)))
    });
    group.bench_function("is_empty/new", |b| {
        b.iter(|| std::hint::black_box(new.is_empty()))
    });

    let expected = digest_ref(old.as_slice().iter());
    assert_eq!(expected, digest_ref(new.entries().iter()));
    group.bench_function("entries-traversal/old", |b| {
        b.iter(|| std::hint::black_box(digest_ref(old.as_slice().iter())))
    });
    group.bench_function("entries-traversal/new", |b| {
        b.iter(|| std::hint::black_box(digest_ref(new.entries().iter())))
    });
    assert_eq!(
        digest_ref(ppvm_tableau::sparsevec::SparseVector::iter(old)),
        digest_ref(new.iter())
    );
    group.bench_function("iter-traversal/old", |b| {
        b.iter(|| {
            std::hint::black_box(digest_ref(ppvm_tableau::sparsevec::SparseVector::iter(old)))
        })
    });
    group.bench_function("iter-traversal/new", |b| {
        b.iter(|| std::hint::black_box(digest_ref(new.iter())))
    });

    assert_eq!(digest_owned(old.clone()), digest_owned(new.clone()));
    group.bench_function("into_iter/old", |b| {
        b.iter_batched(
            || old.clone(),
            |values| std::hint::black_box(digest_owned(values)),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("into_iter/new", |b| {
        b.iter_batched(
            || new.clone(),
            |values| std::hint::black_box(digest_owned(values)),
            BatchSize::SmallInput,
        )
    });
}
