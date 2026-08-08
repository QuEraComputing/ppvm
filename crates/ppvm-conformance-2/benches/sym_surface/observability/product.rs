// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::{Hash, Hasher};

use criterion::{BatchSize, Criterion};
use ppvm_sym::Prod as OldProd;
use ppvm_sym_2::Prod as NewProd;

fn fixtures() -> (OldProd, NewProd) {
    let mut old = OldProd::sin(1);
    let mut new = NewProd::sin(1);
    for var in [1, 3, 5] {
        old.mul_sin(var);
        old.mul_cos(var);
        new.mul_sin(var);
        new.mul_cos(var);
    }
    (old, new)
}

fn hash<T: Hash>(value: &T) -> u64 {
    let mut hasher = fxhash::FxHasher64::default();
    value.hash(&mut hasher);
    hasher.finish()
}

pub(super) fn bench(c: &mut Criterion) {
    let (old, new) = fixtures();
    assert_eq!(old.pow(), new.pow());
    assert_eq!(old.pow(), 7);
    assert_eq!(old.sin_pow(), new.sin_pow());
    assert_eq!(old.sin_pow(), 4);
    assert_eq!(old.cos_pow(), new.cos_pow());
    assert_eq!(old.cos_pow(), 3);
    assert_eq!(old.clone(), old);
    assert_eq!(new.clone(), new);
    assert_eq!(hash(&old), hash(&old.clone()));
    assert_eq!(hash(&new), hash(&new.clone()));
    assert_eq!(old.to_string(), new.to_string());

    let old_equal = old.clone();
    let new_equal = new.clone();
    let mut group = c.benchmark_group("sym/surface/observable/product");
    group.bench_function("new/pow", |b| b.iter(|| new.pow()));
    group.bench_function("old/pow", |b| b.iter(|| old.pow()));
    group.bench_function("new/sin_pow", |b| b.iter(|| new.sin_pow()));
    group.bench_function("old/sin_pow", |b| b.iter(|| old.sin_pow()));
    group.bench_function("new/cos_pow", |b| b.iter(|| new.cos_pow()));
    group.bench_function("old/cos_pow", |b| b.iter(|| old.cos_pow()));
    group.bench_function("new/clone", |b| b.iter(|| new.clone()));
    group.bench_function("old/clone", |b| b.iter(|| old.clone()));
    group.bench_function("new/equality", |b| b.iter(|| new == new_equal));
    group.bench_function("old/equality", |b| b.iter(|| old == old_equal));
    group.bench_function("new/hash", |b| {
        b.iter_batched(
            fxhash::FxHasher64::default,
            |mut hasher| {
                new.hash(&mut hasher);
                hasher.finish()
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old/hash", |b| {
        b.iter_batched(
            fxhash::FxHasher64::default,
            |mut hasher| {
                old.hash(&mut hasher);
                hasher.finish()
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("new/display", |b| b.iter(|| new.to_string()));
    group.bench_function("old/display", |b| b.iter(|| old.to_string()));
    group.finish();
}
