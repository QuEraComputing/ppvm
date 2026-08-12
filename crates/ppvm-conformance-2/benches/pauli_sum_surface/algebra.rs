// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hint::black_box;

use criterion::{BatchSize, Criterion};
use num::Complex;
use ppvm_pauli_sum::config::fxhash::Byte;
use ppvm_pauli_sum::sum::PauliSum;
use ppvm_pauli_sum_2::{HashMapStore, NoPolicy, Sum};
use ppvm_traits::traits::NoStrategy;

use super::{CAPACITY, N, NewKey, OldKey, assert_pair, bench_mut, build_new, build_old, terms};

type OldComplex = PauliSum<Byte<8, Complex<f64>, NoStrategy, OldKey>>;
type NewComplex = Sum<HashMapStore<NewKey, Complex<f64>>, NoPolicy>;

pub fn bench(c: &mut Criterion) {
    bench_mut(
        c,
        "algebra/scale",
        |s| *s *= 1.000_001,
        |s| s.scale(&1.000_001),
    );
    overlap(c);
    mul_word(c);
}

fn overlap(c: &mut Criterion) {
    let a = terms(0, 192);
    let b = terms(96, 192);
    let (old_a, old_b) = (build_old(&a), build_old(&b));
    let (new_a, new_b) = (build_new(&a), build_new(&b));
    assert_pair(&old_a, &new_a);
    assert_pair(&old_b, &new_b);
    let old_value = old_a.overlap(&old_b);
    let new_value = new_a.overlap(&new_b);
    assert!((old_value - new_value).abs() <= 1e-9_f64.max(old_value.abs() * 1e-10));

    let mut group = c.benchmark_group("pauli_sum_surface/algebra/overlap");
    group.bench_function("old", |b| {
        b.iter(|| black_box(old_a.overlap(black_box(&old_b))))
    });
    group.bench_function("new", |b| {
        b.iter(|| black_box(new_a.overlap(black_box(&new_b))))
    });
    group.finish();
}

fn complex_terms() -> Vec<(String, Complex<f64>)> {
    terms(0, 192)
        .into_iter()
        .enumerate()
        .map(|(i, (word, coeff))| (word, Complex::new(coeff, (i % 13) as f64 / 17.0 - 0.3)))
        .collect()
}

fn old_complex(data: &[(String, Complex<f64>)]) -> OldComplex {
    let mut sum = OldComplex::builder()
        .n_qubits(N)
        .strategy(NoStrategy)
        .capacity(CAPACITY)
        .build();
    for (word, coeff) in data {
        sum += (word.as_str(), *coeff);
    }
    sum
}

fn new_complex(data: &[(String, Complex<f64>)]) -> NewComplex {
    let mut sum = NewComplex::with_capacity(N, NoPolicy, CAPACITY);
    for (word, coeff) in data {
        sum += (NewKey::from(word.as_str()), *coeff);
    }
    sum
}

fn complex_support_old(sum: &OldComplex) -> Vec<(String, Complex<f64>)> {
    let mut out: Vec<_> = sum.iter().map(|(w, c)| (w.to_string(), *c)).collect();
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

fn complex_support_new(sum: &NewComplex) -> Vec<(String, Complex<f64>)> {
    let mut out: Vec<_> = sum.iter().map(|(w, c)| (w.to_string(), c)).collect();
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

fn mul_word(c: &mut Criterion) {
    let data = complex_terms();
    let (old, new) = (old_complex(&data), new_complex(&data));
    let (old_word, new_word) = (OldKey::from("XYZXYZXY"), NewKey::from("XYZXYZXY"));
    let (mut op, mut np) = (old.clone(), new.clone());
    op *= old_word;
    np.mul_word_assign(&new_word);
    let old_support = complex_support_old(&op);
    let new_support = complex_support_new(&np);
    assert_eq!(old_support.len(), new_support.len());
    for ((ok, oc), (nk, nc)) in old_support.iter().zip(&new_support) {
        assert_eq!(ok, nk);
        assert!((*oc - *nc).norm() <= 1e-12);
    }

    let mut group = c.benchmark_group("pauli_sum_surface/algebra/mul_word");
    group.bench_function("old", |b| {
        b.iter_batched_ref(
            || old.clone(),
            |s| *s *= black_box(old_word),
            BatchSize::LargeInput,
        )
    });
    group.bench_function("new", |b| {
        b.iter_batched_ref(
            || new.clone(),
            |s| s.mul_word_assign(black_box(&new_word)),
            BatchSize::LargeInput,
        )
    });
    group.finish();
}
