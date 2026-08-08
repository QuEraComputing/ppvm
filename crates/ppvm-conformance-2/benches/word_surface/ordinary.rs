// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_traits::traits::{PauliIter, PauliWordTrait};
use ppvm_traits_2::{IdentityBuildHasher, Indexable, KeyProduct, PauliBits, Word};
use std::{hash::BuildHasher, hint::black_box};

use super::common::*;

pub fn bench(c: &mut Criterion) {
    reads(c);
    hash(c);
    product(c);
}

fn reads(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/ordinary/read/256");
    let s = ordinary_string(WIDTH);
    let new = NewWord::from(s.as_str());
    let old = OldWord::from(s.as_str());
    assert_eq!(new.weight(), old.weight());
    assert_eq!(
        PauliBits::loss_weight(&new),
        PauliWordTrait::loss_weight(&old)
    );
    assert_eq!(old_pauli(new.get(SITE)), old.get(SITE));
    assert_eq!(new.x_bit(SITE), old.get_xbit(SITE));
    assert_eq!(new.z_bit(SITE), old.get_zbit(SITE));
    let ni: String = new.iter().map(|p| format!("{p:?}")).collect();
    let oi: String = PauliIter::iter(&old).map(|p| p.to_string()).collect();
    assert_eq!(ni, oi);
    g.bench_function("new/width", |b| {
        b.iter(|| black_box(black_box(&new).n_sites()))
    });
    g.bench_function("old/width", |b| {
        b.iter(|| black_box(black_box(&old).n_qubits()))
    });
    g.bench_function("new/iter_traverse", |b| {
        b.iter(|| {
            black_box(&new).iter().for_each(|p| {
                black_box(p);
            })
        })
    });
    g.bench_function("old/iter_traverse", |b| {
        b.iter(|| {
            PauliIter::iter(black_box(&old)).for_each(|p| {
                black_box(p);
            })
        })
    });
    g.bench_function("new/weight", |b| {
        b.iter(|| black_box(black_box(&new).weight()))
    });
    g.bench_function("old/weight", |b| {
        b.iter(|| black_box(black_box(&old).weight()))
    });
    g.bench_function("new/loss_weight", |b| {
        b.iter(|| black_box(PauliBits::loss_weight(black_box(&new))))
    });
    g.bench_function("old/loss_weight", |b| {
        b.iter(|| black_box(PauliWordTrait::loss_weight(black_box(&old))))
    });
    g.bench_function("new/get", |b| {
        b.iter(|| black_box(black_box(&new).get(black_box(SITE))))
    });
    g.bench_function("old/get", |b| {
        b.iter(|| black_box(black_box(&old).get(black_box(SITE))))
    });
    g.bench_function("new/x_bit", |b| {
        b.iter(|| black_box(black_box(&new).x_bit(black_box(SITE))))
    });
    g.bench_function("old/x_bit", |b| {
        b.iter(|| black_box(black_box(&old).get_xbit(black_box(SITE))))
    });
    g.bench_function("new/z_bit", |b| {
        b.iter(|| black_box(black_box(&new).z_bit(black_box(SITE))))
    });
    g.bench_function("old/z_bit", |b| {
        b.iter(|| black_box(black_box(&old).get_zbit(black_box(SITE))))
    });
    g.finish();
}

fn hash(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/ordinary/hash_protocol/256");
    let bh = IdentityBuildHasher;
    let cold = NewWord::from(ordinary_string(WIDTH).as_str());
    let warm = cold.clone();
    assert_eq!(bh.hash_one(&warm), warm.key_hash());
    let old = OldWord::from(ordinary_string(WIDTH).as_str());
    assert_eq!(bh.hash_one(old), bh.hash_one(old));
    g.bench_function("new/new_only_cold_first_compute", |b| {
        b.iter_batched(
            || cold.clone(),
            |w| black_box(bh.hash_one(black_box(&w))),
            BatchSize::SmallInput,
        )
    });
    black_box(warm.key_hash());
    g.bench_function("new/warm", |b| {
        b.iter(|| black_box(bh.hash_one(black_box(&warm))))
    });
    // Old construction eagerly computes its digest, so no honest cold case exists.
    g.bench_function("old/warm", |b| {
        b.iter(|| black_box(bh.hash_one(black_box(&old))))
    });
    g.finish();
}

fn product(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/ordinary/product/256");
    let a = ordinary_string(WIDTH);
    let b: String = "ZYXI".chars().cycle().take(WIDTH).collect();
    let (na, nb) = (NewWord::from(a.as_str()), NewWord::from(b.as_str()));
    let oa = OldPhased::from(format!("+{a}").as_str());
    let ob = OldPhased::from(format!("+{b}").as_str());
    let (nw, np) = na.key_mul(&nb);
    let op = oa * ob;
    assert_eq!(nw.to_string(), op.word.to_string());
    assert_eq!(np.exponent(), op.phase);
    g.bench_function("new/product", |x| {
        x.iter(|| black_box(black_box(&na).key_mul(black_box(&nb))))
    });
    g.bench_function("old/product", |x| {
        x.iter(|| black_box(black_box(oa) * black_box(ob)))
    });
    g.finish();
}
