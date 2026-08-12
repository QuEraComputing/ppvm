// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_traits::traits::{PauliIter, PauliWordTrait};
use ppvm_traits_2::{IdentityBuildHasher, Indexable, PauliBits, Word};
use std::{hash::BuildHasher, hint::black_box};

use super::common::*;

pub fn bench(c: &mut Criterion) {
    reads(c);
    hash(c);
    // Deliberately no product group: neither lossy type defines a native
    // product. Projecting loss to identity would benchmark an ordinary word and
    // falsely label it as a lossy product.
}

fn reads(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/lossy/read/256");
    let s = lossy_string(WIDTH);
    let new = NewLossy::from(s.as_str());
    let old = OldLossy::from(s.as_str());
    assert_eq!(new.weight(), old.weight());
    assert_eq!(new.loss_weight(), old.loss_weight());
    assert_eq!(new.is_lost(SITE), old.get_lbit(SITE));
    let ni: String = new
        .iter()
        .map(|site| match site {
            ppvm_traits_2::LossySite::Present(p) => format!("{p:?}"),
            ppvm_traits_2::LossySite::Lost => "L".to_owned(),
        })
        .collect();
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
            black_box(&new).iter().for_each(|site| {
                black_box(site);
            })
        })
    });
    g.bench_function("old/iter_traverse", |b| {
        b.iter(|| {
            PauliIter::iter(black_box(&old)).for_each(|site| {
                black_box(site);
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
        b.iter(|| black_box(black_box(&new).loss_weight()))
    });
    g.bench_function("old/loss_weight", |b| {
        b.iter(|| black_box(black_box(&old).loss_weight()))
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
    g.bench_function("new/is_lost", |b| {
        b.iter(|| black_box(black_box(&new).is_lost(black_box(SITE))))
    });
    g.bench_function("old/is_lost", |b| {
        b.iter(|| black_box(black_box(&old).get_lbit(black_box(SITE))))
    });
    g.finish();
}

fn hash(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/lossy/hash_protocol/256");
    let bh = IdentityBuildHasher;
    let cold = NewLossy::from(lossy_string(WIDTH).as_str());
    let warm = cold.clone();
    assert_eq!(bh.hash_one(&warm), warm.key_hash());
    let old = OldLossy::from(lossy_string(WIDTH).as_str());
    assert_eq!(bh.hash_one(&old), bh.hash_one(old.clone()));
    g.bench_function("new/new_only_cold_all_components", |b| {
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
    // Old construction eagerly hashes all three planes.
    g.bench_function("old/warm", |b| {
        b.iter(|| black_box(bh.hash_one(black_box(&old))))
    });
    g.finish();
}
