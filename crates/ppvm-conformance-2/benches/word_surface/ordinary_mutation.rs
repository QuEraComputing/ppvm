// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_traits::traits::PauliWordTrait;
use ppvm_traits_2::{Indexable, PauliBits};
use std::hint::black_box;

use super::common::*;

pub fn bench(c: &mut Criterion) {
    mutations(c);
    branches(c);
}

fn mutations(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/ordinary/mutate/256");
    let cold = NewWord::from(ordinary_string(WIDTH).as_str());
    let warm = cold.clone();
    black_box(warm.key_hash());
    let old = OldWord::from(ordinary_string(WIDTH).as_str());
    let mut nc = warm.clone();
    let mut oc = old;
    nc.set_x_bit(SITE, !nc.x_bit(SITE));
    oc.set_xbit(SITE, !oc.get_xbit(SITE));
    assert_eq!(nc.x_bit(SITE), oc.get_xbit(SITE));
    nc = warm.clone();
    oc = old;
    nc.set_z_bit(SITE, !nc.z_bit(SITE));
    oc.set_zbit(SITE, !oc.get_zbit(SITE));
    assert_eq!(nc.z_bit(SITE), oc.get_zbit(SITE));
    g.bench_function("new/new_only_set_x_bit_cold_cache", |b| {
        b.iter_batched(
            || cold.clone(),
            |mut w| {
                let v = !w.x_bit(SITE);
                w.set_x_bit(black_box(SITE), black_box(v));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("new/set_x_bit", |b| {
        b.iter_batched(
            || warm.clone(),
            |mut w| {
                let v = !w.x_bit(SITE);
                w.set_x_bit(black_box(SITE), black_box(v));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    // Old's public bit setter does not maintain its eager hash cache.
    g.bench_function("old/set_x_bit", |b| {
        b.iter_batched(
            || old,
            |mut w| {
                let v = !w.get_xbit(SITE);
                w.set_xbit(black_box(SITE), black_box(v));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("new/set_z_bit", |b| {
        b.iter_batched(
            || warm.clone(),
            |mut w| {
                let v = !w.z_bit(SITE);
                w.set_z_bit(black_box(SITE), black_box(v));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/set_z_bit", |b| {
        b.iter_batched(
            || old,
            |mut w| {
                let v = !w.get_zbit(SITE);
                w.set_zbit(black_box(SITE), black_box(v));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    // New has no whole-site setter; matched public mutation is the X/Z plane API.
    g.finish();
}

fn branches(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/ordinary/branch_key/256");
    let new = NewWord::from(ordinary_string(WIDTH).as_str());
    let old = OldWord::from(ordinary_string(WIDTH).as_str());
    let n1 = new.with_bits_toggled(SITE, true, true);
    let mut o1 = old;
    o1.set_xbit(SITE, !o1.get_xbit(SITE));
    o1.set_zbit(SITE, !o1.get_zbit(SITE));
    assert_eq!(n1.to_string(), o1.to_string());
    g.bench_function("new/one_site", |b| {
        b.iter(|| black_box(black_box(&new).with_bits_toggled(SITE, true, true)))
    });
    g.bench_function("old/one_site", |b| {
        b.iter(|| {
            let mut w = *black_box(&old);
            w.set_xbit(SITE, !w.get_xbit(SITE));
            w.set_zbit(SITE, !w.get_zbit(SITE));
            black_box(w)
        })
    });
    let n2 = new.with_bits_toggled2(SITE, true, true, SITE2, true, false);
    let mut o2 = old;
    o2.set_xbit(SITE, !o2.get_xbit(SITE));
    o2.set_zbit(SITE, !o2.get_zbit(SITE));
    o2.set_xbit(SITE2, !o2.get_xbit(SITE2));
    assert_eq!(n2.to_string(), o2.to_string());
    g.bench_function("new/two_site", |b| {
        b.iter(|| {
            black_box(black_box(&new).with_bits_toggled2(SITE, true, true, SITE2, true, false))
        })
    });
    g.bench_function("old/two_site", |b| {
        b.iter(|| {
            let mut w = *black_box(&old);
            w.set_xbit(SITE, !w.get_xbit(SITE));
            w.set_zbit(SITE, !w.get_zbit(SITE));
            w.set_xbit(SITE2, !w.get_xbit(SITE2));
            black_box(w)
        })
    });
    g.finish();
}
