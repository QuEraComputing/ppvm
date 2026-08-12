// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_traits::traits::PauliWordTrait;
use ppvm_traits_2::{Indexable, PauliBits};
use std::hint::black_box;

use super::common::*;

pub fn bench(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/lossy/branch_key/256");
    // A branch key is immediately consumed by a hash-map probe. Include that
    // required hash on both sides: new computes it lazily with `key_hash`, while
    // old's setters require the caller to finish the key with `rehash`.
    let new = NewLossy::from(lossy_string(WIDTH).as_str());
    let old = OldLossy::from(lossy_string(WIDTH).as_str());
    let n1 = new.toggled_bits(SITE, true, true);
    let mut o1 = old.clone();
    o1.set_xbit(SITE, !o1.get_xbit(SITE));
    o1.set_zbit(SITE, !o1.get_zbit(SITE));
    assert_eq!(n1.to_string(), o1.to_string());
    g.bench_function("new/one_site_clone_then_bits", |b| {
        b.iter(|| {
            let w = black_box(&new).toggled_bits(SITE, true, true);
            black_box(w.key_hash());
            black_box(w)
        })
    });
    g.bench_function("old/one_site_clone_then_bits", |b| {
        b.iter(|| {
            let mut w = black_box(&old).clone();
            w.set_xbit(SITE, !w.get_xbit(SITE));
            w.set_zbit(SITE, !w.get_zbit(SITE));
            w.rehash();
            black_box(w)
        })
    });

    let n2 = new.toggled_bits2(SITE, true, true, SITE2, true, false);
    let mut o2 = old.clone();
    o2.set_xbit(SITE, !o2.get_xbit(SITE));
    o2.set_zbit(SITE, !o2.get_zbit(SITE));
    o2.set_xbit(SITE2, !o2.get_xbit(SITE2));
    assert_eq!(n2.to_string(), o2.to_string());
    g.bench_function("new/two_site_clone_then_bits", |b| {
        b.iter(|| {
            let w = black_box(&new).toggled_bits2(SITE, true, true, SITE2, true, false);
            black_box(w.key_hash());
            black_box(w)
        })
    });
    g.bench_function("old/two_site_clone_then_bits", |b| {
        b.iter(|| {
            let mut w = black_box(&old).clone();
            w.set_xbit(SITE, !w.get_xbit(SITE));
            w.set_zbit(SITE, !w.get_zbit(SITE));
            w.set_xbit(SITE2, !w.get_xbit(SITE2));
            w.rehash();
            black_box(w)
        })
    });
    g.finish();
}
