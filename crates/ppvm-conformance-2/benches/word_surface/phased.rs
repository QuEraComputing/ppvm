// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_traits::traits::{PauliIter, PauliWordTrait};
use ppvm_traits_2::Word;
use std::hint::black_box;

use super::common::*;

pub fn bench(c: &mut Criterion) {
    reads(c);
    phase_mutation(c);
    product(c);
    // New `Phased<W>` is deliberately non-Hash/non-Indexable. There is no
    // comparable hash benchmark; hashing its inner word would measure PauliWord.
    // New also intentionally exposes no whole-site setter on the wrapper, so
    // old `set`/`set_new` are old-only and omitted rather than faked.
    // Branch-key builders belong to indexable bare/lossy keys; neither phased
    // public surface has a comparable one-site or two-site branch operation.
}

fn reads(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/phased/read/256");
    let s = phased_string(WIDTH);
    let new = NewPhased::from(s.as_str());
    let old = OldPhased::from(s.as_str());
    assert_eq!(new.phase().exponent(), old.phase);
    assert_eq!(new.is_positive(), old.is_positive());
    assert_eq!(new.weight(), old.word.weight());
    let ni: String = new.iter().map(|p| format!("{p:?}")).collect();
    let oi: String = PauliIter::iter(&old.word).map(|p| p.to_string()).collect();
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
            PauliIter::iter(black_box(&old.word)).for_each(|p| {
                black_box(p);
            })
        })
    });
    g.bench_function("new/get", |b| {
        b.iter(|| black_box(black_box(&new).get(black_box(SITE))))
    });
    g.bench_function("old/get", |b| {
        b.iter(|| black_box(black_box(&old).get(black_box(SITE))))
    });
    g.bench_function("new/weight_delegate", |b| {
        b.iter(|| black_box(black_box(&new).weight()))
    });
    g.bench_function("old/weight_delegate", |b| {
        b.iter(|| black_box(black_box(&old).word.weight()))
    });
    g.bench_function("new/phase", |b| {
        b.iter(|| black_box(black_box(&new).phase()))
    });
    g.bench_function("old/phase", |b| b.iter(|| black_box(black_box(&old).phase)));
    g.bench_function("new/is_positive", |b| {
        b.iter(|| black_box(black_box(&new).is_positive()))
    });
    g.bench_function("old/is_positive", |b| {
        b.iter(|| black_box(black_box(&old).is_positive()))
    });
    g.finish();
}

fn phase_mutation(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/phased/add_phase/256");
    let new = NewPhased::from(phased_string(WIDTH).as_str());
    let old = OldPhased::from(phased_string(WIDTH).as_str());
    let mut nc = new.clone();
    let mut oc = old;
    nc.add_phase(ppvm_traits_2::Phase::NegI);
    oc.add_phase(3);
    assert_eq!(nc.phase().exponent(), oc.phase);
    g.bench_function("new/add_phase", |b| {
        b.iter_batched(
            || new.clone(),
            |mut w| {
                w.add_phase(black_box(ppvm_traits_2::Phase::NegI));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/add_phase", |b| {
        b.iter_batched(
            || old,
            |mut w| {
                w.add_phase(black_box(3));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.finish();
}

fn product(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/phased/product/256");
    let a = phased_string(WIDTH);
    let b = format!(
        "-{}",
        "ZYXI".chars().cycle().take(WIDTH).collect::<String>()
    );
    let (na, nb) = (NewPhased::from(a.as_str()), NewPhased::from(b.as_str()));
    let (oa, ob) = (OldPhased::from(a.as_str()), OldPhased::from(b.as_str()));
    let np = &na * &nb;
    let op = oa * ob;
    assert_eq!(np.to_string(), op.to_string());
    g.bench_function("new/product", |x| {
        x.iter(|| black_box(black_box(&na) * black_box(&nb)))
    });
    g.bench_function("old/product", |x| {
        x.iter(|| black_box(black_box(oa) * black_box(ob)))
    });
    g.finish();
}
