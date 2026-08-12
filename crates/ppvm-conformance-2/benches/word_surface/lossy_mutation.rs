// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_traits::traits::PauliWordTrait;
use ppvm_traits_2::{Indexable, PauliBits};
use std::hint::black_box;

use super::common::*;

pub fn bench(c: &mut Criterion) {
    mutations(c);
}

fn mutations(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/lossy/mutate/256");
    let cold = NewLossy::from(lossy_string(WIDTH).as_str());
    let warm = cold.clone();
    black_box(warm.key_hash());
    let old = OldLossy::from(lossy_string(WIDTH).as_str());
    let mut nc = warm.clone();
    let mut oc = old.clone();
    nc.set(
        SITE,
        ppvm_traits_2::LossySite::Present(ppvm_traits_2::Pauli::X),
    );
    oc.set(SITE, ppvm_traits::char::Pauli::X);
    assert_eq!(nc.to_string(), oc.to_string());
    nc = warm.clone();
    oc = old.clone();
    nc.set_lost(SITE);
    oc.set(SITE, ppvm_traits::char::Pauli::L);
    assert_eq!(nc.to_string(), oc.to_string());
    const LOST: usize = 124;
    nc = warm.clone();
    oc = old.clone();
    nc.clear_loss(LOST);
    oc.set(LOST, ppvm_traits::char::Pauli::I);
    assert_eq!(nc.to_string(), oc.to_string());
    nc = warm.clone();
    oc = old.clone();
    nc.set_x_bit(SITE, false);
    oc.set_xbit(SITE, false);
    assert_eq!(nc.to_string(), oc.to_string());
    nc = warm.clone();
    oc = old.clone();
    nc.set_z_bit(SITE, false);
    oc.set_zbit(SITE, false);
    assert_eq!(nc.to_string(), oc.to_string());

    g.bench_function("new/set_present", |b| {
        b.iter_batched(
            || warm.clone(),
            |mut w| {
                w.set(
                    black_box(SITE),
                    black_box(ppvm_traits_2::LossySite::Present(ppvm_traits_2::Pauli::X)),
                );
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/set_present", |b| {
        b.iter_batched(
            || old.clone(),
            |mut w| {
                w.set(black_box(SITE), black_box(ppvm_traits::char::Pauli::X));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("new/set_lost", |b| {
        b.iter_batched(
            || warm.clone(),
            |mut w| {
                w.set_lost(black_box(SITE));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/set_lost", |b| {
        b.iter_batched(
            || old.clone(),
            |mut w| {
                w.set(black_box(SITE), black_box(ppvm_traits::char::Pauli::L));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("new/clear_loss", |b| {
        b.iter_batched(
            || {
                let w = warm.clone();
                black_box(w.key_hash());
                w
            },
            |mut w| {
                w.clear_loss(black_box(LOST));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/clear_loss", |b| {
        b.iter_batched(
            || old.clone(),
            |mut w| {
                w.set(black_box(LOST), black_box(ppvm_traits::char::Pauli::I));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("new/set_x_bit", |b| {
        b.iter_batched(
            || cold.clone(),
            |mut w| {
                w.set_x_bit(black_box(SITE), black_box(false));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/set_x_bit", |b| {
        b.iter_batched(
            || old.clone(),
            |mut w| {
                w.set_xbit(black_box(SITE), black_box(false));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("new/set_z_bit", |b| {
        b.iter_batched(
            || warm.clone(),
            |mut w| {
                w.set_z_bit(black_box(SITE), black_box(false));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/set_z_bit", |b| {
        b.iter_batched(
            || old.clone(),
            |mut w| {
                w.set_zbit(black_box(SITE), black_box(false));
                black_box(w)
            },
            BatchSize::SmallInput,
        )
    });
    g.finish();
}
