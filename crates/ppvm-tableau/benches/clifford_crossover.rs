// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Crossover bench: at fixed tableau size, sweep the *batch size* `k` and
//! compare a flat loop of `k` single-gate calls against the batched call.
//!
//! Goal is to find the smallest `k` at which the bitmask / fused batch pays
//! off versus the loop. The fixed tableau size is 128 (two `u64` words, so
//! batches up to `k = 64` stay in one storage word; `k > 64` spans both).
//!
//! Three representative boolean-op patterns:
//! - `h`: bit swap + phase from `x & z`
//! - `x`: phase from `z` only, no bit update (simplest)
//! - `s`: bit update `z ^= x` + phase from `x & z`

use std::time::Duration;

use criterion::{BatchSize, BenchmarkId, Criterion, criterion_group, criterion_main};
use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_pauli_sum::prelude::*;
use ppvm_tableau::prelude::*;

type Tab = Tableau<Byte8F64<2>>;

const N: usize = 128;
const BATCH_SIZES: &[usize] = &[1, 2, 4, 8];

fn configure() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50)
}

fn setup() -> Tab {
    let mut tab = Tab::new(N);
    for q in 0..8 {
        tab.h(q);
    }
    for q in (0..8).step_by(2) {
        tab.s(q);
    }
    tab
}

macro_rules! crossover_group {
    ($c:expr, $tab:expr, $label:literal, $single:ident, $batch:ident) => {{
        let mut group = $c.benchmark_group(concat!("clifford-crossover/", $label));
        for &k in BATCH_SIZES {
            let indices: Vec<usize> = (0..k).collect();
            group.bench_with_input(BenchmarkId::new("loop", k), &indices, |b, idx| {
                b.iter_batched_ref(
                    || $tab.clone(),
                    |t| {
                        for &q in idx {
                            t.$single(q);
                        }
                    },
                    BatchSize::SmallInput,
                );
            });
            group.bench_with_input(BenchmarkId::new("batch", k), &indices, |b, idx| {
                b.iter_batched_ref(|| $tab.clone(), |t| t.$batch(idx), BatchSize::SmallInput);
            });
        }
        group.finish();
    }};
}

fn bench_crossover(c: &mut Criterion) {
    let tab = setup();
    crossover_group!(c, tab, "h", h, h_many);
    crossover_group!(c, tab, "x", x, x_many);
    crossover_group!(c, tab, "s", s, s_many);
}

criterion_group! {
    name = benches;
    config = configure();
    targets = bench_crossover
}
criterion_main!(benches);
