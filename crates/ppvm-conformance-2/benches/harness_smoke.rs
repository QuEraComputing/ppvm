// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Phase 0 smoke benchmark: validates the comparative-benchmark harness
//! compiles and runs under Criterion. Real old-vs-new benchmarks (with a
//! printed ratio) are added per crate as the `-2` types come online.

use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_conformance_2::{random_circuit, seeded_rng};
use std::hint::black_box;

fn bench_generator(c: &mut Criterion) {
    c.bench_function("random_circuit/5q/1024", |b| {
        b.iter(|| {
            let mut rng = seeded_rng(0);
            black_box(random_circuit(&mut rng, 5, 1024))
        })
    });
}

criterion_group!(benches, bench_generator);
criterion_main!(benches);
