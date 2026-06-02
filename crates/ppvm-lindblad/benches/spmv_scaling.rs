// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Parallel SpMV scaling on a tridiagonal CSR matrix. One bench per thread
//! count, with each invocation running inside a freshly-built rayon pool of
//! that size to isolate the SpMV path from outer-loop work.

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use ppvm_lindblad::CsrMatrix;

fn build_tridiagonal(n: usize) -> CsrMatrix {
    let mut trips = Vec::with_capacity(3 * n);
    for i in 0..n {
        trips.push((i, i, -2.0));
        if i > 0 {
            trips.push((i, i - 1, 1.0));
        }
        if i + 1 < n {
            trips.push((i, i + 1, 1.0));
        }
    }
    CsrMatrix::from_triplets(n, &trips)
}

fn bench_spmv_scaling(c: &mut Criterion) {
    let n: usize = 100_000;
    let m = build_tridiagonal(n);
    let x: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();

    let mut group = c.benchmark_group("spmv_parallel");
    group.throughput(Throughput::Elements(m.nnz() as u64));

    for &threads in &[1usize, 2, 3, 4, 6, 8, 10] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        group.bench_with_input(BenchmarkId::from_parameter(threads), &threads, |b, _| {
            b.iter_batched_ref(
                || vec![0f64; n],
                |y| pool.install(|| m.spmv_parallel(&x, y)),
                criterion::BatchSize::SmallInput,
            );
        });
    }

    group.finish();
}

criterion_group!(benches, bench_spmv_scaling);
criterion_main!(benches);
