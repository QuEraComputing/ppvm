// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Parallel SpMV scaling on a tridiagonal CSR matrix. Sweeps both matrix
//! size and thread count, so we can see where parallelism actually pays
//! off (compute-bound, L3-fitting) vs where DRAM bandwidth ceilings out
//! the speedup (large matrices). Each invocation runs inside a freshly
//! built rayon pool of the requested size.

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use ppvm_lindblad::{Csr, csr_from_triplets, spmv_parallel};

fn build_tridiagonal(n: usize) -> Csr {
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
    csr_from_triplets(n, &trips)
}

fn bench_spmv_scaling(c: &mut Criterion) {
    // (n, label) pairs: nnz ≈ 3n. ~16 bytes per nonzero (usize indices + f64)
    // + 8 bytes per row.
    //
    //   n =  100_000  → nnz ≈   300_000  →   ~5 MB  (fits L2/L3, compute-bound)
    //   n =  500_000  → nnz ≈ 1_500_000  →  ~24 MB  (L3 boundary on M-series)
    //   n = 2_000_000 → nnz ≈ 6_000_000  →  ~96 MB  (busts L3, bandwidth-bound)
    //   n = 8_000_000 → nnz ≈ 24_000_000 → ~384 MB  (well into DRAM)
    let sizes = [100_000usize, 500_000, 2_000_000, 8_000_000];

    for &n in &sizes {
        let m = build_tridiagonal(n);
        let x: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();

        let mut group = c.benchmark_group(format!("spmv_parallel/n={}", n));
        group.throughput(Throughput::Elements(m.nnz() as u64));
        group.sample_size(30);
        group.measurement_time(std::time::Duration::from_secs(3));
        group.warm_up_time(std::time::Duration::from_secs(1));

        for &threads in &[1usize, 2, 4, 8] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            group.bench_with_input(BenchmarkId::from_parameter(threads), &threads, |b, _| {
                b.iter_batched_ref(
                    || vec![0f64; n],
                    |y| pool.install(|| spmv_parallel(&m, &x, y)),
                    criterion::BatchSize::SmallInput,
                );
            });
        }

        group.finish();
    }
}

criterion_group!(benches, bench_spmv_scaling);
criterion_main!(benches);
