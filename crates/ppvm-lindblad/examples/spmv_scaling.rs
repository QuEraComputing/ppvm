// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Tiny SpMV scaling check: builds a tridiagonal CSR matrix of a target
//! `nnz`, times `spmv_parallel` repeatedly inside rayon pools of varying
//! thread counts, and prints the wall time and speedup. Useful to confirm
//! the SpMV path itself scales before debugging higher-level paths.

use ppvm_lindblad::CsrMatrix;
use std::time::Instant;

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

fn main() {
    let n: usize = 100_000;
    let m = build_tridiagonal(n);
    let x: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
    let mut y = vec![0f64; n];
    let n_iters = 1000;
    println!(
        "n = {n}, nnz = {}, {n_iters} SpMVs per measurement",
        m.nnz()
    );
    println!();
    println!("{:>8} {:>14} {:>10}", "threads", "wall (ms)", "speedup");
    println!("{}", "-".repeat(36));

    let mut baseline_ms: Option<f64> = None;
    for nt in [1usize, 2, 3, 4, 6, 8, 10] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(nt)
            .build()
            .unwrap();
        // warm up
        pool.install(|| {
            for _ in 0..50 {
                m.spmv_parallel(&x, &mut y);
            }
        });
        let t0 = Instant::now();
        pool.install(|| {
            for _ in 0..n_iters {
                m.spmv_parallel(&x, &mut y);
            }
        });
        let wall_ms = t0.elapsed().as_secs_f64() * 1000.0;
        let speedup = baseline_ms.unwrap_or(wall_ms) / wall_ms;
        if baseline_ms.is_none() {
            baseline_ms = Some(wall_ms);
        }
        println!("{nt:>8} {wall_ms:>14.1} {speedup:>9.2}x");
    }
}
