// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Layout baseline: the cost of the hot tableau paths, reported per operation.
//!
//! This exists to give the column-major / inverse-tableau rewrite a before/after
//! on the two things the layout controls:
//!
//! 1. **Clifford gates.** On the row-major `Vec<Row<A>>` layout every gate walks
//!    all `2n` rows and touches one machine word per plane per row, so the cost
//!    is `O(n)` *strided* — and the stride is `size_of::<Row<A>>()`, which is set
//!    by the compile-time storage width `A`, not by `n`. The `n = 85` rows are
//!    reported at both `U256` and `U2048` to price that padding directly.
//! 2. **`measure_all`.** Dominated by `compute_decomposition`, which is `O(n)`
//!    anticommutation probes plus `O(n)` full row multiplies per measured qubit.
//!
//! Run with `cargo run --release -p ppvm-tableau-2 --example layout_baseline`.

use std::hint::black_box;
use std::time::Instant;

use bnum::types::{U256, U512, U2048};

use ppvm_tableau_2::prelude::*;
use rand::SeedableRng;
use rand::rngs::SmallRng;

/// Time `f` over enough repetitions to clear timer noise, returning ns/op.
fn bench(reps: usize, mut f: impl FnMut()) -> f64 {
    // Warm up so the first-touch page faults land outside the measurement.
    for _ in 0..(reps / 10).max(1) {
        f();
    }
    let start = Instant::now();
    for _ in 0..reps {
        f();
    }
    start.elapsed().as_secs_f64() * 1e9 / reps as f64
}

/// One `n` of the report.
fn row<I: Bitstring>(n: usize) {
    let mut tab: GeneralizedTableau<I> = GeneralizedTableau::new(n, 1e-12);
    let h = bench(2000, || {
        for q in 0..n {
            tab.h(black_box(q));
        }
    }) / n as f64;

    let mut tab: GeneralizedTableau<I> = GeneralizedTableau::new(n, 1e-12);
    let s = bench(2000, || {
        for q in 0..n {
            tab.s(black_box(q));
        }
    }) / n as f64;

    let mut tab: GeneralizedTableau<I> = GeneralizedTableau::new(n, 1e-12);
    let cnot = bench(2000, || {
        for q in 0..n - 1 {
            tab.cnot(black_box(q), black_box(q + 1));
        }
    }) / (n - 1) as f64;

    // A GHZ-like frame so `measure_all` hits the random (case-a) branch rather
    // than measuring a product state.
    let seed_frame = |tab: &mut GeneralizedTableau<I>| {
        tab.h(0);
        for q in 0..n - 1 {
            tab.cnot(q, q + 1);
        }
    };
    let mut rng = SmallRng::seed_from_u64(0);
    let measure_all = bench(3, || {
        let mut tab: GeneralizedTableau<I> = GeneralizedTableau::new(n, 1e-12);
        seed_frame(&mut tab);
        black_box(tab.measure_all(&mut rng));
    }) / n as f64;

    println!("  n={n:<5}  h={h:>9.1}  s={s:>9.1}  cnot={cnot:>9.1}  measure={measure_all:>11.1}");
}

fn main() {
    println!("ns per operation (h/s/cnot are per gate; measure is per measured qubit)\n");
    // The frame is runtime-sized now, so there is no storage-width axis left to
    // sweep: the `n = 85` rows that used to differ by 4x between 128-bit and
    // 2048-bit storage are one row.
    row::<u128>(85);
    row::<U256>(200);
    row::<U512>(500);
    row::<U2048>(1889);
}
