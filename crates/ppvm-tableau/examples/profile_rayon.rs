// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Profile the rayon vs sequential coefficient branching at various scales.
//!
//! This measures only the T-gate application phase (where coefficients branch),
//! isolating the exact code path that rayon parallelizes.
//!
//! Run without rayon:  cargo run --release -p ppvm-tableau --example profile_rayon
//! Run with rayon:     cargo run --release -p ppvm-tableau --example profile_rayon --features rayon

use std::time::Instant;

use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<Byte8F64<2>, u128>;

fn main() {
    // Use 128 qubits so we have enough room for many T gates on distinct qubits
    let n_qubits = 128;

    println!("=== T-gate branching: sequential vs rayon ({n_qubits} qubits) ===");
    println!(
        "{:>4}  {:>8}  {:>12}  {:>12}",
        "T#", "coeffs", "total (µs)", "last T (µs)"
    );

    for n_tgates in [1, 4, 8, 12, 14, 16, 18, 20] {
        // Prepare: H on each qubit that will get a T gate (ensures branching)
        let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
        for i in 0..n_tgates {
            tab.h(i);
        }

        // Warm up: one run to prime caches
        {
            let mut t = tab.fork(Some(42));
            for i in 0..n_tgates {
                t.t(i);
            }
        }

        // Measure: time all T gates together
        let n_runs = if n_tgates <= 16 { 5 } else { 3 };
        let mut total_all_ns = 0u128;
        let mut total_last_ns = 0u128;
        let mut coeff_count = 0;

        for _ in 0..n_runs {
            let mut t = tab.fork(Some(42));

            // Apply T gates 0..n-1 (setup, not timed individually)
            let t0 = Instant::now();
            for i in 0..n_tgates.saturating_sub(1) {
                t.t(i);
            }
            let setup_ns = t0.elapsed().as_nanos();

            // Time only the last T gate (largest coefficient set)
            let t1 = Instant::now();
            if n_tgates > 0 {
                t.t(n_tgates - 1);
            }
            let last_ns = t1.elapsed().as_nanos();

            total_all_ns += setup_ns + last_ns;
            total_last_ns += last_ns;
            coeff_count = t.coefficients.len();
        }

        println!(
            "{:>4}  {:>8}  {:>12.1}  {:>12.1}",
            n_tgates,
            coeff_count,
            total_all_ns as f64 / n_runs as f64 / 1000.0,
            total_last_ns as f64 / n_runs as f64 / 1000.0,
        );
    }
}
