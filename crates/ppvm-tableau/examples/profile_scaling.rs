// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::time::Instant;

use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<Byte8F64<2>, usize>;

fn main() {
    let n_runs = 10;

    // Profile 1: Scaling with qubit count (Clifford-heavy)
    println!("=== Clifford gate scaling (CNOT chain) ===");
    for n_qubits in [32, 64, 96, 128] {
        let mut total_ns = 0u128;
        for _ in 0..n_runs {
            let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
            tab.h(0);
            let t0 = Instant::now();
            for i in 0..n_qubits - 1 {
                tab.cnot(i, i + 1);
            }
            total_ns += t0.elapsed().as_nanos();
        }
        println!(
            "  n={:3}: {:.1} µs ({} cnots, {} rows)",
            n_qubits,
            total_ns as f64 / n_runs as f64 / 1000.0,
            n_qubits - 1,
            2 * n_qubits
        );
    }

    // Profile 2: T-gate coefficient branching scaling
    println!("\n=== T-gate branching scaling (85 qubits) ===");
    for n_tgates in [1, 4, 8, 12, 16] {
        let mut tab: Tab = GeneralizedTableau::new(85, 1e-10);
        for i in 0..n_tgates.min(85) {
            tab.h(i); // ensure non-deterministic for branching
        }

        let mut total_ns = 0u128;
        let mut coeff_count = 0;
        for _ in 0..n_runs {
            let mut t = tab.fork(Some(42));
            let t0 = Instant::now();
            for i in 0..n_tgates {
                t.t(i % 85);
            }
            total_ns += t0.elapsed().as_nanos();
            coeff_count = t.coefficients.len();
        }
        println!(
            "  t={:2}: {:.1} µs (coeffs={})",
            n_tgates,
            total_ns as f64 / n_runs as f64 / 1000.0,
            coeff_count,
        );
    }

    // Profile 3: Measurement after T gates
    println!("\n=== Measurement scaling after T gates (85 qubits) ===");
    for n_tgates in [1, 4, 8, 12] {
        let mut tab: Tab = GeneralizedTableau::new(85, 1e-10);
        for i in 0..n_tgates.min(85) {
            tab.h(i);
        }
        for i in 0..n_tgates {
            tab.t(i % 85);
        }

        let mut total_ns = 0u128;
        for _ in 0..n_runs {
            let mut t = tab.fork(Some(42));
            let t0 = Instant::now();
            for i in 0..85 {
                t.measure(i);
            }
            total_ns += t0.elapsed().as_nanos();
        }
        println!(
            "  t={:2}: {:.1} µs (85 measures, initial_coeffs={})",
            n_tgates,
            total_ns as f64 / n_runs as f64 / 1000.0,
            tab.coefficients.len(),
        );
    }

    // Profile 4: Individual gate cost at different qubit counts
    println!("\n=== Single gate cost (per-call overhead test) ===");
    for n_qubits in [32, 64, 128] {
        let tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
        let n_calls = 1000;
        let mut total_ns = 0u128;
        for _ in 0..n_runs {
            let mut t = tab.fork(Some(42));
            let t0 = Instant::now();
            for _ in 0..n_calls {
                t.h(0);
            }
            total_ns += t0.elapsed().as_nanos();
        }
        println!(
            "  n={:3}: {:.0} ns/call (H gate, {} rows)",
            n_qubits,
            total_ns as f64 / (n_runs * n_calls) as f64,
            2 * n_qubits
        );
    }
}
