// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// Ad-hoc profiling binary to measure measurement time scaling.
/// Run: cargo run -p ppvm-tableau --example profile_measure --release
use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;
use std::time::Instant;

type Tab = GeneralizedTableau<Byte8F64<2>, u128>;

fn build_circuit(n_qubits: usize, n_t_gates: usize) -> Tab {
    let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
    for i in 0..n_t_gates {
        tab.h(i);
        tab.t(i);
    }
    for i in 0..n_qubits - 1 {
        tab.cz([i, i + 1]);
    }
    tab
}

fn main() {
    let n_qubits = 85;

    for n_t in [5, 8, 10, 12, 14] {
        let tab = build_circuit(n_qubits, n_t);
        let n_coeffs = tab.coefficients.len();

        let n_runs = if n_t <= 10 { 20 } else { 5 };
        let mut total_times = Vec::with_capacity(n_runs);

        for _ in 0..n_runs {
            let mut t = tab.fork(Some(42));
            let start = Instant::now();
            for q in 0..n_qubits {
                let _ = t.measure(q);
            }
            total_times.push(start.elapsed());
        }

        total_times.sort();
        let median = total_times[total_times.len() / 2];
        println!(
            "n_t={:2}, coeffs={:6}, measurement({}q): {:?}",
            n_t, n_coeffs, n_qubits, median
        );
    }
}
