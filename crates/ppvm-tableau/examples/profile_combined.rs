// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Profile the time breakdown of a fused circuit with variable T gates.

use std::time::Instant;

use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<Byte8F64<2>, u128>;

fn main() {
    let n_qubits = 85;

    println!(
        "{:>4}  {:>10}  {:>10}  {:>10}  {:>10}  {:>6}  {:>6}  {:>6}",
        "T#", "clifford", "t-gates", "measure", "total", "%clif", "%t", "%meas"
    );

    for n_tgates in [4, 8, 12, 14, 16] {
        let block1: Vec<usize> = (0..17).collect();
        let block2: Vec<usize> = (17..34).collect();

        let n_runs = if n_tgates <= 12 { 5 } else { 2 };
        let mut clif_ns = 0u128;
        let mut tgate_ns = 0u128;
        let mut meas_ns = 0u128;

        for _ in 0..n_runs {
            let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);

            // Clifford layer 1 (fused)
            let t0 = Instant::now();
            block1.iter().for_each(|&i| tab.sqrt_y(i));
            block2.iter().for_each(|&i| tab.sqrt_x(i));
            tab.cz_block_pairs(0, 17, 17);
            for i in 0..n_tgates {
                tab.h(i);
            }
            clif_ns += t0.elapsed().as_nanos();

            // T gates
            let t1 = Instant::now();
            for i in 0..n_tgates {
                tab.t(i);
            }
            tgate_ns += t1.elapsed().as_nanos();

            // Clifford layer 2 (fused)
            let t2 = Instant::now();
            block1.iter().for_each(|&i| tab.sqrt_x_dag(i));
            block2.iter().for_each(|&i| tab.sqrt_y_dag(i));
            clif_ns += t2.elapsed().as_nanos();

            // Measure
            let t3 = Instant::now();
            for i in 0..n_qubits {
                tab.measure(i);
            }
            meas_ns += t3.elapsed().as_nanos();
        }

        let d = n_runs as f64;
        let c = clif_ns as f64 / d / 1000.0;
        let t = tgate_ns as f64 / d / 1000.0;
        let m = meas_ns as f64 / d / 1000.0;
        let total = c + t + m;
        println!(
            "{:>4}  {:>9.0}µ  {:>9.0}µ  {:>9.0}µ  {:>9.0}µ  {:>5.1}%  {:>5.1}%  {:>5.1}%",
            n_tgates,
            c,
            t,
            m,
            total,
            c / total * 100.0,
            t / total * 100.0,
            m / total * 100.0,
        );
    }
}
