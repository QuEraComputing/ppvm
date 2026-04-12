use std::time::Instant;

use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<Byte8F64<2>, u128>;

fn encode(tab: &mut Tab, qubits: &[usize]) {
    if qubits.len() == 17 {
        for i in [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16] {
            tab.sqrt_y(qubits[i]);
        }
        for [i, j] in [[1, 3], [7, 10], [12, 14], [13, 16]] {
            tab.cz(qubits[i], qubits[j]);
        }
        for i in [7, 16] {
            tab.sqrt_y_adj(qubits[i]);
        }
        for [i, j] in [[4, 7], [8, 10], [11, 14], [15, 16]] {
            tab.cz(qubits[i], qubits[j]);
        }
        for i in [4, 10, 14, 16] {
            tab.sqrt_y_adj(qubits[i]);
        }
        for [i, j] in [[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]] {
            tab.cz(qubits[i], qubits[j]);
        }
        for i in [3, 6, 9, 10, 12, 13] {
            tab.sqrt_y(qubits[i]);
        }
        for [i, j] in [[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]] {
            tab.cz(qubits[i], qubits[j]);
        }
        for i in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14] {
            tab.sqrt_y(qubits[i]);
        }
        for [i, j] in [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]] {
            tab.cz(qubits[i], qubits[j]);
        }
        for i in [0, 2, 5, 6, 8, 10, 12] {
            tab.sqrt_y_adj(qubits[i]);
        }
    }
}

fn main() {
    let n_qubits = 85;
    let n_runs = 20;

    let mut encoding_total = 0u128;
    let mut middle_gates_total = 0u128;
    let mut measure_total = 0u128;
    let mut full_total = 0u128;

    for _ in 0..n_runs {
        let full_start = Instant::now();

        let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
        let qubit_addrs: Vec<usize> = (0..n_qubits).collect();
        let ql: Vec<&[usize]> = qubit_addrs.chunks_exact(17).collect();

        // Phase 1: Encoding (H + T + encode per block)
        let t0 = Instant::now();
        for q in ql.iter() {
            let encoding_qubit = q[7];
            tab.h(encoding_qubit);
            tab.t(encoding_qubit);
            encode(&mut tab, q);
        }
        encoding_total += t0.elapsed().as_nanos();

        // Phase 2: Middle gates (sqrt_x, cz, sqrt_y, sqrt_x_adj layers)
        let t1 = Instant::now();
        for i in [0, 1, 4] {
            for q in ql[i] {
                tab.sqrt_x(*q);
            }
        }
        for (control, target) in ql[0].iter().zip(ql[1]) {
            tab.cz(*control, *target);
        }
        for (control, target) in ql[2].iter().zip(ql[3]) {
            tab.cz(*control, *target);
        }
        for q in ql[0] {
            tab.sqrt_y(*q);
        }
        for q in ql[3] {
            tab.sqrt_y(*q);
        }
        for (control, target) in ql[0].iter().zip(ql[2]) {
            tab.cz(*control, *target);
        }
        for (control, target) in ql[3].iter().zip(ql[4]) {
            tab.cz(*control, *target);
        }
        for q in ql[0] {
            tab.sqrt_x_adj(*q);
        }
        for (control, target) in ql[0].iter().zip(ql[4]) {
            tab.cz(*control, *target);
        }
        for (control, target) in ql[1].iter().zip(ql[3]) {
            tab.cz(*control, *target);
        }
        for i in 0..5 {
            for q in ql[i] {
                tab.sqrt_x_adj(*q);
            }
        }
        middle_gates_total += t1.elapsed().as_nanos();

        // Phase 3: Measurement
        let t2 = Instant::now();
        for i in 0..n_qubits {
            tab.measure(i);
        }
        measure_total += t2.elapsed().as_nanos();

        full_total += full_start.elapsed().as_nanos();
    }

    let div = n_runs as f64;
    println!("=== MSD 85-qubit profile ({n_runs} runs) ===");
    println!(
        "Encoding (H+T+encode): {:.1} µs",
        encoding_total as f64 / div / 1000.0
    );
    println!(
        "Middle gates:          {:.1} µs",
        middle_gates_total as f64 / div / 1000.0
    );
    println!(
        "Measurement:           {:.1} µs",
        measure_total as f64 / div / 1000.0
    );
    println!(
        "Full run:              {:.1} µs",
        full_total as f64 / div / 1000.0
    );
    println!();
    let total = encoding_total as f64;
    let mid = middle_gates_total as f64;
    let meas = measure_total as f64;
    let full = full_total as f64;
    println!("Encoding:    {:.1}%", total / full * 100.0);
    println!("Middle:      {:.1}%", mid / full * 100.0);
    println!("Measurement: {:.1}%", meas / full * 100.0);
}
