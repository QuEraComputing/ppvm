// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Per-qubit timing breakdown of `measure_all` on the MSD workload.
//!
//! Builds the same 85-qubit MSD state used by the `measure-all` bench, then on
//! a fresh clone times each `measure(q)` individually. Repeats `n_runs` times
//! and prints per-qubit median time. Also reports total time and the
//! coefficient count after the circuit (so the case-a HashMap cost is visible
//! in context).
//!
//! Run: `cargo run -p ppvm-tableau --example profile_measure_all --release`

use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::measure_all::LossyMeasureAll;
use ppvm_tableau::prelude::*;
use ppvm_traits::traits::LossyMeasure;
use std::time::{Duration, Instant};

type Tab = GeneralizedTableau<Byte8F64<2>, u128>;

fn build_msd_state() -> Tab {
    let qubits_per_code_block = 17;
    let n_qubits = qubits_per_code_block * 5;
    let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
    let qubit_addrs: Vec<usize> = (0..n_qubits).collect();
    let ql: Vec<&[usize]> = qubit_addrs.chunks_exact(qubits_per_code_block).collect();

    for q in ql.iter() {
        let encoding_qubit = q[7];
        tab.h(encoding_qubit);
        tab.t(encoding_qubit);
        encode(&mut tab, q);
    }

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
    for block in ql.iter().take(5) {
        for q in *block {
            tab.sqrt_x_adj(*q);
        }
    }

    tab
}

fn encode(tab: &mut Tab, qubits: &[usize]) {
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

fn median(mut xs: Vec<Duration>) -> Duration {
    xs.sort();
    xs[xs.len() / 2]
}

fn main() {
    let base = build_msd_state();
    let n = base.n_qubits();
    let n_coeffs = base.coefficients.len();
    let n_runs = 200;

    println!("MSD state: {} qubits, {} coefficients", n, n_coeffs);
    println!();

    // ---- End-to-end per-shot timing (fork + measure_all) ----
    // This matches how samplers actually use the API: one shot = one fork
    // followed by one measure_all on the clone.
    let mut shot_times = Vec::with_capacity(n_runs);
    for _ in 0..n_runs {
        let start = Instant::now();
        let mut t = base.fork(Some(42));
        let _ = t.measure_all();
        shot_times.push(start.elapsed());
    }
    let shot_median = median(shot_times.clone());

    // ---- Fork-only timing (so we can read off fork's share) ----
    let mut fork_times = Vec::with_capacity(n_runs);
    for _ in 0..n_runs {
        let start = Instant::now();
        let t = base.fork(Some(42));
        fork_times.push(start.elapsed());
        std::hint::black_box(t);
    }
    let fork_median = median(fork_times);

    let measure_only_median = shot_median.saturating_sub(fork_median);
    println!(
        "per-shot (fork + measure_all) median over {n_runs} runs: {:?}",
        shot_median
    );
    println!(
        "    fork only:           {:?}  ({:5.1}%)",
        fork_median,
        100.0 * fork_median.as_secs_f64() / shot_median.as_secs_f64()
    );
    println!(
        "    measure_all only:    {:?}  ({:5.1}%)",
        measure_only_median,
        100.0 * measure_only_median.as_secs_f64() / shot_median.as_secs_f64()
    );
    println!(
        "    avg per qubit (measure_all): {:?}",
        measure_only_median / n as u32
    );
    println!();

    // ---- Per-qubit breakdown ----
    // For each qubit index, record the time to measure THAT qubit when
    // measuring sequentially 0..n on a fresh clone. Each run uses the same
    // seed so the case-a vs case-b classification is deterministic.
    let mut per_qubit_runs: Vec<Vec<Duration>> = vec![Vec::with_capacity(n_runs); n];
    for _ in 0..n_runs {
        let mut t = base.fork(Some(42));
        for q in 0..n {
            let start = Instant::now();
            let _ = LossyMeasure::measure(&mut t, q);
            per_qubit_runs[q].push(start.elapsed());
        }
    }
    let per_qubit_medians: Vec<Duration> = per_qubit_runs.into_iter().map(median).collect();

    // Bucket by speed to surface the bimodal distribution.
    let mut sorted = per_qubit_medians.clone();
    sorted.sort();
    let p50 = sorted[n / 2];
    let p90 = sorted[(n * 9) / 10];
    let pmin = *sorted.first().unwrap();
    let pmax = *sorted.last().unwrap();

    println!("per-qubit median time  (sorted ascending):");
    println!("    min  : {:?}", pmin);
    println!("    p50  : {:?}", p50);
    println!("    p90  : {:?}", p90);
    println!("    max  : {:?}", pmax);
    println!();

    // ---- Classify qubits into "cheap" and "expensive" buckets ----
    // Case-b (Z is already a stabilizer, deterministic outcome) is the cheap
    // path; case-a is the expensive path with HashMap + tableau row updates.
    // Use 2× the minimum as a heuristic cutoff.
    let cutoff = pmin * 2;
    let cheap: Vec<usize> = (0..n).filter(|&q| per_qubit_medians[q] <= cutoff).collect();
    let expensive: Vec<usize> = (0..n).filter(|&q| per_qubit_medians[q] > cutoff).collect();

    let cheap_total: Duration = cheap.iter().map(|&q| per_qubit_medians[q]).sum();
    let expensive_total: Duration = expensive.iter().map(|&q| per_qubit_medians[q]).sum();
    let grand_total: Duration = per_qubit_medians.iter().sum();

    println!(
        "cheap path (likely case-b): {:3} qubits, total {:?} ({:5.1}% of work)",
        cheap.len(),
        cheap_total,
        100.0 * cheap_total.as_secs_f64() / grand_total.as_secs_f64()
    );
    println!(
        "slow path  (likely case-a): {:3} qubits, total {:?} ({:5.1}% of work)",
        expensive.len(),
        expensive_total,
        100.0 * expensive_total.as_secs_f64() / grand_total.as_secs_f64()
    );
    println!();

    println!("per-qubit median times (ns), in measurement order:");
    for (q, t) in per_qubit_medians.iter().enumerate() {
        let bar_units = (t.as_nanos() / 200).min(60) as usize;
        let bar: String = std::iter::repeat_n('#', bar_units).collect();
        println!("    q={:3}  {:>7} ns  {}", q, t.as_nanos(), bar);
    }
}
