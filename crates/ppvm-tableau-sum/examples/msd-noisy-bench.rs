// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

// Deterministic timing harness for the msd-noisy build + sample workload.
// Mirrors examples/msd-noisy.rs but uses a fixed seed, runs the build several
// times (median), and asserts the final branch count so an optimization that
// silently changes the math is caught. Used by the autotune experiment
// `docs/autotune/2026-06-23-tableau-sum-build`.

use std::time::Instant;

use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;
use ppvm_tableau_sum::data::GeneralizedTableauSum;
use ppvm_tableau_sum::storage::EntryStore;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

type GTabSum = GeneralizedTableauSum<Byte8F64<2>, u128>;

fn encode(tab: &mut GTabSum, qubits: &[usize], p_loss: f64, p_depolarize: f64) {
    if qubits.len() == 17 {
        for i in [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16] {
            tab.sqrt_y(qubits[i]);
            tab.loss_channel(qubits[i], p_loss);
            tab.depolarize1(qubits[i], p_depolarize);
        }
        for [i, j] in [[1, 3], [7, 10], [12, 14], [13, 16]] {
            tab.cz(qubits[i], qubits[j]);
            tab.loss_channel(qubits[i], p_loss);
            tab.loss_channel(qubits[j], p_loss);
            tab.depolarize1(qubits[i], p_depolarize);
            tab.depolarize1(qubits[j], p_depolarize);
        }
        for i in [7, 16] {
            tab.sqrt_y_dag(qubits[i]);
            tab.loss_channel(qubits[i], p_loss);
            tab.depolarize1(qubits[i], p_depolarize);
        }
        for [i, j] in [[4, 7], [8, 10], [11, 14], [15, 16]] {
            tab.cz(qubits[i], qubits[j]);
            tab.loss_channel(qubits[i], p_loss);
            tab.loss_channel(qubits[j], p_loss);
            tab.depolarize1(qubits[i], p_depolarize);
            tab.depolarize1(qubits[j], p_depolarize);
        }
        for i in [4, 10, 14, 16] {
            tab.sqrt_y_dag(qubits[i]);
            tab.loss_channel(qubits[i], p_loss);
            tab.depolarize1(qubits[i], p_depolarize);
        }
        for [i, j] in [[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]] {
            tab.cz(qubits[i], qubits[j]);
            tab.loss_channel(qubits[i], p_loss);
            tab.loss_channel(qubits[j], p_loss);
            tab.depolarize1(qubits[i], p_depolarize);
            tab.depolarize1(qubits[j], p_depolarize);
        }
        for i in [3, 6, 9, 10, 12, 13] {
            tab.sqrt_y(qubits[i]);
            tab.loss_channel(qubits[i], p_loss);
            tab.depolarize1(qubits[i], p_depolarize);
        }
        for [i, j] in [[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]] {
            tab.cz(qubits[i], qubits[j]);
            tab.loss_channel(qubits[i], p_loss);
            tab.loss_channel(qubits[j], p_loss);
            tab.depolarize1(qubits[i], p_depolarize);
            tab.depolarize1(qubits[j], p_depolarize);
        }
        for i in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14] {
            tab.sqrt_y(qubits[i]);
            tab.loss_channel(qubits[i], p_loss);
            tab.depolarize1(qubits[i], p_depolarize);
        }
        for [i, j] in [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]] {
            tab.cz(qubits[i], qubits[j]);
            tab.loss_channel(qubits[i], p_loss);
            tab.loss_channel(qubits[j], p_loss);
            tab.depolarize1(qubits[i], p_depolarize);
            tab.depolarize1(qubits[j], p_depolarize);
        }
        for i in [0, 2, 5, 6, 8, 10, 12] {
            tab.sqrt_y_dag(qubits[i]);
            tab.loss_channel(qubits[i], p_loss);
            tab.depolarize1(qubits[i], p_depolarize);
        }
    }
}

fn build(seed: u64) -> GTabSum {
    let n_qubits = 85;
    let p_loss = 1e-4;
    let p_depolarize = 1e-4;
    let sum_cutoff = 1e-7;

    let mut tab: GTabSum = GeneralizedTableauSum::new_with_seed(n_qubits, 1e-10, sum_cutoff, seed);
    let qubit_addrs: Vec<usize> = (0..n_qubits).collect();
    let ql: Vec<&[usize]> = qubit_addrs.chunks_exact(17).collect();

    for q in ql.iter() {
        let encoding_qubit = q[7];
        tab.h(encoding_qubit);
        tab.loss_channel(encoding_qubit, p_loss);
        tab.depolarize1(encoding_qubit, p_depolarize);
        tab.t(encoding_qubit);
        tab.loss_channel(encoding_qubit, p_loss);
        tab.depolarize1(encoding_qubit, p_depolarize);
        encode(&mut tab, q, p_loss, p_depolarize);
    }

    for i in [0, 1, 4] {
        for q in ql[i] {
            tab.sqrt_x(*q);
            tab.loss_channel(*q, p_loss);
            tab.depolarize1(*q, p_depolarize);
        }
    }
    for (control, target) in ql[0].iter().zip(ql[1]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize1(*control, p_depolarize);
        tab.depolarize1(*target, p_depolarize);
    }
    for (control, target) in ql[2].iter().zip(ql[3]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize1(*control, p_depolarize);
        tab.depolarize1(*target, p_depolarize);
    }
    for q in ql[0] {
        tab.sqrt_y(*q);
        tab.loss_channel(*q, p_loss);
        tab.depolarize1(*q, p_depolarize);
    }
    for q in ql[3] {
        tab.sqrt_y(*q);
        tab.loss_channel(*q, p_loss);
        tab.depolarize1(*q, p_depolarize);
    }
    for (control, target) in ql[0].iter().zip(ql[2]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize1(*control, p_depolarize);
        tab.depolarize1(*target, p_depolarize);
    }
    for (control, target) in ql[3].iter().zip(ql[4]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize1(*control, p_depolarize);
        tab.depolarize1(*target, p_depolarize);
    }
    for q in ql[0] {
        tab.sqrt_x_dag(*q);
        tab.loss_channel(*q, p_loss);
        tab.depolarize1(*q, p_depolarize);
    }
    for (control, target) in ql[0].iter().zip(ql[4]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize1(*control, p_depolarize);
        tab.depolarize1(*target, p_depolarize);
    }
    for (control, target) in ql[1].iter().zip(ql[3]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize1(*control, p_depolarize);
        tab.depolarize1(*target, p_depolarize);
    }
    for block in ql.iter().take(5) {
        for q in *block {
            tab.sqrt_x_dag(*q);
            tab.loss_channel(*q, p_loss);
            tab.depolarize1(*q, p_depolarize);
        }
    }

    tab
}

fn main() {
    const BUILD_RUNS: usize = 5;
    const N_SHOTS: usize = 20000;
    const SEED: u64 = 12345;

    let mut build_times_ms: Vec<f64> = Vec::new();
    let mut branches = 0usize;
    for r in 0..BUILD_RUNS {
        let now = Instant::now();
        let tab = build(SEED + r as u64);
        let ms = now.elapsed().as_secs_f64() * 1e3;
        build_times_ms.push(ms);
        branches = tab.len();
    }
    build_times_ms.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_ms = build_times_ms[build_times_ms.len() / 2];
    let min_ms = build_times_ms[0];

    // Accuracy fingerprint: the optimizations under test must not change the
    // math, so the multiset of branch probabilities must be invariant (up to
    // float summation-order noise). Capture sum(p), sum(p^2) (participation
    // ratio), and the top-5 probabilities from a fresh deterministic build.
    let tab_acc = build(SEED);
    let mut probs: Vec<f64> = tab_acc.entries.iter().map(|(_, p)| *p).collect();
    probs.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let sum_p: f64 = probs.iter().sum();
    let sum_p2: f64 = probs.iter().map(|p| p * p).sum();
    let top5: Vec<f64> = probs.iter().take(5).copied().collect();

    // Sample timing on a fresh build.
    let mut tab = build(SEED);
    let mut sampler = tab.sampler();
    let now = Instant::now();
    sampler.sample_shots(N_SHOTS);
    let per_shot_ns = now.elapsed().as_nanos() as f64 / N_SHOTS as f64;

    println!("branches       = {}", branches);
    println!("build_min_ms   = {:.1}", min_ms);
    println!("build_median_ms= {:.1}", median_ms);
    println!("per_shot_ns    = {:.1}", per_shot_ns);
    println!("sum_p          = {:.12}", sum_p);
    println!("sum_p2         = {:.12}", sum_p2);
    println!("top5_p         = {:?}", top5);
    println!("all_build_ms   = {:?}", build_times_ms);

    // Accuracy guard: the optimizations under test must not change the math,
    // so the final branch count must stay at the baseline value.
    const EXPECTED_BRANCHES: usize = 2025;
    if branches != EXPECTED_BRANCHES {
        eprintln!(
            "WARNING: branch count {} != baseline {} — accuracy/structure changed!",
            branches, EXPECTED_BRANCHES
        );
    }
}
