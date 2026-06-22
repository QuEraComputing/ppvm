// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::time::Instant;

use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;
use ppvm_tableau_sum::data::GeneralizedTableauSum;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

type GTabSum = GeneralizedTableauSum<Byte8F64<2>, u128>;

fn encode(tab: &mut GTabSum, qubits: &[usize], p_loss: f64, p_depolarize: f64) {
    if qubits.len() == 17 {
        for i in [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16] {
            tab.sqrt_y(qubits[i]);
            tab.loss_channel(qubits[i], p_loss);
            tab.depolarize(qubits[i], p_depolarize);
        }
        println!("Branches: {}", tab.len());

        for [i, j] in [[1, 3], [7, 10], [12, 14], [13, 16]] {
            tab.cz(qubits[i], qubits[j]);
            tab.loss_channel(qubits[i], p_loss);
            tab.loss_channel(qubits[j], p_loss);
            tab.depolarize(qubits[i], p_depolarize);
            tab.depolarize(qubits[j], p_depolarize);
        }
        for i in [7, 16] {
            tab.sqrt_y_adj(qubits[i]);
            tab.loss_channel(qubits[i], p_loss);
            tab.depolarize(qubits[i], p_depolarize);
        }
        for [i, j] in [[4, 7], [8, 10], [11, 14], [15, 16]] {
            tab.cz(qubits[i], qubits[j]);
            tab.loss_channel(qubits[i], p_loss);
            tab.loss_channel(qubits[j], p_loss);
            tab.depolarize(qubits[i], p_depolarize);
            tab.depolarize(qubits[j], p_depolarize);
        }
        for i in [4, 10, 14, 16] {
            tab.sqrt_y_adj(qubits[i]);
            tab.loss_channel(qubits[i], p_loss);
            tab.depolarize(qubits[i], p_depolarize);
        }
        for [i, j] in [[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]] {
            tab.cz(qubits[i], qubits[j]);
            tab.loss_channel(qubits[i], p_loss);
            tab.loss_channel(qubits[j], p_loss);
            tab.depolarize(qubits[i], p_depolarize);
            tab.depolarize(qubits[j], p_depolarize);
        }
        for i in [3, 6, 9, 10, 12, 13] {
            tab.sqrt_y(qubits[i]);
            tab.loss_channel(qubits[i], p_loss);
            tab.depolarize(qubits[i], p_depolarize);
        }
        for [i, j] in [[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]] {
            tab.cz(qubits[i], qubits[j]);
            tab.loss_channel(qubits[i], p_loss);
            tab.loss_channel(qubits[j], p_loss);
            tab.depolarize(qubits[i], p_depolarize);
            tab.depolarize(qubits[j], p_depolarize);
        }
        for i in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14] {
            tab.sqrt_y(qubits[i]);
            tab.loss_channel(qubits[i], p_loss);
            tab.depolarize(qubits[i], p_depolarize);
        }
        for [i, j] in [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]] {
            tab.cz(qubits[i], qubits[j]);
            tab.loss_channel(qubits[i], p_loss);
            tab.loss_channel(qubits[j], p_loss);
            tab.depolarize(qubits[i], p_depolarize);
            tab.depolarize(qubits[j], p_depolarize);
        }
        for i in [0, 2, 5, 6, 8, 10, 12] {
            tab.sqrt_y_adj(qubits[i]);
            tab.loss_channel(qubits[i], p_loss);
            tab.depolarize(qubits[i], p_depolarize);
        }
    }
}

fn main() {
    let n_qubits = 85;
    let p_loss = 1e-4;
    let p_depolarize = 1e-4;
    let sum_cutoff = 1e-7;
    let n_shots = 1000;
    let build_start = Instant::now();

    let mut tab: GTabSum = GeneralizedTableauSum::new(n_qubits, 1e-10, sum_cutoff);
    let qubit_addrs: Vec<usize> = (0..n_qubits).collect();
    let ql: Vec<&[usize]> = qubit_addrs.chunks_exact(17).collect();

    // Phase 1: Encoding (H + T + encode per block)
    for q in ql.iter() {
        let encoding_qubit = q[7];

        tab.h(encoding_qubit);
        tab.loss_channel(encoding_qubit, p_loss);
        tab.depolarize(encoding_qubit, p_depolarize);

        tab.t(encoding_qubit);
        tab.loss_channel(encoding_qubit, p_loss);
        tab.depolarize(encoding_qubit, p_depolarize);

        encode(&mut tab, q, p_loss, p_depolarize);
    }

    println!("Branches: {}", tab.len());

    // Phase 2: Middle gates (sqrt_x, cz, sqrt_y, sqrt_x_adj layers)
    for i in [0, 1, 4] {
        for q in ql[i] {
            tab.sqrt_x(*q);
            tab.loss_channel(*q, p_loss);
            tab.depolarize(*q, p_depolarize);
        }
    }

    println!("Branches: {}", tab.len());

    for (control, target) in ql[0].iter().zip(ql[1]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize(*control, p_depolarize);
        tab.depolarize(*target, p_depolarize);
    }
    for (control, target) in ql[2].iter().zip(ql[3]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize(*control, p_depolarize);
        tab.depolarize(*target, p_depolarize);
    }
    for q in ql[0] {
        tab.sqrt_y(*q);
        tab.loss_channel(*q, p_loss);
        tab.depolarize(*q, p_depolarize);
    }
    for q in ql[3] {
        tab.sqrt_y(*q);
        tab.loss_channel(*q, p_loss);
        tab.depolarize(*q, p_depolarize);
    }
    for (control, target) in ql[0].iter().zip(ql[2]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize(*control, p_depolarize);
        tab.depolarize(*target, p_depolarize);
    }
    for (control, target) in ql[3].iter().zip(ql[4]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize(*control, p_depolarize);
        tab.depolarize(*target, p_depolarize);
    }
    for q in ql[0] {
        tab.sqrt_x_adj(*q);
        tab.loss_channel(*q, p_loss);
        tab.depolarize(*q, p_depolarize);
    }
    for (control, target) in ql[0].iter().zip(ql[4]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize(*control, p_depolarize);
        tab.depolarize(*target, p_depolarize);
    }
    for (control, target) in ql[1].iter().zip(ql[3]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize(*control, p_depolarize);
        tab.depolarize(*target, p_depolarize);
    }
    for block in ql.iter().take(5) {
        for q in *block {
            tab.sqrt_x_adj(*q);
            tab.loss_channel(*q, p_loss);
            tab.depolarize(*q, p_depolarize);
        }
    }

    println!("Branches: {}", tab.len());
    let total_build_time = Instant::now() - build_start;
    println!("Build time: {} ms", total_build_time.as_millis());

    let mut sampler = tab.sampler();
    let now = Instant::now();
    sampler.sample_shots(n_shots);
    let sample_time = Instant::now() - now;
    let per_shot_us = sample_time.as_nanos() / n_shots as u128;
    println!(
        "Time to {} samples: {} us",
        n_shots,
        sample_time.as_micros()
    );
    println!("Time per shot: {} ns", per_shot_us);
}
