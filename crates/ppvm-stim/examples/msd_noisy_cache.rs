// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Noisy MSD-shaped STIM sampling with and without trajectory caching.
//!
//! This is intentionally STIM rather than `.sst`: the circuit is explicit,
//! compact to generate, and it exercises the same shared trajectory-cache path
//! used by bytecode execution. Run for example:
//!
//! ```text
//! cargo run --release -p ppvm-stim --example msd_noisy_cache -- 512 32768
//! ```
//!
//! The first argument is the shot count. The second is the maximum number of
//! cached continuation states.

use std::time::{Duration, Instant};

use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;
use ppvm_stim::{parse_extended, sample_cached, sample_serial};
use ppvm_tableau::prelude::*;
use ppvm_trajectory_cache::{CacheConfig, CacheStats};

type Tab = GeneralizedTableau<ByteFxHashF64<11>, u128>;

const QUBITS_PER_CODE_BLOCK: usize = 17;
const N_BLOCKS: usize = 5;
const N_QUBITS: usize = QUBITS_PER_CODE_BLOCK * N_BLOCKS;
const DEFAULT_SHOTS: usize = 256;
const DEFAULT_CACHE_NODES: usize = 16_384;
const BASE_SEED: u64 = 0x5eed;

// Low-probability noise gives the cache a useful hot path while still making
// the run genuinely stochastic. These mirror the scale used by msd-noisy.
const P_LOSS: f64 = 1e-4;
const P_DEPOLARIZE: f64 = 1e-4;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let shots = parse_arg(1, DEFAULT_SHOTS, "shots")?;
    let cache_nodes = parse_arg(2, DEFAULT_CACHE_NODES, "cache nodes")?;

    let circuit = msd_noisy_stim_string();
    let prog = parse_extended(&circuit)?;

    println!(
        "noisy MSD STIM: {N_QUBITS} qubits, {} measurements, {shots} shots",
        prog.measurement_count()
    );
    println!("noise: loss={P_LOSS}, depolarize={P_DEPOLARIZE}");

    let start = Instant::now();
    let uncached = sample_serial::<_, _, _, _>(&prog, shots, seeded_tableau)?;
    let uncached_elapsed = start.elapsed();
    print_run("uncached", uncached_elapsed, &uncached, None);

    let start = Instant::now();
    let cached = sample_cached::<_, _, _, _>(
        &prog,
        shots,
        seeded_tableau,
        CacheConfig::bounded(cache_nodes),
    )?;
    let cached_elapsed = start.elapsed();
    print_run(
        "cached",
        cached_elapsed,
        &cached.output,
        Some(cached.cache_stats),
    );

    if cached_elapsed.as_nanos() > 0 {
        println!(
            "speedup: {:.2}x",
            uncached_elapsed.as_secs_f64() / cached_elapsed.as_secs_f64()
        );
    }

    Ok(())
}

fn parse_arg(index: usize, default: usize, label: &str) -> Result<usize, String> {
    match std::env::args().nth(index) {
        Some(raw) => raw
            .parse()
            .map_err(|_| format!("invalid {label} argument: {raw}")),
        None => Ok(default),
    }
}

fn seeded_tableau(shot: usize) -> Tab {
    Tab::new_with_seed(N_QUBITS, 1e-10, BASE_SEED.wrapping_add(shot as u64))
}

fn print_run(
    label: &str,
    elapsed: Duration,
    shots: &[Vec<Option<bool>>],
    cache_stats: Option<CacheStats>,
) {
    let (zeros, ones, lost) = shot_counts(shots);
    println!(
        "{label}: {elapsed:.2?} ({:.2?}/shot), zeros={zeros}, ones={ones}, lost={lost}",
        elapsed / shots.len().max(1) as u32
    );
    if let Some(stats) = cache_stats {
        println!(
            "cache: states={} hits={} misses={} evictions={} terminal_hits={}",
            stats.nodes, stats.hits, stats.misses, stats.evictions, stats.terminal_hits
        );
    }
}

fn shot_counts(shots: &[Vec<Option<bool>>]) -> (usize, usize, usize) {
    let mut zeros = 0;
    let mut ones = 0;
    let mut lost = 0;
    for outcome in shots.iter().flatten() {
        match outcome {
            Some(false) => zeros += 1,
            Some(true) => ones += 1,
            None => lost += 1,
        }
    }
    (zeros, ones, lost)
}

fn msd_noisy_stim_string() -> String {
    let qubit_addrs: Vec<usize> = (0..N_QUBITS).collect();
    let ql: Vec<&[usize]> = qubit_addrs.chunks_exact(QUBITS_PER_CODE_BLOCK).collect();
    let all_qubits: Vec<usize> = (0..N_QUBITS).collect();

    let mut lines = Vec::new();

    for q in &ql {
        let encoding_qubit = q[7];
        push_1q_gate(&mut lines, "H", &[encoding_qubit]);
        push_1q_gate(&mut lines, "S[T]", &[encoding_qubit]);
        encode_stim(&mut lines, q);
    }

    for i in [0, 1, 4] {
        push_1q_gate(&mut lines, "SQRT_X", ql[i]);
    }

    push_cz_pairs(&mut lines, ql[0], ql[1]);
    push_cz_pairs(&mut lines, ql[2], ql[3]);

    push_1q_gate(&mut lines, "SQRT_Y", ql[0]);
    push_1q_gate(&mut lines, "SQRT_Y", ql[3]);

    push_cz_pairs(&mut lines, ql[0], ql[2]);
    push_cz_pairs(&mut lines, ql[3], ql[4]);

    push_1q_gate(&mut lines, "SQRT_X_DAG", ql[0]);

    push_cz_pairs(&mut lines, ql[0], ql[4]);
    push_cz_pairs(&mut lines, ql[1], ql[3]);

    for q in ql.iter().take(N_BLOCKS) {
        push_1q_gate(&mut lines, "SQRT_X_DAG", q);
    }

    lines.push(fmt_gate("M", &all_qubits));
    lines.join("\n")
}

fn encode_stim(lines: &mut Vec<String>, q: &[usize]) {
    assert_eq!(q.len(), QUBITS_PER_CODE_BLOCK);

    push_1q_gate_indices(
        lines,
        "SQRT_Y",
        q,
        &[0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16],
    );

    push_cz_index_pairs(lines, q, &[[1, 3], [7, 10], [12, 14], [13, 16]]);
    push_1q_gate_indices(lines, "SQRT_Y_DAG", q, &[7, 16]);

    push_cz_index_pairs(lines, q, &[[4, 7], [8, 10], [11, 14], [15, 16]]);
    push_1q_gate_indices(lines, "SQRT_Y_DAG", q, &[4, 10, 14, 16]);

    push_cz_index_pairs(lines, q, &[[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]]);
    push_1q_gate_indices(lines, "SQRT_Y", q, &[3, 6, 9, 10, 12, 13]);

    push_cz_index_pairs(lines, q, &[[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]]);
    push_1q_gate_indices(lines, "SQRT_Y", q, &[1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14]);

    push_cz_index_pairs(
        lines,
        q,
        &[[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]],
    );
    push_1q_gate_indices(lines, "SQRT_Y_DAG", q, &[0, 2, 5, 6, 8, 10, 12]);
}

fn push_1q_gate(lines: &mut Vec<String>, gate: &str, qubits: &[usize]) {
    lines.push(fmt_gate(gate, qubits));
    push_1q_noise(lines, qubits);
}

fn push_1q_gate_indices(lines: &mut Vec<String>, gate: &str, q: &[usize], indices: &[usize]) {
    let targets: Vec<usize> = indices.iter().map(|&i| q[i]).collect();
    push_1q_gate(lines, gate, &targets);
}

fn push_cz_pairs(lines: &mut Vec<String>, controls: &[usize], targets: &[usize]) {
    let flat: Vec<usize> = controls
        .iter()
        .zip(targets)
        .flat_map(|(&c, &t)| [c, t])
        .collect();
    lines.push(fmt_gate("CZ", &flat));
    push_2q_noise(lines, &flat);
}

fn push_cz_index_pairs(lines: &mut Vec<String>, q: &[usize], pairs: &[[usize; 2]]) {
    let flat: Vec<usize> = pairs.iter().flat_map(|[i, j]| [q[*i], q[*j]]).collect();
    lines.push(fmt_gate("CZ", &flat));
    push_2q_noise(lines, &flat);
}

fn push_1q_noise(lines: &mut Vec<String>, qubits: &[usize]) {
    lines.push(fmt_arg_gate("I_ERROR[loss]", P_LOSS, qubits));
    lines.push(fmt_arg_gate("DEPOLARIZE1", P_DEPOLARIZE, qubits));
}

fn push_2q_noise(lines: &mut Vec<String>, flat_pairs: &[usize]) {
    lines.push(fmt_arg_gate("I_ERROR[loss]", P_LOSS, flat_pairs));
    lines.push(fmt_arg_gate("DEPOLARIZE2", P_DEPOLARIZE, flat_pairs));
}

fn fmt_gate(gate: &str, qubits: &[usize]) -> String {
    let targets: Vec<String> = qubits.iter().map(|q| q.to_string()).collect();
    format!("{gate} {}", targets.join(" "))
}

fn fmt_arg_gate(gate: &str, arg: f64, qubits: &[usize]) -> String {
    let targets: Vec<String> = qubits.iter().map(|q| q.to_string()).collect();
    format!("{gate}({arg}) {}", targets.join(" "))
}
