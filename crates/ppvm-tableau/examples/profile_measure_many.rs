// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Profiling driver for `measure_many` across coefficient-count regimes.
//!
//! `measure_many` has two cost components per measured qubit:
//!   - tableau work  — `compute_decomposition` (scans all stabilizers +
//!     destabilizers) and, on case-a, `update_tableau_according_to_outcome`.
//!     This is O(n_qubits) per qubit and *independent of the coefficient count*.
//!   - coefficient work — the case-a HashMap build / overlap / merge / normalize
//!     passes, which scale with the number of stored coefficients.
//!
//! The three workloads bracket the crossover between the two regimes. Every
//! state is a full-width entangled Clifford state (so measurements take the
//! case-a path), with `n_t` T-gates controlling the coefficient count:
//!   - `few`   — ~1 coefficient    (tableau-bound)
//!   - `mid`   — ~100 coefficients
//!   - `large` — ~1000 coefficients (coefficient-bound)
//!
//! Usage:
//!   # quick timing + achieved coefficient count (use this to calibrate n_t)
//!   cargo run -p ppvm-tableau --example profile_measure_many --release -- mid
//!
//!   # sustained hot loop (~FLAME_SECS seconds, default 6) for a sampling
//!   # profiler. scripts/profile_measure_many.sh drives this with samply
//!   # (no sudo on macOS) and summarizes via scripts/samply_top.py:
//!   samply record --save-only --unstable-presymbolicate -o p.json.gz -- \
//!     target/release/examples/profile_measure_many mid flame

use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;
use std::time::{Duration, Instant};

type Tab = GeneralizedTableau<Byte8F64<2>, u128>;

const N_QUBITS: usize = 85;

/// Number of T-gates per workload, tuned so the post-circuit coefficient count
/// lands near the target. Recalibrate with the quick-timing mode if the state
/// builder changes.
fn n_t_for(workload: &str) -> usize {
    match workload {
        "few" => 0,
        "mid" => 7,
        "large" => 10,
        other => panic!("unknown workload {other:?}; use few|mid|large"),
    }
}

/// Full-width entangled state: H on every qubit (Clifford superposition, so
/// measurements hit the case-a path), `n_t` T-gates to branch the coefficient
/// vector, then a CZ chain to entangle.
fn build_state(n_t: usize) -> Tab {
    let mut tab: Tab = GeneralizedTableau::new(N_QUBITS, 1e-10);
    for q in 0..N_QUBITS {
        tab.h(q);
    }
    for q in 0..n_t {
        tab.t(q);
    }
    for q in 0..N_QUBITS - 1 {
        tab.cz(q, q + 1);
    }
    tab
}

fn median(mut xs: Vec<Duration>) -> Duration {
    xs.sort();
    xs[xs.len() / 2]
}

fn main() {
    let mut args = std::env::args().skip(1);
    let workload = args.next().unwrap_or_else(|| "mid".to_string());
    let mode = args.next().unwrap_or_else(|| "quick".to_string());

    let n_t = n_t_for(&workload);
    let base = build_state(n_t);
    let n = base.n_qubits();
    let n_coeffs = base.coefficients.len();
    let all: Vec<usize> = (0..n).collect();

    // Setup info goes to stderr so it never pollutes a flamegraph's sampled window.
    eprintln!("workload={workload}  n_qubits={n}  n_t={n_t}  coefficients={n_coeffs}  mode={mode}");

    match mode.as_str() {
        "quick" => {
            let n_runs = 200;
            let mut shot = Vec::with_capacity(n_runs);
            for _ in 0..n_runs {
                let start = Instant::now();
                let mut t = base.fork(Some(42));
                std::hint::black_box(t.measure_many(&all));
                shot.push(start.elapsed());
            }
            let mut forks = Vec::with_capacity(n_runs);
            for _ in 0..n_runs {
                let start = Instant::now();
                let t = base.fork(Some(42));
                forks.push(start.elapsed());
                std::hint::black_box(t);
            }
            let shot_m = median(shot);
            let fork_m = median(forks);
            eprintln!("per-shot (fork + measure_many) median over {n_runs}: {shot_m:?}");
            eprintln!("    fork only:    {fork_m:?}");
            eprintln!("    measure_many: {:?}", shot_m.saturating_sub(fork_m));
            eprintln!(
                "    per qubit:    {:?}",
                shot_m.saturating_sub(fork_m) / n as u32
            );
        }
        "flame" => {
            // Size the hot loop to ~10s of sampling regardless of workload cost.
            let mut probe = Vec::with_capacity(50);
            for _ in 0..50 {
                let start = Instant::now();
                let mut t = base.fork(Some(42));
                std::hint::black_box(t.measure_many(&all));
                probe.push(start.elapsed());
            }
            let per = median(probe).max(Duration::from_nanos(1));
            let secs: u64 = std::env::var("FLAME_SECS")
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(6);
            let iters = (Duration::from_secs(secs).as_nanos() / per.as_nanos()).max(2_000) as u64;
            eprintln!("flame: ~{per:?}/shot -> {iters} iterations (~{secs}s)");
            // fork feeds measure_many and the result is black-boxed, so neither
            // call can be hoisted out of the loop.
            for _ in 0..iters {
                let mut t = std::hint::black_box(base.fork(Some(42)));
                std::hint::black_box(t.measure_many(&all));
            }
        }
        other => panic!("unknown mode {other:?}; use quick|flame"),
    }
}
