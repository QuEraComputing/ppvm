// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Reproducer for the storage-tier performance cliff observed from Python
//! benchmarks/ppvm/run_bench.py.
//!
//! For N in 33..=64, the standard Python dispatch routes to
//! `ByteFxHashF64<8>` (Hash3, 8-byte storage), but the **same circuit** on
//! the same N runs ~4-5x slower than the next-tier-up `ByteFxHashF64<16>`
//! (Hash4, 16-byte storage). The work is identical -- only the const-generic
//! storage size differs -- so this is a pure code-gen / data-layout effect.
//!
//! Standalone numbers from Python (J=1.0, n_steps=20, truncation 1e-6):
//! ```text
//!   N    storage          wall time
//!   32   Hash2 ([u8;4])    223 ms
//!   40   Hash3 ([u8;8])    544 ms
//!   48   Hash3 ([u8;8])   1044 ms
//!   56   Hash3 ([u8;8])   1842 ms     <- bottom of the cliff
//!   64   Hash4 ([u8;16])   401 ms     <- 4.6x faster on the next tier
//! ```
//! At N=56 both tiers can hold the circuit (Hash3 caps at 8*8=64 qubits, Hash4
//! at 8*16=128), so the same N can be run on both configs side by side. That's
//! what this binary does.
//!
//! Usage:
//! ```bash
//! # Time both back-to-back:
//! cargo run --release --example trotter_storage_cliff
//!
//! # Profile one tier at a time (set TIER env var):
//! TIER=8  cargo flamegraph --release --example trotter_storage_cliff
//! TIER=16 cargo flamegraph --release --example trotter_storage_cliff
//!
//! # Tweakables:
//! N_QUBITS=56 ITERS=5 TIER=8 cargo run --release --example trotter_storage_cliff
//! ```

use std::time::Instant;

use ppvm_pauli_sum::prelude::*;
use ppvm_pauli_sum::strategy::CoefficientThreshold;

const H: f64 = 1.0;
const DT: f64 = 0.1 / H;
const TIME: f64 = 2.0 / H;
const J: f64 = 1.0 * H;
const MIN_ABS_COEFF: f64 = 1e-6;
const NOISE: [f64; 3] = [1e-4 / 4.0; 3];

macro_rules! run_tier {
    ($cfg:ty, $n:expr, $iters:expr, $label:expr) => {{
        let n = $n;
        let strat = CoefficientThreshold(MIN_ABS_COEFF);
        let mut seed: PauliSum<$cfg> = PauliSum::builder()
            .n_qubits(n)
            .strategy(strat)
            .capacity(n.pow(2))
            .build();
        for i in 0..n {
            let term: String = (0..n).map(|j| if j == i { 'Z' } else { 'I' }).collect();
            seed += (term.as_str(), 1.0);
        }

        let steps = (TIME / DT) as usize;
        let theta_zz = DT * J;
        let theta_x = DT * H;

        let mut times_ms = Vec::with_capacity($iters);
        let mut final_len = 0_usize;
        for _ in 0..$iters {
            let mut state = seed.clone();
            let t0 = Instant::now();
            for _ in 0..steps {
                for i in 0..n {
                    state.rx(i, theta_x);
                    state.truncate();
                    state.pauli_error(i, NOISE);
                    state.truncate();
                }
                for i in 0..n - 1 {
                    state.rzz(i, i + 1, theta_zz);
                    state.truncate();
                    state.pauli_error(i, NOISE);
                    state.truncate();
                    state.pauli_error(i + 1, NOISE);
                    state.truncate();
                }
            }
            times_ms.push(t0.elapsed().as_secs_f64() * 1000.0);
            final_len = state.len();
        }
        let min = times_ms.iter().cloned().fold(f64::INFINITY, f64::min);
        let mean = times_ms.iter().sum::<f64>() / times_ms.len() as f64;
        println!(
            "{:>20}  N={:3}  iters={:2}  min={:8.1} ms  mean={:8.1} ms  |state|={}",
            $label, n, $iters, min, mean, final_len
        );
    }};
}

fn main() {
    let n: usize = std::env::var("N_QUBITS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(56);
    let iters: usize = std::env::var("ITERS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);
    let tier = std::env::var("TIER").ok();

    assert!(n <= 64, "Hash3 (8-byte) tier can only hold up to 64 qubits");

    println!(
        "Trotter storage-cliff reproducer: N={}, iters={}, J={}, h={}, dt={}, n_steps={}",
        n,
        iters,
        J,
        H,
        DT,
        (TIME / DT) as usize
    );

    type Cfg4 = config::indexmap::ByteFxHashF64<4, CoefficientThreshold>;
    type Cfg8 = config::indexmap::ByteFxHashF64<8, CoefficientThreshold>;
    type Cfg16 = config::indexmap::ByteFxHashF64<16, CoefficientThreshold>;
    type Cfg32 = config::indexmap::ByteFxHashF64<32, CoefficientThreshold>;

    match tier.as_deref() {
        Some("4") if n <= 32 => run_tier!(Cfg4, n, iters, "ByteFxHashF64<4>"),
        Some("8") => run_tier!(Cfg8, n, iters, "ByteFxHashF64<8>"),
        Some("16") => run_tier!(Cfg16, n, iters, "ByteFxHashF64<16>"),
        Some("32") => run_tier!(Cfg32, n, iters, "ByteFxHashF64<32>"),
        _ => {
            if n <= 32 {
                run_tier!(Cfg4, n, iters, "ByteFxHashF64<4>");
            }
            run_tier!(Cfg8, n, iters, "ByteFxHashF64<8>");
            run_tier!(Cfg16, n, iters, "ByteFxHashF64<16>");
            run_tier!(Cfg32, n, iters, "ByteFxHashF64<32>");
        }
    }
}
