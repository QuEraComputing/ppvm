// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Tight-loop driver for `cargo flamegraph` of the per-shot sampling cost
//! (fork + measure_all) on the 85-qubit MSD state.
//!
//! Run:
//!   sudo cargo flamegraph --example profile_measure_all_flame -p ppvm-tableau --release
//!
//! For function names (not raw addresses) to appear, the workspace
//! `[profile.release]` should have `debug = "line-tables-only"` (or
//! `debug = true`). cargo-flamegraph tries to set this for you via env vars,
//! but if you see hex addresses, add it manually to /Users/david/git/ppvm/Cargo.toml:
//!
//!     [profile.release]
//!     debug = "line-tables-only"
//!
//! 200k iterations × ~40 µs/shot ≈ 8 s of dtrace sampling. The setup runs
//! once, so it'll be a thin slice of the flamegraph that's easy to ignore.
//! Look for the `measure_all` subtree to read off:
//!   - `compute_decomposition` cost
//!   - `update_tableau_according_to_outcome` cost
//!   - HashMap traffic in the case-a path
//!
//! And a small adjacent subtree for `fork` (expect ~1% based on the
//! instrumented run).
//!
//! NOTE: the MSD setup is duplicated from `profile_measure_all.rs` — fine for
//! examples, hoist into a shared module if we keep growing these.

use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::measure_all::LossyMeasureAll;
use ppvm_tableau::prelude::*;

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
        tab.sqrt_x_dag(*q);
    }
    for (control, target) in ql[0].iter().zip(ql[4]) {
        tab.cz(*control, *target);
    }
    for (control, target) in ql[1].iter().zip(ql[3]) {
        tab.cz(*control, *target);
    }
    for block in ql.iter().take(5) {
        for q in *block {
            tab.sqrt_x_dag(*q);
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
        tab.sqrt_y_dag(qubits[i]);
    }
    for [i, j] in [[4, 7], [8, 10], [11, 14], [15, 16]] {
        tab.cz(qubits[i], qubits[j]);
    }
    for i in [4, 10, 14, 16] {
        tab.sqrt_y_dag(qubits[i]);
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
        tab.sqrt_y_dag(qubits[i]);
    }
}

fn main() {
    let base = build_msd_state();

    // Warm-up: lets caches / branch predictors settle so the recorded window
    // reflects steady-state cost.
    for _ in 0..2_000 {
        let mut t = base.fork(Some(42));
        std::hint::black_box(t.measure_all());
    }

    // Hot loop. The compiler can't hoist either call out because fork's
    // result feeds measure_all and measure_all's result is fed to black_box.
    for _ in 0..200_000 {
        let mut t = std::hint::black_box(base.fork(Some(42)));
        std::hint::black_box(t.measure_all());
    }
}
