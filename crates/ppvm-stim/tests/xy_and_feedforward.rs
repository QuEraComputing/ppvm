// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! End-to-end semantics for the X/Y-basis measurement & reset gates
//! (`MX`/`MY`/`MRX`/`MRY`, `RX`/`RY`) and measurement-record controlled
//! feed-forward gates (`CX`/`CY`/`CZ rec[-k]`).

use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;
use ppvm_stim::{execute, parse_extended};
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<ByteFxHashF64<1>, usize>;

fn run(src: &str, qubits: usize) -> Vec<Option<bool>> {
    let prog = parse_extended(src).expect("parse_extended");
    let mut t: Tab = GeneralizedTableau::new(qubits, 1e-12);
    execute(&prog, &mut t).expect("execute")
}

// ---- X/Y-basis measurement & reset -----------------------------------------

#[test]
fn rx_prepares_plus_state_measured_zero_in_x() {
    // RX -> |+>, which is the +1 eigenstate of X, so MX reads 0.
    assert_eq!(run("RX 0\nMX 0", 1), vec![Some(false)]);
}

#[test]
fn z_flips_plus_to_minus_read_one_in_x() {
    // Z|+> = |->, the -1 eigenstate of X, so MX reads 1.
    assert_eq!(run("RX 0\nZ 0\nMX 0", 1), vec![Some(true)]);
}

#[test]
fn ry_prepares_i_state_measured_zero_in_y() {
    // RY -> |i>, the +1 eigenstate of Y, so MY reads 0.
    assert_eq!(run("RY 0\nMY 0", 1), vec![Some(false)]);
}

#[test]
fn z_flips_i_to_minus_i_read_one_in_y() {
    // Z|i> = |-i>, the -1 eigenstate of Y, so MY reads 1.
    assert_eq!(run("RY 0\nZ 0\nMY 0", 1), vec![Some(true)]);
}

#[test]
fn mx_leaves_qubit_in_measured_eigenstate() {
    // After MX the qubit stays in the X eigenstate it was projected onto, so a
    // second MX agrees with the first.
    assert_eq!(run("H 0\nMX 0\nMX 0", 1), vec![Some(false), Some(false)]);
}

#[test]
fn mrx_resets_to_plus_regardless_of_outcome() {
    // Start in |-> (read 1), MRX records 1 and resets to |+>, so the following
    // MX reads 0.
    assert_eq!(
        run("RX 0\nZ 0\nMRX 0\nMX 0", 1),
        vec![Some(true), Some(false)]
    );
}

#[test]
fn mry_resets_to_i_regardless_of_outcome() {
    assert_eq!(
        run("RY 0\nZ 0\nMRY 0\nMY 0", 1),
        vec![Some(true), Some(false)]
    );
}

#[test]
fn mx_readout_noise_always_flips() {
    // p=1 readout flip turns the deterministic 0 into a reported 1.
    assert_eq!(run("RX 0\nMX(1.0) 0", 1), vec![Some(true)]);
}

// ---- measurement-record controlled feed-forward ----------------------------

#[test]
fn cx_rec_applies_x_when_control_bit_is_one() {
    // q0 = |1> -> measured 1 -> CX rec[-1] 1 applies X to q1 -> M1 reads 1.
    assert_eq!(
        run("X 0\nM 0\nCX rec[-1] 1\nM 1", 2),
        vec![Some(true), Some(true)]
    );
}

#[test]
fn cx_rec_is_noop_when_control_bit_is_zero() {
    // q0 = |0> -> measured 0 -> CX rec[-1] 1 does nothing -> M1 reads 0.
    assert_eq!(
        run("M 0\nCX rec[-1] 1\nM 1", 2),
        vec![Some(false), Some(false)]
    );
}

#[test]
fn cy_rec_applies_y_when_control_bit_is_one() {
    // Y|0> = i|1>, so M1 reads 1.
    assert_eq!(
        run("X 0\nM 0\nCY rec[-1] 1\nM 1", 2),
        vec![Some(true), Some(true)]
    );
}

#[test]
fn cz_rec_applies_z_when_control_bit_is_one() {
    // Z on |+> gives |->; H maps it back to |1>, so M1 reads 1.
    assert_eq!(
        run("X 0\nM 0\nH 1\nCZ rec[-1] 1\nH 1\nM 1", 2),
        vec![Some(true), Some(true)]
    );
}

#[test]
fn deeper_lookback_selects_the_right_bit() {
    // Record before the CX: [1 (q0), 0 (q1)]. rec[-1] is q1's bit (0), rec[-2]
    // is q0's bit (1) -> the CX rec[-2] applies X to q2, so M2 reads 1.
    assert_eq!(
        run("X 0\nM 0\nM 1\nCX rec[-2] 2\nM 2", 3),
        vec![Some(true), Some(false), Some(true)]
    );
}

#[test]
fn classic_feedforward_correction_in_repeat() {
    // rec controls resolve against the live record even inside REPEAT bodies.
    let out = run("X 0\nREPEAT 1 {\n    M 0\n    CX rec[-1] 1\n}\nM 1", 2);
    assert_eq!(out, vec![Some(true), Some(true)]);
}
