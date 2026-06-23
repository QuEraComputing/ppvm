// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! End-to-end semantics for multi-qubit Pauli-product measurements (`MPP`),
//! implemented via the basis-change + CX-ladder gadget.

use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;
use ppvm_stim::{execute, parse_extended};
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<ByteFxHashF64<1>, usize>;

fn run(src: &str, qubits: usize) -> Vec<Option<bool>> {
    let prog = parse_extended(src).expect("parse_extended");
    let mut t: Tab = GeneralizedTableau::new(qubits, 1e-12);
    execute(&prog, &mut t).expect("execute")
}

// ---- single-factor products (reduce to MX/MY/MZ) ---------------------------

#[test]
fn mpp_single_z_reads_computational_basis() {
    assert_eq!(run("MPP Z0", 1), vec![Some(false)]);
    assert_eq!(run("X 0\nMPP Z0", 1), vec![Some(true)]);
}

#[test]
fn mpp_single_x_reads_plus_minus() {
    // |+> is +1 of X, Z|+> = |-> is -1 of X.
    assert_eq!(run("H 0\nMPP X0", 1), vec![Some(false)]);
    assert_eq!(run("H 0\nZ 0\nMPP X0", 1), vec![Some(true)]);
}

#[test]
fn mpp_single_y_reads_i_eigenstates() {
    // RY -> |i>, the +1 eigenstate of Y; Z flips it to -1.
    assert_eq!(run("RY 0\nMPP Y0", 1), vec![Some(false)]);
    assert_eq!(run("RY 0\nZ 0\nMPP Y0", 1), vec![Some(true)]);
}

// ---- two-qubit parity products ---------------------------------------------

#[test]
fn mpp_zz_measures_computational_parity() {
    // Even parity -> +1 -> 0; odd parity -> -1 -> 1.
    assert_eq!(run("MPP Z0*Z1", 2), vec![Some(false)]);
    assert_eq!(run("X 0\nMPP Z0*Z1", 2), vec![Some(true)]);
    assert_eq!(run("X 0\nX 1\nMPP Z0*Z1", 2), vec![Some(false)]);
}

#[test]
fn mpp_xx_on_bell_state_is_plus_one() {
    // Bell state (|00>+|11>)/sqrt2 is a +1 eigenstate of both X0*X1 and Z0*Z1.
    assert_eq!(run("H 0\nCX 0 1\nMPP X0*X1", 2), vec![Some(false)]);
    assert_eq!(run("H 0\nCX 0 1\nMPP Z0*Z1", 2), vec![Some(false)]);
}

#[test]
fn mpp_xx_on_odd_bell_state_is_minus_one() {
    // Z on the Bell state gives (|00>-|11>)/sqrt2, the -1 eigenstate of X0*X1.
    assert_eq!(run("H 0\nCX 0 1\nZ 0\nMPP X0*X1", 2), vec![Some(true)]);
}

// ---- multiple products on one line -----------------------------------------

#[test]
fn mpp_multiple_products_yield_one_result_each() {
    // Two space-separated products -> two results.
    assert_eq!(run("X 0\nMPP Z0 Z1", 2), vec![Some(true), Some(false)]);
}

#[test]
fn mpp_three_qubit_product_parity() {
    // GHZ (|000>+|111>)/sqrt2: Z0*Z1*Z2 has even parity on both branches -> +1.
    assert_eq!(
        run("H 0\nCX 0 1\nCX 0 2\nMPP Z0*Z1*Z2", 3),
        vec![Some(false)]
    );
}

// ---- non-destructive: repeated measurement agrees --------------------------

#[test]
fn mpp_is_non_destructive() {
    // Measuring the same product twice on the same state gives the same result
    // and leaves the computational-basis readout consistent.
    assert_eq!(
        run("X 0\nMPP Z0*Z1\nMPP Z0*Z1\nM 0\nM 1", 2),
        vec![Some(true), Some(true), Some(true), Some(false)]
    );
}

// ---- mixed-basis product ---------------------------------------------------

#[test]
fn mpp_mixed_xz_product() {
    // Prepare |+>_0 |0>_1: X0 = +1, Z1 = +1, so X0*Z1 = +1 -> 0.
    assert_eq!(run("H 0\nMPP X0*Z1", 2), vec![Some(false)]);
    // Flip qubit 1: Z1 = -1, so X0*Z1 = -1 -> 1.
    assert_eq!(run("H 0\nX 1\nMPP X0*Z1", 2), vec![Some(true)]);
}
