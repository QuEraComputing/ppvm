// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! End-to-end smoke test for the `cultivation_d5.stim` magic-state cultivation
//! circuit. It exercises the full extended pipeline added for this feature:
//! `RX`/`MX` X-basis prep & readout, native `T`/`T_DAG` rotations, and the
//! `MPP` multi-qubit Pauli-product detectors — alongside the existing
//! depolarizing / Pauli-error noise, resets, and CX/CZ Cliffords.

use std::path::PathBuf;

use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;
use ppvm_stim::{execute, parse_extended};
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<ByteFxHashF64<8>, usize>;

fn cultivation_src() -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("cultivation_d5.stim");
    std::fs::read_to_string(path).expect("read cultivation_d5.stim")
}

#[test]
fn cultivation_d5_parses_validates_and_runs() {
    let src = cultivation_src();
    let prog = parse_extended(&src).expect("cultivation_d5 must parse");

    // 36 Z-basis (M) + 57 X-basis (MX) single-qubit readouts, plus the two MPP
    // lines (1 + 18 Pauli-product detectors) = 112 recorded measurements.
    let mut tab: Tab = GeneralizedTableau::new_with_seed(64, 1e-10, 1);
    let results = execute(&prog, &mut tab).expect("cultivation_d5 must execute");
    assert_eq!(results.len(), 112);
    // Every measurement resolves to a concrete bit (no lost qubits in this file).
    assert!(results.iter().all(|r| r.is_some()));
}
