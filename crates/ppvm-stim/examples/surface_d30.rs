// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Run the distance-30 surface-code memory circuit in `surface_d30.stim`
//! against a [`GeneralizedTableau`].
//!
//! The circuit is pure Clifford (resets, H, CX, measurements, plus
//! depolarizing / X noise), so no non-Clifford branching occurs and a plain
//! `usize` branch index suffices. Its 1889 qubits need `ceil(1889 / 8) = 237`
//! bytes per Pauli word, hence the `<237>` in the config.

use std::time::Instant;

use bnum::types::U2048;
use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;
use ppvm_stim::{parse_extended, sample};
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<ByteFxHashF64<237>, U2048>;

const N_QUBITS: usize = 1889;
const SHOTS: usize = 12;

fn main() -> Result<(), ppvm_stim::Error> {
    // The .stim file sits next to this example; embed it at compile time.
    let prog = parse_extended(include_str!("surface_d30.stim"))?;
    println!(
        "parsed {} measurement records over {N_QUBITS} qubits",
        prog.measurement_count()
    );

    // Parse once, build a fresh tableau per shot via the factory closure.
    let start = Instant::now();
    let shots = sample(&prog, SHOTS, |_| Tab::new(N_QUBITS, 1e-10))?;
    let elapsed = start.elapsed();

    for (i, shot) in shots.iter().enumerate() {
        let ones = shot.iter().filter(|&&b| b == Some(true)).count();
        println!("shot {i}: {ones} ones / {} measurements", shot.len());
    }
    println!(
        "{SHOTS} shots in {elapsed:.2?} ({:.2?}/shot)",
        elapsed / SHOTS as u32
    );
    Ok(())
}
