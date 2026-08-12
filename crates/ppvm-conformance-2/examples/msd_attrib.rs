// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Profiling harness for the MSD measurement sweep, both engines.

use ppvm_conformance_2::tableau::{Driver, MSD_QUBITS, NewWide, OldWide, msd_state};
use std::time::Instant;

fn main() {
    let side = std::env::args().nth(1).unwrap_or_else(|| "new".into());
    let reps: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(200);

    let start = Instant::now();
    if side == "old" {
        let base: OldWide = msd_state(Some(3));
        for _ in 0..reps {
            let mut t = base.fork(None);
            for i in 0..MSD_QUBITS {
                std::hint::black_box(t.measure(i));
            }
        }
    } else if side == "new_guard" {
        let base: NewWide = msd_state(Some(3));
        let mut scratch = ppvm_tableau_2::MeasureScratch::new();
        for _ in 0..reps {
            let mut t = base.fork(None);
            std::hint::black_box(t.measure_all_with_scratch(&mut scratch));
        }
    } else if side == "new_all" {
        let base: NewWide = msd_state(Some(3));
        for _ in 0..reps {
            let mut t = base.fork(None);
            std::hint::black_box(t.measure_all());
        }
    } else {
        let base: NewWide = msd_state(Some(3));
        for _ in 0..reps {
            let mut t = base.fork(None);
            for i in 0..MSD_QUBITS {
                std::hint::black_box(t.measure(i));
            }
        }
    }
    let elapsed = start.elapsed();
    println!(
        "{side}: {reps} sweeps in {elapsed:.2?} ({:.2?}/sweep, {:.1} ns/measurement)",
        elapsed / reps as u32,
        elapsed.as_nanos() as f64 / (reps * MSD_QUBITS) as f64
    );
}
