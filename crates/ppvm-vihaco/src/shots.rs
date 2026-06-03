// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Running a compiled program for many shots, optionally across threads.
//!
//! Each shot runs on a fresh [`PPVM`] so shots are fully independent; the
//! module is compiled once and shared. With the `rayon` feature, [`run_shots`]
//! parallelizes across shots when asked for more than one thread and there are
//! enough shots to amortize the overhead.

use crate::PPVMModule;
use crate::composite::PPVM;
use crate::measurements::MeasurementResult;

/// Below this many shots, parallelism's overhead outweighs its benefit and we
/// always run serially. Provisional — tune with benchmarks.
pub const PARALLEL_SHOT_THRESHOLD: usize = 128;

/// Per-shot seed derived from the base seed and the shot index, so every shot
/// gets a distinct RNG stream (a shared seed would make all shots identical).
/// Depends only on `(base, index)`, so serial and parallel runs are bit-for-bit
/// identical for a given seed regardless of thread count.
fn shot_seed(base: Option<u64>, index: usize) -> Option<u64> {
    base.map(|b| b.wrapping_add(index as u64))
}

/// Run a single shot on a fresh machine and return its measurement record.
fn run_one_shot(module: &PPVMModule, seed: Option<u64>) -> eyre::Result<Vec<MeasurementResult>> {
    let mut machine = PPVM::default();
    machine.load(module)?;
    machine.run_with_seed(seed)?;
    Ok(machine.measurement_record())
}

/// Run `shots` shots serially. One entry per shot, in order.
pub fn run_shots_serial(
    module: &PPVMModule,
    shots: usize,
    seed: Option<u64>,
) -> eyre::Result<Vec<Vec<MeasurementResult>>> {
    (0..shots)
        .map(|i| run_one_shot(module, shot_seed(seed, i)))
        .collect()
}

/// Run `shots` shots across a scoped rayon pool of `threads` threads. One entry
/// per shot, in order (preserved by the indexed parallel iterator).
#[cfg(feature = "rayon")]
pub fn run_shots_parallel(
    module: &PPVMModule,
    shots: usize,
    threads: usize,
    seed: Option<u64>,
) -> eyre::Result<Vec<Vec<MeasurementResult>>> {
    use rayon::prelude::*;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;
    pool.install(|| {
        (0..shots)
            .into_par_iter()
            .map(|i| run_one_shot(module, shot_seed(seed, i)))
            .collect()
    })
}

/// Run `shots` shots, choosing serial or parallel execution. Goes parallel only
/// when built with `rayon`, more than one thread is requested, and there are
/// enough shots to be worth it; otherwise runs serially.
pub fn run_shots(
    module: &PPVMModule,
    shots: usize,
    threads: usize,
    seed: Option<u64>,
) -> eyre::Result<Vec<Vec<MeasurementResult>>> {
    #[cfg(feature = "rayon")]
    if threads > 1 && shots >= PARALLEL_SHOT_THRESHOLD {
        return run_shots_parallel(module, shots, threads, seed);
    }
    #[cfg(not(feature = "rayon"))]
    let _ = threads;

    run_shots_serial(module, shots, seed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compile_program;
    use crate::measurements::MeasurementOutcome;

    /// Measures q0 in |0>: every shot is deterministically `0`.
    const DETERMINISTIC: &str =
        "device circuit.n_qubits 1;\nfn @main() { const.u64 0\n gate measure\n ret }\n";

    /// Prepares |+> with H, then measures q0: each shot is a random 0/1.
    const RANDOM: &str = "device circuit.n_qubits 1;\nfn @main() { const.u64 0\n gate h\n const.u64 0\n gate measure\n ret }\n";

    fn module(src: &str) -> PPVMModule {
        compile_program(src).unwrap()
    }

    #[test]
    fn serial_runs_one_record_per_shot() {
        let m = module(DETERMINISTIC);
        let records = run_shots_serial(&m, 5, None).unwrap();
        assert_eq!(records.len(), 5);
        for shot in &records {
            // One measurement event, one qubit, deterministically |0>.
            assert_eq!(shot.len(), 1);
            assert_eq!(shot[0].as_slice(), [MeasurementOutcome::Zero]);
        }
    }

    #[test]
    fn dispatcher_runs_all_shots() {
        let m = module(DETERMINISTIC);
        // threads = 1 forces the serial path regardless of the rayon feature.
        let records = run_shots(&m, 10, 1, None).unwrap();
        assert_eq!(records.len(), 10);
    }

    #[test]
    fn same_seed_is_reproducible() {
        let m = module(RANDOM);
        let a = run_shots_serial(&m, 20, Some(42)).unwrap();
        let b = run_shots_serial(&m, 20, Some(42)).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn per_shot_seeds_differ() {
        // If every shot shared one seed, all 20 outcomes would be identical.
        // Distinct per-shot seeds must produce a mix of 0s and 1s.
        let m = module(RANDOM);
        let records = run_shots_serial(&m, 20, Some(42)).unwrap();
        let first = &records[0];
        assert!(
            records.iter().any(|r| r != first),
            "expected varied outcomes across shots, got {records:?}"
        );
    }

    #[cfg(feature = "rayon")]
    #[test]
    fn serial_and_parallel_match_for_same_seed() {
        let m = module(RANDOM);
        let serial = run_shots_serial(&m, 64, Some(7)).unwrap();
        let parallel = run_shots_parallel(&m, 64, 4, Some(7)).unwrap();
        assert_eq!(serial, parallel);
    }
}
