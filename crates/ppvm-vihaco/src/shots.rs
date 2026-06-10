// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Running a compiled program for many shots, optionally across threads.
//!
//! Each shot runs on a fresh [`PPVM`] so shots are fully independent; the
//! module is compiled once and shared. With the `rayon` feature, [`run_shots`]
//! parallelizes across shots when the global pool (sized once by
//! [`set_global_threads`]) has more than one thread and there are enough shots
//! to amortize the overhead.

use crate::PPVMModule;
use crate::composite::PPVM;
use crate::measurements::MeasurementResult;

/// One shot's full output: the measurement record and the trace-instruction
/// record. Either may be empty depending on what the program emits.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ShotRecord {
    pub measurements: Vec<MeasurementResult>,
    pub traces: Vec<f64>,
}

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

/// Run a single shot on a fresh machine and return both records.
fn run_one_shot(module: &PPVMModule, seed: Option<u64>) -> eyre::Result<ShotRecord> {
    let mut machine = PPVM::default();
    machine.load(module)?;
    machine.run_with_seed(seed)?;
    Ok(ShotRecord {
        measurements: machine.measurement_record(),
        traces: machine.trace_record(),
    })
}

/// Run `shots` shots serially. One entry per shot, in order.
pub fn run_shots_serial(
    module: &PPVMModule,
    shots: usize,
    seed: Option<u64>,
) -> eyre::Result<Vec<ShotRecord>> {
    (0..shots)
        .map(|i| run_one_shot(module, shot_seed(seed, i)))
        .collect()
}

/// Run `shots` shots across the global rayon pool. One entry per shot, in order
/// (preserved by the indexed parallel iterator). The pool size is whatever
/// [`set_global_threads`] configured at startup; each shot runs on a worker
/// thread, so the intra-shot parallelism guard keeps a single shot serial and
/// the pool is never oversubscribed.
#[cfg(feature = "rayon")]
pub fn run_shots_parallel(
    module: &PPVMModule,
    shots: usize,
    seed: Option<u64>,
) -> eyre::Result<Vec<ShotRecord>> {
    use rayon::prelude::*;

    (0..shots)
        .into_par_iter()
        .map(|i| run_one_shot(module, shot_seed(seed, i)))
        .collect()
}

/// Decide whether to spread shots across the rayon pool. Worth it only with a
/// multi-thread pool and enough shots to amortize the overhead; a single-thread
/// pool always takes the serial path, keeping results deterministic.
#[cfg(feature = "rayon")]
fn should_parallelize(num_threads: usize, shots: usize) -> bool {
    num_threads > 1 && shots >= PARALLEL_SHOT_THRESHOLD
}

/// Run `shots` shots, choosing serial or parallel execution. Goes parallel only
/// when built with `rayon`, the global pool has more than one thread, and there
/// are enough shots to be worth it; otherwise runs serially. The pool size is
/// set once at startup by [`set_global_threads`].
pub fn run_shots(
    module: &PPVMModule,
    shots: usize,
    seed: Option<u64>,
) -> eyre::Result<Vec<ShotRecord>> {
    #[cfg(feature = "rayon")]
    if should_parallelize(rayon::current_num_threads(), shots) {
        return run_shots_parallel(module, shots, seed);
    }

    run_shots_serial(module, shots, seed)
}

/// Configure the process-wide rayon thread pool. Call once, before any parallel
/// work runs. A count of `1` forces fully serial, deterministic execution — both
/// across shots and within a single machine's coefficient propagation.
#[cfg(feature = "rayon")]
pub fn set_global_threads(threads: usize) -> eyre::Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?;
    Ok(())
}

/// Without the `rayon` feature there is no pool to size; anything but a single
/// thread is meaningless, so reject it rather than silently run serially.
#[cfg(not(feature = "rayon"))]
pub fn set_global_threads(threads: usize) -> eyre::Result<()> {
    if threads > 1 {
        eyre::bail!("this build has no parallelism support; --threads must be 1");
    }
    Ok(())
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
            assert_eq!(shot.measurements.len(), 1);
            assert_eq!(shot.measurements[0].as_slice(), [MeasurementOutcome::Zero]);
            assert!(shot.traces.is_empty(), "Tableau backend emits no traces");
        }
    }

    #[test]
    fn dispatcher_runs_all_shots() {
        let m = module(DETERMINISTIC);
        // 10 shots is below the parallel threshold, so this takes the serial path.
        let records = run_shots(&m, 10, None).unwrap();
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
        let parallel = run_shots_parallel(&m, 64, Some(7)).unwrap();
        assert_eq!(serial, parallel);
    }

    #[cfg(feature = "rayon")]
    #[test]
    fn parallelizes_only_with_multiple_threads_above_threshold() {
        // A single-thread pool is always serial, no matter how many shots.
        assert!(!should_parallelize(1, 100_000));
        // Multiple threads, but too few shots to amortize the overhead: serial.
        assert!(!should_parallelize(8, PARALLEL_SHOT_THRESHOLD - 1));
        // Multiple threads at the threshold: parallel.
        assert!(should_parallelize(8, PARALLEL_SHOT_THRESHOLD));
        // Multiple threads, plenty of shots: parallel.
        assert!(should_parallelize(8, 100_000));
    }
}
