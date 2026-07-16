// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Running a compiled program for many shots, optionally across threads.
//!
//! Each shot runs on a fresh [`PPVM`] so shots are fully independent; the
//! module is compiled once and cloned into each shot's machine. With the `rayon` feature,
//! [`run_shots`] parallelizes across shots when the global pool (sized once by
//! [`set_global_threads`]) has more than one thread and there are enough shots
//! to amortize the overhead.

use crate::PPVMModule;
use crate::composite::{BackendKind, PPVM, PPVMInstruction, PPVMSnapshot, StepOutcome};
use crate::measurements::MeasurementResult;
use ppvm_trajectory_cache::{
    CacheConfig, CacheStats, CachedRun, TrajectoryEvent, TrajectoryProgram, random_base_seed,
    run_cached_shots,
};
use vihaco_circuit_isa::CircuitInstruction;

/// One shot's full output: the measurement record and the trace-instruction
/// record. Either may be empty depending on what the program emits.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ShotRecord {
    pub measurements: Vec<MeasurementResult>,
    pub traces: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ShotOptions {
    pub seed: Option<u64>,
    pub cache: CacheConfig,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ShotBatch {
    pub records: Vec<ShotRecord>,
    pub cache_stats: Option<CacheStats>,
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

pub fn run_shots_with_options(
    module: &PPVMModule,
    shots: usize,
    options: ShotOptions,
) -> eyre::Result<ShotBatch> {
    if !options.cache.enabled {
        return Ok(ShotBatch {
            records: run_shots(module, shots, options.seed)?,
            cache_stats: None,
        });
    }

    if module.extra.backend != BackendKind::Tableau {
        eyre::bail!("trajectory cache currently supports only the tableau backend");
    }
    reject_deferred_side_effects(module)?;

    let base_seed = options.seed.unwrap_or_else(random_base_seed);
    let mut program = VihacoTrajectoryProgram::new(module, base_seed);
    let CachedRun {
        output,
        cache_stats,
    } = run_cached_shots(&mut program, shots, options.cache, base_seed)?;
    Ok(ShotBatch {
        records: output,
        cache_stats: Some(cache_stats),
    })
}

fn reject_deferred_side_effects(module: &PPVMModule) -> eyre::Result<()> {
    for inst in &module.code {
        if matches!(inst, PPVMInstruction::Cpu(vihaco_cpu::Instruction::Print)) {
            eyre::bail!("trajectory cache does not yet support print side effects");
        }
    }
    Ok(())
}

struct VihacoTrajectoryProgram<'a> {
    module: &'a PPVMModule,
    machine: PPVM,
    base_seed: u64,
}

impl<'a> VihacoTrajectoryProgram<'a> {
    fn new(module: &'a PPVMModule, base_seed: u64) -> Self {
        Self {
            module,
            machine: PPVM::default(),
            base_seed,
        }
    }

    fn shot_seed(&self, shot: usize) -> u64 {
        self.base_seed.wrapping_add(shot as u64)
    }

    fn shot_record(&self) -> ShotRecord {
        ShotRecord {
            measurements: self.machine.measurement_record(),
            traces: self.machine.trace_record(),
        }
    }
}

impl TrajectoryProgram for VihacoTrajectoryProgram<'_> {
    type Snapshot = PPVMSnapshot;
    type Choice = Vec<u8>;
    type Output = ShotRecord;
    type Error = eyre::Report;

    fn reset_for_shot(&mut self, shot: usize) -> eyre::Result<()> {
        self.machine = PPVM::default();
        self.machine.load(self.module)?;
        self.machine.init_with_seed(self.shot_seed(shot))?;
        Ok(())
    }

    fn snapshot(&self) -> Self::Snapshot {
        self.machine.cache_snapshot()
    }

    fn restore(&mut self, snapshot: &Self::Snapshot) -> eyre::Result<()> {
        self.machine.restore_cache_snapshot(snapshot);
        Ok(())
    }

    fn reseed(&mut self, seed: u64) -> eyre::Result<()> {
        self.machine.reseed_cache_rng(seed)
    }

    fn run_until_boundary(&mut self) -> eyre::Result<TrajectoryEvent<Self::Output>> {
        loop {
            let Some(inst) = self.machine.current_instruction() else {
                return Ok(TrajectoryEvent::Terminal(self.shot_record()));
            };
            if is_supported_boundary(&inst) {
                return Ok(TrajectoryEvent::Boundary);
            }
            reject_unsupported_stochastic(&inst)?;
            match self.machine.step_once()? {
                StepOutcome::Continue | StepOutcome::Breakpoint => {}
                StepOutcome::Return | StepOutcome::Halt => {
                    return Ok(TrajectoryEvent::Terminal(self.shot_record()));
                }
            }
        }
    }

    fn execute_boundary(&mut self) -> eyre::Result<Self::Choice> {
        self.machine.step_cache_boundary()
    }
}

fn is_supported_boundary(inst: &PPVMInstruction) -> bool {
    matches!(
        inst,
        PPVMInstruction::Circuit(
            CircuitInstruction::Measure
                | CircuitInstruction::Depolarize
                | CircuitInstruction::Depolarize2
                | CircuitInstruction::PauliError
                | CircuitInstruction::TwoQubitPauliError
                | CircuitInstruction::Loss
                | CircuitInstruction::CorrelatedLoss
        )
    )
}

fn reject_unsupported_stochastic(inst: &PPVMInstruction) -> eyre::Result<()> {
    match inst {
        PPVMInstruction::Circuit(op @ CircuitInstruction::Reset) => {
            eyre::bail!("trajectory cache does not yet support hidden stochastic instruction {op}");
        }
        _ => Ok(()),
    }
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
        "device circuit.n_qubits 1;\nfn @main() { const.u64 0\n circuit.measure\n ret }\n";

    /// Prepares |+> with H, then measures q0: each shot is a random 0/1.
    const RANDOM: &str = "device circuit.n_qubits 1;\nfn @main() { const.u64 0\n circuit.h\n const.u64 0\n circuit.measure\n ret }\n";

    const NOISY: &str = "device circuit.n_qubits 1;\n\
        fn @main() {\n\
            const.u64 0\n\
            const.f64 0.25\n\
            const.f64 0.0\n\
            const.f64 0.0\n\
            circuit.paulierror\n\
            const.u64 0\n\
            circuit.measure\n\
            ret\n\
        }\n";

    const LOSSY: &str = "device circuit.n_qubits 1;\n\
        fn @main() {\n\
            const.u64 0\n\
            const.f64 0.5\n\
            circuit.loss\n\
            const.u64 0\n\
            circuit.measure\n\
            ret\n\
        }\n";

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

    #[test]
    fn cached_runner_reuses_measurement_paths() {
        let m = module(RANDOM);
        let batch = run_shots_with_options(
            &m,
            64,
            ShotOptions {
                seed: Some(42),
                cache: ppvm_trajectory_cache::CacheConfig::bounded(16),
            },
        )
        .unwrap();

        assert_eq!(batch.records.len(), 64);
        let stats = batch.cache_stats.expect("cache enabled");
        assert!(stats.hits > 0, "expected cache hits, got {stats:?}");
        assert!(stats.nodes <= 16, "bounded cache exceeded limit: {stats:?}");
    }

    #[test]
    fn cached_runner_supports_noise_boundaries() {
        let m = module(NOISY);
        let batch = run_shots_with_options(
            &m,
            256,
            ShotOptions {
                seed: Some(42),
                cache: ppvm_trajectory_cache::CacheConfig::bounded(16),
            },
        )
        .unwrap();

        assert_eq!(batch.records.len(), 256);
        let stats = batch.cache_stats.expect("cache enabled");
        assert!(stats.hits > 0, "expected cache hits, got {stats:?}");
        assert!(
            batch
                .records
                .iter()
                .any(|r| r.measurements[0].as_slice() == [MeasurementOutcome::Zero]),
            "expected some no-error shots"
        );
        assert!(
            batch
                .records
                .iter()
                .any(|r| r.measurements[0].as_slice() == [MeasurementOutcome::One]),
            "expected some Pauli-error shots"
        );
    }

    #[test]
    fn cached_runner_supports_loss_boundaries() {
        let m = module(LOSSY);
        let batch = run_shots_with_options(
            &m,
            256,
            ShotOptions {
                seed: Some(42),
                cache: ppvm_trajectory_cache::CacheConfig::bounded(16),
            },
        )
        .unwrap();

        assert_eq!(batch.records.len(), 256);
        let stats = batch.cache_stats.expect("cache enabled");
        assert!(stats.hits > 0, "expected cache hits, got {stats:?}");
        assert!(
            batch
                .records
                .iter()
                .any(|r| r.measurements[0].as_slice() == [MeasurementOutcome::Zero]),
            "expected some surviving shots"
        );
        assert!(
            batch
                .records
                .iter()
                .any(|r| r.measurements[0].as_slice() == [MeasurementOutcome::Lost]),
            "expected some lost shots"
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
