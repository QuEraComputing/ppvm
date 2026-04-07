/// Superradiance flamegraph target.
///
/// Identical physics to `superradiance`, but runs the baseline (untruncated)
/// solve in a tight loop so the profiler accumulates enough samples for a
/// meaningful flamegraph.
///
/// N=6 gives 36 Lindblad terms and a state of ~1056 Pauli terms at T=0.5 —
/// large enough to give meaningful profile samples without the exponential cost
/// of n=8. Each solve takes ~200ms; 30 loops gives ~6s of profiling time.
///
/// Run with:
///   cargo flamegraph --example superradiance_flame --release
use ppvm_runtime::{
    config::indexmap::ByteFxHashF64,
    prelude::*,
    strategy::{CoefficientThreshold, CombinedStrategy, MaxPauliWeight},
};
use ppvm_timeevolve::{
    JumpOp, LadderDirection, LadderOp, LindbladOp, RateMatrix, SolverConfig,
    solve::solve,
};

const N:      usize = 6;
const NBYTES: usize = 2;
const GAMMA0: f64   = 1.0;
const D:      f64   = 0.1;
const TMAX:   f64   = 0.5;
const TSTEPS: usize = 20;
const LOOPS:  usize = 30;

// CombinedStrategy activates both the fused weight filter (MaxPauliWeight) and
// coefficient pruning (CoefficientThreshold). MaxPauliWeight caps the Pauli weight
// of generated terms, cutting the cross-site kernel work at the weight ceiling.
type SB = ByteFxHashF64<NBYTES, CombinedStrategy<MaxPauliWeight, CoefficientThreshold>>;

fn rate_matrix() -> RateMatrix {
    RateMatrix::Dense(
        (0..N)
            .map(|i| {
                (0..N)
                    .map(|j| GAMMA0 / (1.0 + D * (i as f64 - j as f64).abs()))
                    .collect()
            })
            .collect(),
    )
}

fn build_ops() -> LindbladOp<SB> {
    let ops: Vec<JumpOp<SB>> = (0..N)
        .map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Raise }))
        .collect();
    LindbladOp::new(ops, rate_matrix())
}

fn initial_state() -> PauliSum<SB> {
    let mut p = PauliSum::builder()
        .n_qubits(N)
        .strategy(CombinedStrategy(MaxPauliWeight(4), CoefficientThreshold(1e-8)))
        .build();
    let t = vec!['I'; N];
    for i in 0..N {
        let mut zi = t.clone();
        zi[i] = 'Z';
        p += (zi.iter().collect::<String>(), 1.0_f64);
    }
    p
}

fn main() {
    let lindblad  = build_ops();
    let save_at: Vec<f64> = (1..=TSTEPS).map(|i| i as f64 * TMAX / TSTEPS as f64).collect();
    let pattern: PauliPattern = "Z?*".into();

    // Warm-up: one solve to populate caches before profiling.
    let _ = solve(
        None, &lindblad, &initial_state(),
        (0.0, TMAX), &save_at,
        |_, p: &PauliSum<SB>| p.trace(&pattern),
        SolverConfig::default(),
    );

    // Hot loop: this is what the flamegraph captures.
    for _ in 0..LOOPS {
        let _ = solve(
            None, &lindblad, &initial_state(),
            (0.0, TMAX), &save_at,
            |_, p: &PauliSum<SB>| p.trace(&pattern),
            SolverConfig::default(),
        );
    }
}
