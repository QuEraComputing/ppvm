/// Criterion benchmarks for the fused cross-site Lindblad kernel (Phase 5).
///
/// `bench_rhs_fused_n6`        — n=6, no weight cap; establishes untruncated baseline.
/// `bench_rhs_fused_n6_wmax3`  — n=6, MaxPauliWeight(3)+CoefficientThreshold; filter benefit.
/// `bench_rhs_fused_n10_wmax4` — n=10, MaxPauliWeight(4); larger system.
/// `bench_solve_superradiance_n6` — full solve n=6; end-to-end integration cost.
use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;
use ppvm_runtime::{
    config::indexmap::ByteFxHashF64,
    prelude::*,
    strategy::{CoefficientThreshold, CombinedStrategy, MaxPauliWeight},
};
use ppvm_timeevolve::{
    JumpOp, LadderDirection, LadderOp, LindbladOp, RateMatrix, SolverConfig,
    rhs, solve::solve,
};

type S6    = ByteFxHashF64<1, CoefficientThreshold>;
type S6W3  = ByteFxHashF64<1, CombinedStrategy<MaxPauliWeight, CoefficientThreshold>>;
type S10W4 = ByteFxHashF64<2, CombinedStrategy<MaxPauliWeight, CoefficientThreshold>>;

fn dense_rate_matrix(n: usize) -> RateMatrix {
    let rates = (0..n)
        .map(|i| (0..n).map(|j| 1.0 / (1.0 + (i as f64 - j as f64).abs())).collect())
        .collect();
    RateMatrix::Dense(rates)
}

// ── bench_rhs_fused_n6 ────────────────────────────────────────────────────────

fn bench_rhs_fused_n6(c: &mut Criterion) {
    const N: usize = 6;
    let lop: LindbladOp<S6> = LindbladOp::new(
        (0..N).map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower })).collect(),
        dense_rate_matrix(N),
    );
    let mut p: PauliSum<S6> = PauliSum::builder().n_qubits(N).strategy(CoefficientThreshold(1e-6)).build();
    for i in 0..N { let mut z = vec!['I'; N]; z[i] = 'Z'; p += (z.iter().collect::<String>(), 1.0_f64); }
    let (_, states) = solve(None, &lop, &p, (0.0, 0.1), &[0.1], |_, s| s.clone(), SolverConfig::default());
    let p_warm = states.into_iter().next().unwrap();
    c.bench_function("bench_rhs_fused_n6", |b| {
        b.iter(|| rhs(None, &lop, &p_warm));
    });
}

// ── bench_rhs_fused_n6_wmax3 ─────────────────────────────────────────────────

fn bench_rhs_fused_n6_wmax3(c: &mut Criterion) {
    const N: usize = 6;
    let lop: LindbladOp<S6W3> = LindbladOp::new(
        (0..N).map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower })).collect(),
        dense_rate_matrix(N),
    );
    let strat = CombinedStrategy(MaxPauliWeight(3), CoefficientThreshold(1e-6));
    let mut p: PauliSum<S6W3> = PauliSum::builder().n_qubits(N).strategy(strat).build();
    for i in 0..N { let mut z = vec!['I'; N]; z[i] = 'Z'; p += (z.iter().collect::<String>(), 1.0_f64); }
    let (_, states) = solve(None, &lop, &p, (0.0, 0.1), &[0.1], |_, s| s.clone(), SolverConfig::default());
    let p_warm = states.into_iter().next().unwrap();
    c.bench_function("bench_rhs_fused_n6_wmax3", |b| {
        b.iter(|| rhs(None, &lop, &p_warm));
    });
}

// ── bench_rhs_fused_n10_wmax4 ────────────────────────────────────────────────

fn bench_rhs_fused_n10_wmax4(c: &mut Criterion) {
    const N: usize = 10;
    let lop: LindbladOp<S10W4> = LindbladOp::new(
        (0..N).map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower })).collect(),
        dense_rate_matrix(N),
    );
    let strat = CombinedStrategy(MaxPauliWeight(4), CoefficientThreshold(1e-6));
    let mut p: PauliSum<S10W4> = PauliSum::builder().n_qubits(N).strategy(strat).build();
    for i in 0..N { let mut z = vec!['I'; N]; z[i] = 'Z'; p += (z.iter().collect::<String>(), 1.0_f64); }
    let (_, states) = solve(None, &lop, &p, (0.0, 0.1), &[0.1], |_, s| s.clone(), SolverConfig::default());
    let p_warm = states.into_iter().next().unwrap();
    c.bench_function("bench_rhs_fused_n10_wmax4", |b| {
        b.iter(|| rhs(None, &lop, &p_warm));
    });
}

// ── bench_solve_superradiance_n6 ─────────────────────────────────────────────

fn bench_solve_superradiance_n6(c: &mut Criterion) {
    const N: usize = 6;
    const GAMMA0: f64 = 1.0;
    const D: f64 = 0.1;
    const TMAX: f64 = 0.5;

    let rates: Vec<Vec<f64>> = (0..N)
        .map(|i| (0..N).map(|j| GAMMA0 / (1.0 + D * (i as f64 - j as f64).abs())).collect())
        .collect();
    let lop: LindbladOp<S6> = LindbladOp::new(
        (0..N).map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Raise })).collect(),
        RateMatrix::Dense(rates),
    );
    let mut p0: PauliSum<S6> = PauliSum::builder().n_qubits(N).strategy(CoefficientThreshold(1e-8)).build();
    for i in 0..N { let mut z = vec!['I'; N]; z[i] = 'Z'; p0 += (z.iter().collect::<String>(), 1.0_f64); }
    let save_at: Vec<f64> = (1..=20).map(|i| i as f64 * TMAX / 20.0).collect();
    let pattern: PauliPattern = "Z?*".into();

    let mut group = c.benchmark_group("superradiance");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(60));
    group.bench_function("bench_solve_superradiance_n6", |b| {
        b.iter(|| {
            solve(
                None, &lop, &p0,
                (0.0, TMAX), &save_at,
                |_, p: &PauliSum<S6>| p.trace(&pattern),
                SolverConfig::default(),
            )
        });
    });
    group.finish();
}

// ── bench_solve_superradiance_n6_wmax4 (flame example configuration) ─────────

fn bench_solve_superradiance_n6_wmax4(c: &mut Criterion) {
    const N: usize = 6;
    const GAMMA0: f64 = 1.0;
    const D: f64 = 0.1;
    const TMAX: f64 = 0.5;

    let rates: Vec<Vec<f64>> = (0..N)
        .map(|i| (0..N).map(|j| GAMMA0 / (1.0 + D * (i as f64 - j as f64).abs())).collect())
        .collect();
    let lop: LindbladOp<S6W3> = LindbladOp::new(
        (0..N).map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Raise })).collect(),
        RateMatrix::Dense(rates),
    );
    let strat = CombinedStrategy(MaxPauliWeight(4), CoefficientThreshold(1e-8));
    let mut p0: PauliSum<S6W3> = PauliSum::builder().n_qubits(N).strategy(strat).build();
    for i in 0..N { let mut z = vec!['I'; N]; z[i] = 'Z'; p0 += (z.iter().collect::<String>(), 1.0_f64); }
    let save_at: Vec<f64> = (1..=20).map(|i| i as f64 * TMAX / 20.0).collect();
    let pattern: PauliPattern = "Z?*".into();

    let mut group = c.benchmark_group("superradiance");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(60));
    group.bench_function("bench_solve_superradiance_n6_wmax4", |b| {
        b.iter(|| {
            solve(
                None, &lop, &p0,
                (0.0, TMAX), &save_at,
                |_, p: &PauliSum<S6W3>| p.trace(&pattern),
                SolverConfig::default(),
            )
        });
    });
    group.finish();
}

// ── bench_solve_superradiance_n6_generic (pre-Phase-5 baseline) ──────────────

fn bench_solve_superradiance_n6_generic(c: &mut Criterion) {
    const N: usize = 6;
    const GAMMA0: f64 = 1.0;
    const D: f64 = 0.1;
    const TMAX: f64 = 0.5;

    let rates: Vec<Vec<f64>> = (0..N)
        .map(|i| (0..N).map(|j| GAMMA0 / (1.0 + D * (i as f64 - j as f64).abs())).collect())
        .collect();
    // Generic-expanded Ladder ops: what LindbladOp used before Phase 5.
    let lop: LindbladOp<S6> = LindbladOp::new(
        (0..N).map(|i| JumpOp::Generic(
            LadderOp { qubit: i, direction: LadderDirection::Raise }.expand::<S6>(N)
        )).collect(),
        RateMatrix::Dense(rates),
    );
    let mut p0: PauliSum<S6> = PauliSum::builder().n_qubits(N).strategy(CoefficientThreshold(1e-8)).build();
    for i in 0..N { let mut z = vec!['I'; N]; z[i] = 'Z'; p0 += (z.iter().collect::<String>(), 1.0_f64); }
    let save_at: Vec<f64> = (1..=20).map(|i| i as f64 * TMAX / 20.0).collect();
    let pattern: PauliPattern = "Z?*".into();

    let mut group = c.benchmark_group("superradiance");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(60));
    group.bench_function("bench_solve_superradiance_n6_generic", |b| {
        b.iter(|| {
            solve(
                None, &lop, &p0,
                (0.0, TMAX), &save_at,
                |_, p: &PauliSum<S6>| p.trace(&pattern),
                SolverConfig::default(),
            )
        });
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_rhs_fused_n6,
    bench_rhs_fused_n6_wmax3,
    bench_rhs_fused_n10_wmax4,
    bench_solve_superradiance_n6,
    bench_solve_superradiance_n6_wmax4,
    bench_solve_superradiance_n6_generic,
);
criterion_main!(benches);
