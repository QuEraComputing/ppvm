/// Benchmarks for a single `step()` call across all three truncation strategies.
///
/// These measure the per-step cost of truncation at each of the 6 truncation
/// points introduced in Task 26 (stages 2–6 + y5).  The fixture matches
/// `bench_rhs`: n=5, lowering operators, dense 5×5 γ.
///
/// Results are used in Task 28 to decide whether `Budget` is recommended as
/// a default strategy (see the comment above `Budget` in `src/strategy.rs`).
use criterion::{criterion_group, criterion_main, Criterion};
use ppvm_runtime::{config::fxhash::ByteF64, prelude::*, strategy::{CoefficientThreshold, MaxPauliWeight}};
use ppvm_timeevolve::{
    Budget, JumpOp, LadderDirection, LadderOp, LindbladOp, RateMatrix,
    dopri5::step,
    rhs_into,
    solve::{SolverCache, SolverConfig},
};

const N: usize = 5;

fn rate_matrix() -> RateMatrix {
    let rates = (0..N)
        .map(|i| (0..N).map(|j| 1.0 / (1.0 + (i as f64 - j as f64).abs())).collect())
        .collect();
    RateMatrix::Dense(rates)
}

// ── CoefficientThreshold ──────────────────────────────────────────────────────

fn bench_step_ct(c: &mut Criterion) {
    type S = ByteF64<1, CoefficientThreshold>;
    let lop: LindbladOp<S> = LindbladOp::new(
        (0..N).map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower })).collect(),
        rate_matrix(),
    );
    let mut y: PauliSum<S> = PauliSum::builder()
        .n_qubits(N).strategy(CoefficientThreshold(1e-6)).build();
    for i in 0..N {
        let mut zi = vec!['I'; N]; zi[i] = 'Z';
        y += (zi.into_iter().collect::<String>(), 1.0_f64);
    }
    let cfg = SolverConfig::default();
    let mut cache = SolverCache::new(&y);
    rhs_into(None, &lop, &y, &mut cache.k[0]);

    c.bench_function("bench_step_ct", |b| {
        b.iter(|| step(None, &lop, &y, 0.01, &cfg, &mut cache));
    });
}

// ── MaxPauliWeight ────────────────────────────────────────────────────────────

fn bench_step_mpw(c: &mut Criterion) {
    type S = ByteF64<1, MaxPauliWeight>;
    let lop: LindbladOp<S> = LindbladOp::new(
        (0..N).map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower })).collect(),
        rate_matrix(),
    );
    let mut y: PauliSum<S> = PauliSum::builder()
        .n_qubits(N).strategy(MaxPauliWeight(2)).build();
    for i in 0..N {
        let mut zi = vec!['I'; N]; zi[i] = 'Z';
        y += (zi.into_iter().collect::<String>(), 1.0_f64);
    }
    let cfg = SolverConfig::default();
    let mut cache = SolverCache::new(&y);
    rhs_into(None, &lop, &y, &mut cache.k[0]);

    c.bench_function("bench_step_mpw", |b| {
        b.iter(|| step(None, &lop, &y, 0.01, &cfg, &mut cache));
    });
}

// ── Budget ────────────────────────────────────────────────────────────────────

fn bench_step_budget(c: &mut Criterion) {
    type S = ByteF64<1, Budget>;
    let lop: LindbladOp<S> = LindbladOp::new(
        (0..N).map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower })).collect(),
        rate_matrix(),
    );
    let mut y: PauliSum<S> = PauliSum::builder()
        .n_qubits(N).strategy(Budget { target: 300 }).build();
    for i in 0..N {
        let mut zi = vec!['I'; N]; zi[i] = 'Z';
        y += (zi.into_iter().collect::<String>(), 1.0_f64);
    }
    let cfg = SolverConfig::default();
    let mut cache = SolverCache::new(&y);
    rhs_into(None, &lop, &y, &mut cache.k[0]);

    c.bench_function("bench_step_budget", |b| {
        b.iter(|| step(None, &lop, &y, 0.01, &cfg, &mut cache));
    });
}

criterion_group!(benches, bench_step_ct, bench_step_mpw, bench_step_budget);
criterion_main!(benches);
