/// HashMap vs IndexMap at n=6 (1056 Pauli terms) — iteration bottleneck probe.
///
/// The flamegraph shows ~43% of solve time in hashbrown iteration.  IndexMap
/// stores entries in a dense Vec, eliminating the ~50% empty-slot overhead of
/// the open-addressed HashMap.  Both variants use fxhash to isolate the
/// data-structure change from the hash-function change.
///
/// Run with:  cargo bench --bench indexmap
use criterion::{criterion_group, criterion_main, Criterion};
use ppvm_runtime::{
    config::{fxhash::ByteF64, indexmap::ByteFxHashF64},
    prelude::*,
    strategy::CoefficientThreshold,
};
use ppvm_timeevolve::{
    JumpOp, LadderDirection, LadderOp, LindbladOp, RateMatrix,
    dopri5::step,
    rhs_into,
    solve::{SolverCache, SolverConfig},
};

const N: usize = 6;

fn rate_matrix() -> RateMatrix {
    let rates = (0..N)
        .map(|i| (0..N).map(|j| 1.0 / (1.0 + 0.1 * (i as f64 - j as f64).abs())).collect())
        .collect();
    RateMatrix::Dense(rates)
}

fn initial_zi_ct(n: usize) -> PauliSum<ByteF64<2, CoefficientThreshold>> {
    let mut p = PauliSum::builder()
        .n_qubits(n)
        .strategy(CoefficientThreshold(1e-8))
        .build();
    for i in 0..n {
        let mut zi = vec!['I'; n]; zi[i] = 'Z';
        p += (zi.into_iter().collect::<String>(), 1.0_f64);
    }
    p
}

fn initial_zi_ix(n: usize) -> PauliSum<ByteFxHashF64<2, CoefficientThreshold>> {
    let mut p = PauliSum::builder()
        .n_qubits(n)
        .strategy(CoefficientThreshold(1e-8))
        .build();
    for i in 0..n {
        let mut zi = vec!['I'; n]; zi[i] = 'Z';
        p += (zi.into_iter().collect::<String>(), 1.0_f64);
    }
    p
}

// ── HashMap + fxhash (current default) ───────────────────────────────────────

fn bench_step_hashmap(c: &mut Criterion) {
    type S = ByteF64<2, CoefficientThreshold>;
    let lop: LindbladOp<S> = LindbladOp::new(
        (0..N).map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower }))
            .collect(),
        rate_matrix(),
    );
    let y = initial_zi_ct(N);
    let cfg = SolverConfig::default();
    let mut cache = SolverCache::new(&y);
    rhs_into(None, &lop, &y, &mut cache.k[0]);

    c.bench_function("step_hashmap_n6", |b| {
        b.iter(|| step(None, &lop, &y, 0.001, &cfg, &mut cache));
    });
}

// ── IndexMap + fxhash ─────────────────────────────────────────────────────────

fn bench_step_indexmap(c: &mut Criterion) {
    type S = ByteFxHashF64<2, CoefficientThreshold>;
    let lop: LindbladOp<S> = LindbladOp::new(
        (0..N).map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower }))
            .collect(),
        rate_matrix(),
    );
    let y = initial_zi_ix(N);
    let cfg = SolverConfig::default();
    let mut cache = SolverCache::new(&y);
    rhs_into(None, &lop, &y, &mut cache.k[0]);

    c.bench_function("step_indexmap_n6", |b| {
        b.iter(|| step(None, &lop, &y, 0.001, &cfg, &mut cache));
    });
}

criterion_group!(benches, bench_step_hashmap, bench_step_indexmap);
criterion_main!(benches);
