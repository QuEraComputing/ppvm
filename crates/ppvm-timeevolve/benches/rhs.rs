use criterion::{criterion_group, criterion_main, Criterion};
use ppvm_timeevolve::{rhs, solve::solve, SolverConfig};
use std::time::Duration;

mod fixture;
use fixture::{build_initial, build_lindblad};

fn bench_rhs(c: &mut Criterion) {
    let lindblad = build_lindblad();
    let initial = build_initial();

    // Warm-up: advance state by a short solve so `p` is non-trivially sparse.
    let (_, warm_states) = solve(
        None,
        &lindblad,
        &initial,
        (0.0, 0.1),
        &[0.1],
        |_, p| p.clone(),
        SolverConfig::default(),
    );
    let p = warm_states.into_iter().next().expect("warm-up produced no save point");

    c.bench_function("bench_rhs", |b| {
        b.iter(|| rhs(None, &lindblad, &p));
    });
}

fn bench_solve(c: &mut Criterion) {
    let lindblad = build_lindblad();
    let initial = build_initial();
    let save_at: Vec<f64> = (1..=10).map(|i| i as f64 * 0.1).collect();

    // Each solve takes ~20–30 s at baseline; 10 samples is sufficient for a stable median.
    let mut group = c.benchmark_group("solve");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(300));
    group.bench_function("bench_solve", |b| {
        b.iter(|| {
            solve(
                None,
                &lindblad,
                &initial,
                (0.0, 1.0),
                &save_at,
                |_, _| (),
                SolverConfig::default(),
            )
        });
    });
    group.finish();
}

criterion_group!(benches, bench_rhs, bench_solve);
criterion_main!(benches);
