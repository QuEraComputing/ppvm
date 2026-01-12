use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use ppvm_tableau::Tableau;
use ppvm_tableau::{config::fxhash::FxComplex, TableauSum};

fn bench_gates(c: &mut Criterion) {
    let mut group = c.benchmark_group("tableau_gates");
    for &n_qubits in &[256usize, 1024] {
        group.bench_with_input(
            BenchmarkId::new("x", n_qubits),
            &n_qubits,
            |b, &n| {
                b.iter(|| {
                    let mut t = Tableau::new(n);
                    for q in 0..n {
                        t.x(q);
                    }
                    black_box(t);
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("h", n_qubits),
            &n_qubits,
            |b, &n| {
                b.iter(|| {
                    let mut t = Tableau::new(n);
                    for q in 0..n {
                        t.h(q);
                    }
                    black_box(t);
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("s", n_qubits),
            &n_qubits,
            |b, &n| {
                b.iter(|| {
                    let mut t = Tableau::new(n);
                    for q in 0..n {
                        t.s(q);
                    }
                    black_box(t);
                })
            },
        );
    }
    group.finish();
}

fn bench_gates_scalar(c: &mut Criterion) {
    let mut group = c.benchmark_group("tableau_gates_scalar");
    for &n_qubits in &[256usize, 1024] {
        group.bench_with_input(
            BenchmarkId::new("x", n_qubits),
            &n_qubits,
            |b, &n| {
                b.iter(|| {
                    let mut t = Tableau::new(n);
                    for q in 0..n {
                        t.x_scalar(q);
                    }
                    black_box(t);
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("h", n_qubits),
            &n_qubits,
            |b, &n| {
                b.iter(|| {
                    let mut t = Tableau::new(n);
                    for q in 0..n {
                        t.h_scalar(q);
                    }
                    black_box(t);
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("s", n_qubits),
            &n_qubits,
            |b, &n| {
                b.iter(|| {
                    let mut t = Tableau::new(n);
                    for q in 0..n {
                        t.s_scalar(q);
                    }
                    black_box(t);
                })
            },
        );
    }
    group.finish();
}

fn bench_t_depth(c: &mut Criterion) {
    let mut group = c.benchmark_group("tableau_t_depth");
    for &t_depth in &[4usize, 6, 8] {
        group.bench_with_input(
            BenchmarkId::new("t_depth", t_depth),
            &t_depth,
            |b, &depth| {
                b.iter(|| {
                    let mut state: TableauSum<FxComplex> = TableauSum::new(1);
                    state.h(0);
                    for _ in 0..depth {
                        state.t(0);
                    }
                    black_box(state.len());
                })
            },
        );
    }
    group.finish();
}

fn bench_t_layers(c: &mut Criterion) {
    let mut group = c.benchmark_group("tableau_t_layers");
    for &n_qubits in &[10usize, 12] {
        let layers = 1usize;
        group.bench_with_input(
            BenchmarkId::new(format!("n{}_d{}", n_qubits, layers), n_qubits),
            &(n_qubits, layers),
            |b, &(n, depth)| {
                b.iter(|| {
                        let mut state: TableauSum<FxComplex> = TableauSum::new(n);
                    for q in 0..n {
                        state.h(q);
                    }
                    for _ in 0..depth {
                        for q in 0..n {
                            state.t(q);
                        }
                    }
                    black_box(state.len());
                })
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_gates, bench_gates_scalar, bench_t_depth, bench_t_layers);
criterion_main!(benches);
