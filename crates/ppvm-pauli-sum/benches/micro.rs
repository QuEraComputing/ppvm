// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_pauli_sum::config::fxhash::ByteF64;
use ppvm_pauli_sum::prelude::*;
use ppvm_pauli_sum::strategy::CoefficientThreshold;

type Cfg = ByteF64<2, CoefficientThreshold>;

fn configure() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50)
}

/// Build a PauliSum with `n_terms` random-ish terms on `n_qubits` qubits.
fn make_state(n_qubits: usize, n_terms: usize) -> PauliSum<Cfg> {
    let strat = CoefficientThreshold(1e-10);
    let mut state: PauliSum<Cfg> = PauliSum::builder()
        .n_qubits(n_qubits)
        .capacity(n_terms * 2)
        .strategy(strat)
        .build();

    // Seed terms with distinct Pauli strings by cycling through X/Y/Z positions
    let paulis = [Pauli::X, Pauli::Y, Pauli::Z];
    for k in 0..n_terms {
        let mut word = <Cfg as Config>::PauliWordType::new(n_qubits);
        // Place non-identity Paulis at 2-3 positions based on k
        word.set(k % n_qubits, paulis[k % 3]);
        word.set((k * 3 + 1) % n_qubits, paulis[(k + 1) % 3]);
        if k % 2 == 0 {
            word.set((k * 7 + 2) % n_qubits, paulis[(k + 2) % 3]);
        }
        state += (word, 1.0 / (k as f64 + 1.0));
    }
    state
}

const N_QUBITS: usize = 8;

// ── Single-qubit Clifford gates ─────────────────────────────────────────────

fn bench_clifford_single(c: &mut Criterion) {
    let state = make_state(N_QUBITS, 100);
    let mut group = c.benchmark_group("clifford/single");

    group.bench_function("x", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.x(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("y", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.y(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("z", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.z(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("h", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.h(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("s", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.s(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("s_dag", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.s_dag(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("sqrt_x", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.sqrt_x(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("sqrt_x_dag", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.sqrt_x_dag(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("sqrt_y", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.sqrt_y(0),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("sqrt_y_dag", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.sqrt_y_dag(0),
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}

// ── Two-qubit Clifford gates ────────────────────────────────────────────────

fn bench_clifford_two(c: &mut Criterion) {
    let state = make_state(N_QUBITS, 100);
    let mut group = c.benchmark_group("clifford/two-qubit");

    group.bench_function("cnot", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.cnot([0, 1]),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("cz", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.cz([0, 1]),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("cy", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.cy([0, 1]),
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}

// ── Rotation gates ──────────────────────────────────────────────────────────

fn bench_rotations(c: &mut Criterion) {
    let state = make_state(N_QUBITS, 100);
    let mut group = c.benchmark_group("rotations");

    group.bench_function("rx", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.rx(0, 0.5),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("ry", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.ry(0, 0.5),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("rz", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.rz(0, 0.5),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("rxx", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.rxx([0, 1], 0.5),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("ryy", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.ryy([0, 1], 0.5),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("rzz", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.rzz([0, 1], 0.5),
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}

// ── Noise channels ──────────────────────────────────────────────────────────

fn bench_noise(c: &mut Criterion) {
    let state = make_state(N_QUBITS, 100);
    let mut group = c.benchmark_group("noise");

    group.bench_function("pauli_error", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.pauli_error(0, [0.1, 0.1, 0.1]),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("two_qubit_pauli_error", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.two_qubit_pauli_error([0, 1], [1.0 / 15.0; 15]),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("depolarize", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.depolarize1(0, 0.1),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("depolarize2", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.depolarize2([0, 1], 0.1),
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("amplitude_damping", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.amplitude_damping(0, 0.1),
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}

// ── Truncation ──────────────────────────────────────────────────────────────

fn bench_truncation(c: &mut Criterion) {
    // Build a state with terms that have varying coefficients so truncation
    // actually removes some entries
    let n_qubits = N_QUBITS;
    let strat = CoefficientThreshold(0.05);
    let mut state: PauliSum<Cfg> = PauliSum::builder()
        .n_qubits(n_qubits)
        .capacity(200)
        .strategy(strat)
        .build();
    let paulis = [Pauli::X, Pauli::Y, Pauli::Z];
    for k in 0..100 {
        let mut word = <Cfg as Config>::PauliWordType::new(n_qubits);
        word.set(k % n_qubits, paulis[k % 3]);
        word.set((k * 3 + 1) % n_qubits, paulis[(k + 1) % 3]);
        // Every 4th term gets a small coefficient that truncation will remove
        let coeff = if k % 4 == 0 {
            0.001
        } else {
            1.0 / (k as f64 + 1.0)
        };
        state += (word, coeff);
    }

    let mut group = c.benchmark_group("truncation");

    group.bench_function("coefficient_threshold", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |s| s.truncate(),
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}

// ── Trace / pattern matching ────────────────────────────────────────────────

fn bench_trace(c: &mut Criterion) {
    let state = make_state(N_QUBITS, 100);
    let mut group = c.benchmark_group("trace");

    let pat_z_star: PauliPattern = "Z?*".into();
    let pat_exact: PauliPattern = PauliPattern::parse("X0Z1").unwrap();

    group.bench_function("overlap_with_zero", |b| {
        b.iter(|| state.trace(&pat_z_star));
    });
    group.bench_function("pattern/exact", |b| {
        b.iter(|| state.trace(&pat_exact));
    });

    group.finish();
}

// ── PauliWord operations ────────────────────────────────────────────────────

fn bench_pauli_word(c: &mut Criterion) {
    let mut group = c.benchmark_group("pauli-word");
    let n = N_QUBITS;

    group.bench_function("create", |b| {
        b.iter(|| {
            let mut w = <Cfg as Config>::PauliWordType::new(n);
            w.set(0, Pauli::X);
            w.set(1, Pauli::Z);
            w.set(2, Pauli::Y);
            w
        });
    });

    let w: <Cfg as Config>::PauliWordType = "XYZIIIII".into();
    group.bench_function("get", |b| {
        b.iter(|| {
            let mut sum = 0u8;
            for i in 0..n {
                sum = sum.wrapping_add(w.get(i) as u8);
            }
            sum
        });
    });
    group.bench_function("weight", |b| {
        b.iter(|| w.weight());
    });

    let w1: <Cfg as Config>::PauliWordType = "XYZIIIII".into();
    let w2: <Cfg as Config>::PauliWordType = "XYZIIIII".into();
    let w3: <Cfg as Config>::PauliWordType = "IIIIIXYZ".into();
    group.bench_function("eq/same", |b| {
        b.iter(|| w1 == w2);
    });
    group.bench_function("eq/different", |b| {
        b.iter(|| w1 == w3);
    });

    group.finish();
}

// ── Scaling: gate cost vs state size ────────────────────────────────────────

fn bench_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("scaling");

    for &n_terms in &[10, 100, 1000] {
        let state = make_state(N_QUBITS, n_terms);
        group.bench_function(format!("h/{n_terms}"), |b| {
            b.iter_batched_ref(
                || state.clone(),
                |s| s.h(0),
                criterion::BatchSize::SmallInput,
            );
        });
        group.bench_function(format!("rx/{n_terms}"), |b| {
            b.iter_batched_ref(
                || state.clone(),
                |s| s.rx(0, 0.5),
                criterion::BatchSize::SmallInput,
            );
        });
        group.bench_function(format!("cnot/{n_terms}"), |b| {
            b.iter_batched_ref(
                || state.clone(),
                |s| s.cnot([0, 1]),
                criterion::BatchSize::SmallInput,
            );
        });
    }

    group.finish();
}

criterion_group! {
    name = benches;
    config = configure();
    targets = bench_clifford_single, bench_clifford_two, bench_rotations,
              bench_noise, bench_truncation, bench_trace, bench_pauli_word,
              bench_scaling
}
criterion_main!(benches);
