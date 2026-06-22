// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Micro-benchmark for the cost of `MaxPauliWeight` truncation on a
//! pre-built `PauliSum`. Compares three strategies:
//!   - `MaxPauliWeight(w)` alone (single retain pass, weight check only)
//!   - `CoefficientThreshold(eps)` alone (single retain pass, coeff check only)
//!   - `CombinedStrategy(threshold, weight)` (two sequential retain passes)
//!
//! Mirrored on the Julia side by
//! `julia-benchmarks/benches/truncation-weight.jl`, which uses
//! PauliPropagation.jl's `truncate!` — a *single* walk with a combined
//! predicate. The point of comparison is the two-passes-vs-one structural
//! difference between Rust's `CombinedStrategy` and PP.jl's `truncate!`.

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_pauli_sum::config::fxhash::ByteF64;
use ppvm_pauli_sum::prelude::*;
use ppvm_pauli_sum::strategy::{CoefficientThreshold, CombinedStrategy, MaxPauliWeight};

const N_QUBITS: usize = 128;
const N: usize = 16;
const N_TERMS: usize = 1000;
const COEFF_EPS: f64 = 1e-12;

type Word = <ByteF64<N, NoStrategy> as Config>::PauliWordType;

fn configure() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50)
}

/// Build `N_TERMS` Pauli terms each with exactly `target_weight`
/// non-identity slots (clamped to `N_QUBITS`). Positions are spread by
/// stride to fit `target_weight` distinct slots; Pauli choices depend
/// on `k` so terms are mostly distinct.
fn make_terms(target_weight: usize) -> Vec<(Word, f64)> {
    let paulis = [Pauli::X, Pauli::Y, Pauli::Z];
    let weight = target_weight.min(N_QUBITS);
    let stride = (N_QUBITS / weight.max(1)).max(1);
    let mut terms = Vec::with_capacity(N_TERMS);
    for k in 0..N_TERMS {
        let mut word = Word::new(N_QUBITS);
        for w in 0..weight {
            let pos = (k + w * stride) % N_QUBITS;
            word.set(pos, paulis[(k.wrapping_mul(31) + w) % 3]);
        }
        terms.push((word, 1.0 / (k as f64 + 1.0)));
    }
    terms
}

fn build_state<St>(strat: St, terms: &[(Word, f64)]) -> PauliSum<ByteF64<N, St>>
where
    St: Strategy + Copy,
{
    let mut state: PauliSum<ByteF64<N, St>> = PauliSum::builder()
        .n_qubits(N_QUBITS)
        .capacity(terms.len() * 2)
        .strategy(strat)
        .build();
    for (w, c) in terms {
        state += (w.clone(), *c);
    }
    state
}

const PROFILES: &[(&str, usize)] = &[
    ("w3", 3),
    ("w50", 50),
    ("w120", 120),
];

const CUTOFFS: &[(&str, usize)] = &[
    ("10", 10),
    ("100", 100),
    ("1000", 1000),
    ("max", usize::MAX),
];

fn bench_max_weight_alone(c: &mut Criterion) {
    let mut group = c.benchmark_group("max-weight-only");
    for &(profile, weight) in PROFILES {
        let terms = make_terms(weight);
        for &(cutoff_name, cutoff) in CUTOFFS {
            let state = build_state(MaxPauliWeight(cutoff), &terms);
            println!(
                "[max-weight-only/{profile}/cut-{cutoff_name}] state.len()={}",
                state.len()
            );
            group.bench_function(format!("{profile}/cut-{cutoff_name}"), |b| {
                b.iter_batched_ref(
                    || state.clone(),
                    |s| s.truncate(),
                    criterion::BatchSize::SmallInput,
                );
            });
        }
    }
    group.finish();
}

fn bench_coeff_threshold_alone(c: &mut Criterion) {
    let mut group = c.benchmark_group("coeff-threshold-only");
    for &(profile, weight) in PROFILES {
        let terms = make_terms(weight);
        let state = build_state(CoefficientThreshold(COEFF_EPS), &terms);
        println!(
            "[coeff-threshold-only/{profile}] state.len()={}",
            state.len()
        );
        group.bench_function(profile, |b| {
            b.iter_batched_ref(
                || state.clone(),
                |s| s.truncate(),
                criterion::BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn bench_combined(c: &mut Criterion) {
    let mut group = c.benchmark_group("combined");
    for &(profile, weight) in PROFILES {
        let terms = make_terms(weight);
        for &(cutoff_name, cutoff) in CUTOFFS {
            let strat = CombinedStrategy(CoefficientThreshold(COEFF_EPS), MaxPauliWeight(cutoff));
            let state = build_state(strat, &terms);
            group.bench_function(format!("{profile}/cut-{cutoff_name}"), |b| {
                b.iter_batched_ref(
                    || state.clone(),
                    |s| s.truncate(),
                    criterion::BatchSize::SmallInput,
                );
            });
        }
    }
    group.finish();
}

/// Clone-only baseline so the user can subtract it from the truncate
/// numbers if desired. `iter_batched_ref` excludes setup time, so this
/// shouldn't be necessary — but worth confirming.
fn bench_clone_baseline(c: &mut Criterion) {
    let mut group = c.benchmark_group("clone-baseline");
    for &(profile, weight) in PROFILES {
        let terms = make_terms(weight);
        let state = build_state(NoStrategy, &terms);
        group.bench_function(profile, |b| {
            b.iter(|| state.clone());
        });
    }
    group.finish();
}

criterion_group! {
    name = benches;
    config = configure();
    targets = bench_max_weight_alone,
              bench_coeff_threshold_alone,
              bench_combined,
              bench_clone_baseline
}
criterion_main!(benches);
