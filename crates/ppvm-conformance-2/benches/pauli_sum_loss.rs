// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Same-build old/new performance gate for the lossy-atom `PauliSum` workload.

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::time::Duration;

use ppvm_lossy_pauli_word_2::LossyPauliWord as NewWord;
use ppvm_pauli_sum::config;
use ppvm_pauli_sum::sum::PauliSum as OldSumT;
use ppvm_pauli_sum_2::{HashMapStore, NoPolicy, Sum};
use ppvm_pauli_word::loss::LossyPauliWord as OldWord;
use ppvm_traits::traits::{
    Clifford as OldClifford, CorrelatedLossChannel as OldCorrelatedLossChannel,
    LossChannel as OldLossChannel, NoStrategy, ResetLossChannel as OldResetLossChannel,
    RotationOne as OldRotationOne,
};
use ppvm_traits_2::{
    Clifford as NewClifford, CorrelatedLossChannel as NewCorrelatedLossChannel,
    LossChannel as NewLossChannel, ResetLossChannel as NewResetLossChannel,
    RotationOne as NewRotationOne,
};

const N: usize = 12;
const CAPACITY: usize = 1 << N;

type OldConfig = config::fxhash::Byte<2, f64, NoStrategy, OldWord<[u8; 2], fxhash::FxBuildHasher>>;
type OldSum = OldSumT<OldConfig>;
type NewKey = NewWord<[u8; 2]>;
type NewSum = Sum<HashMapStore<NewKey, f64>, NoPolicy>;

fn old_support(sum: &OldSum) -> Vec<(String, f64)> {
    let mut support: Vec<_> = sum
        .data()
        .iter()
        .map(|(key, coeff)| (key.to_string(), *coeff))
        .collect();
    support.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    support
}

fn new_support(sum: &NewSum) -> Vec<(String, f64)> {
    let mut support: Vec<_> = sum
        .iter()
        .map(|(key, coeff)| (key.to_string(), coeff))
        .collect();
    support.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    support
}

#[track_caller]
fn assert_same(old: &OldSum, new: &NewSum, label: &str) {
    let old = old_support(old);
    let new = new_support(new);
    assert_eq!(
        old.iter().map(|(key, _)| key).collect::<Vec<_>>(),
        new.iter().map(|(key, _)| key).collect::<Vec<_>>(),
        "[{label}] canonical lossy keys differ"
    );
    for ((key, old_coeff), (_, new_coeff)) in old.iter().zip(&new) {
        // The n=12 workload merges thousands of branches in hash iteration
        // order, so the two engines can differ by floating-point summation ulps.
        let tolerance = 1e-12_f64.max(old_coeff.abs() * 1e-12);
        assert!(
            (old_coeff - new_coeff).abs() <= tolerance,
            "[{label}] coefficient at {key} differs: old={old_coeff}, new={new_coeff}, tol={tolerance}"
        );
    }
}

trait ApplyLoss {
    fn apply_loss(&mut self, q: usize, p: f64);
    fn apply_correlated_loss(&mut self, a: usize, b: usize, p: [f64; 3]);
}

impl ApplyLoss for OldSum {
    fn apply_loss(&mut self, q: usize, p: f64) {
        self.loss_channel(q, p);
    }

    fn apply_correlated_loss(&mut self, a: usize, b: usize, p: [f64; 3]) {
        self.correlated_loss_channel(a, b, p);
    }
}

impl ApplyLoss for NewSum {
    fn apply_loss(&mut self, q: usize, p: f64) {
        self.loss_channel(q, p, &mut ppvm_conformance_2::analytic_rng());
    }

    fn apply_correlated_loss(&mut self, a: usize, b: usize, p: [f64; 3]) {
        self.correlated_loss_channel(a, b, p, &mut ppvm_conformance_2::analytic_rng());
    }
}

fn old_seed() -> OldSum {
    let mut sum = OldSum::builder().n_qubits(N).capacity(CAPACITY).build();
    sum += ("Z".repeat(N).as_str(), 1.0);
    sum
}

fn new_seed() -> NewSum {
    let mut sum = Sum::with_capacity(N, NoPolicy, CAPACITY);
    sum += (NewKey::from("Z".repeat(N).as_str()), 1.0);
    sum
}

fn old_expanded() -> OldSum {
    let mut sum = old_seed();
    for q in 0..N {
        sum.reset_loss_channel(q);
    }
    sum
}

fn new_expanded() -> NewSum {
    let mut sum = new_seed();
    for q in 0..N {
        sum.reset_loss_channel(q);
    }
    sum
}

macro_rules! workload {
    ($sum:ident) => {{
        for q in 0..N {
            $sum.reset_loss_channel(q);
            $sum.apply_loss(q, 0.01);
        }
        for q in 0..N - 1 {
            $sum.apply_correlated_loss(q, q + 1, [0.002, 0.003, 0.004]);
        }
        for q in (0..N - 1).rev() {
            $sum.cnot(q, q + 1);
        }
        $sum.h(0);
        black_box($sum.len())
    }};
}

fn bench_loss(c: &mut Criterion) {
    let mut old_probe = old_seed();
    let mut new_probe = new_seed();
    workload!(old_probe);
    workload!(new_probe);
    assert_same(&old_probe, &new_probe, "loss interleaved");

    let mut group = c.benchmark_group("pauli_sum/loss_interleaved_n12");
    // One n=12 iteration is intentionally heavy (~1 s on the old engine).
    // Ten samples still provide a stable same-build ratio without turning this
    // integration gate into a multi-minute benchmark.
    group
        .sample_size(10)
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(15));
    let old = old_seed();
    group.bench_function("old", |b| {
        b.iter_batched_ref(|| old.clone(), |sum| workload!(sum), BatchSize::LargeInput)
    });
    let new = new_seed();
    group.bench_function("new", |b| {
        b.iter_batched_ref(|| new.clone(), |sum| workload!(sum), BatchSize::LargeInput)
    });
    group.finish();

    macro_rules! stage {
        ($name:literal, $old:expr, $new:expr, |$sum:ident| $body:block) => {{
            let mut old_probe = $old;
            {
                let $sum = &mut old_probe;
                $body
            }
            let mut new_probe = $new;
            {
                let $sum = &mut new_probe;
                $body
            }
            assert_same(&old_probe, &new_probe, concat!("loss stage ", $name));

            let mut group = c.benchmark_group(concat!("pauli_sum/loss_attrib/", $name));
            group
                .sample_size(10)
                .warm_up_time(Duration::from_secs(1))
                .measurement_time(Duration::from_secs(10));
            let old = $old;
            group.bench_function("old", |b| {
                b.iter_batched_ref(
                    || old.clone(),
                    |$sum| {
                        $body
                        black_box($sum.len())
                    },
                    BatchSize::LargeInput,
                )
            });
            let new = $new;
            group.bench_function("new", |b| {
                b.iter_batched_ref(
                    || new.clone(),
                    |$sum| {
                        $body
                        black_box($sum.len())
                    },
                    BatchSize::LargeInput,
                )
            });
            group.finish();
        }};
    }

    stage!("reset", old_seed(), new_seed(), |sum| {
        for q in 0..N {
            sum.reset_loss_channel(q);
        }
    });
    stage!("loss", old_expanded(), new_expanded(), |sum| {
        for q in 0..N {
            sum.apply_loss(q, 0.01);
        }
    });
    stage!("correlated", old_expanded(), new_expanded(), |sum| {
        for q in 0..N - 1 {
            sum.apply_correlated_loss(q, q + 1, [0.002, 0.003, 0.004]);
        }
    });
    stage!("clifford", old_expanded(), new_expanded(), |sum| {
        for q in (0..N - 1).rev() {
            sum.cnot(q, q + 1);
        }
    });
    stage!("rotation", old_expanded(), new_expanded(), |sum| {
        for q in 0..N {
            sum.rx(q, 0.07);
        }
    });
}

criterion_group!(benches, bench_loss);
criterion_main!(benches);
