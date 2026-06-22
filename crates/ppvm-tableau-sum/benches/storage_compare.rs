// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Compares `VecStorage` and `MapStorage` on an msd-noisy-shaped workload:
//! every gate is followed by 1-4 noise channels on the qubits it touched.
//! Mixes Clifford, T, loss, and depolarize so `insert_or_merge_batch` runs
//! against a non-trivial branch count.

use std::time::Duration;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use num::complex::Complex64;
use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau_sum::{
    data::GeneralizedTableauSum,
    storage::{EntryStore, map::MapStorage, vec::VecStorage},
};
use ppvm_traits::traits::{Clifford, CliffordExtensions, Depolarizing, LossChannel, TGate};

type Cfg = Byte8F64<2>;
type Idx = u128;
type Coef = Vec<(Complex64, Idx)>;

type VecBackedSum = GeneralizedTableauSum<Cfg, Idx, Coef, VecStorage<Cfg, Idx, Coef>>;
type MapBackedSum = GeneralizedTableauSum<Cfg, Idx, Coef, MapStorage<Cfg, Idx, Coef>>;

const N_QUBITS: usize = 17;
const P_LOSS: f64 = 1e-4;
const P_DEPOL: f64 = 1e-4;
const SUM_CUTOFF: f64 = 1e-7;
const COEFF_THRESHOLD: f64 = 1e-10;
const SEED: u64 = 42;

/// Gate + 1-4 noise channels per op. Mirrors the msd-noisy pattern.
fn apply_circuit<S>(tab: &mut GeneralizedTableauSum<Cfg, Idx, Coef, S>)
where
    S: EntryStore<Cfg, Idx, Coef>,
    GeneralizedTableauSum<Cfg, Idx, Coef, S>:
        Clifford + CliffordExtensions + TGate<Cfg> + LossChannel<Cfg> + Depolarizing<Cfg>,
{
    // Single-qubit Cliffords, two noise channels each.
    for q in 0..N_QUBITS {
        tab.sqrt_y(q);
        tab.loss_channel(q, P_LOSS);
        tab.depolarize(q, P_DEPOL);
    }

    // Non-Clifford (T) on a few qubits, two noise channels each.
    for q in [0, 7, 12] {
        tab.t(q);
        tab.loss_channel(q, P_LOSS);
        tab.depolarize(q, P_DEPOL);
    }

    // Two-qubit Cliffords, four noise channels each (loss+depolarize on both).
    for [i, j] in [
        [1, 3],
        [7, 10],
        [12, 14],
        [13, 16],
        [4, 7],
        [8, 10],
        [11, 14],
        [15, 16],
    ] {
        tab.cz(i, j);
        tab.loss_channel(i, P_LOSS);
        tab.loss_channel(j, P_LOSS);
        tab.depolarize(i, P_DEPOL);
        tab.depolarize(j, P_DEPOL);
    }

    // A trailing single-qubit Clifford layer with one noise channel each.
    for q in [0, 2, 5, 6, 8, 10, 12] {
        tab.sqrt_y_adj(q);
        tab.loss_channel(q, P_LOSS);
    }
}

fn vec_sum() -> VecBackedSum {
    GeneralizedTableauSum::new_with_seed(N_QUBITS, COEFF_THRESHOLD, SUM_CUTOFF, SEED)
}

fn map_sum() -> MapBackedSum {
    GeneralizedTableauSum::new_with_seed(N_QUBITS, COEFF_THRESHOLD, SUM_CUTOFF, SEED)
}

fn storage_compare(c: &mut Criterion) {
    let mut group = c.benchmark_group("storage-compare");

    group.bench_function("vec", |b| {
        b.iter_batched(
            vec_sum,
            |mut tab| {
                apply_circuit(&mut tab);
                tab
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("map", |b| {
        b.iter_batched(
            map_sum,
            |mut tab| {
                apply_circuit(&mut tab);
                tab
            },
            BatchSize::SmallInput,
        );
    });

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_millis(1500))
        .sample_size(30);
    targets = storage_compare
}
criterion_main!(benches);
