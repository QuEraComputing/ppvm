// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hint::black_box;

use criterion::Criterion;
use ppvm_conformance_2::mixture::{New, Old, new, old};
use ppvm_traits::traits::{
    Clifford as OldClifford, CorrelatedLossChannel as OldCorrelatedLoss,
    Depolarizing as OldDepolarizing, LossChannel as OldLoss, TGate as OldTGate,
};
use ppvm_traits_2::{
    Clifford as NewClifford, CorrelatedLossChannel as NewCorrelatedLoss,
    Depolarizing as NewDepolarizing, LossChannel as NewLoss, TGate as NewTGate,
};

use super::support::assert_same;

pub fn register(c: &mut Criterion) {
    let old_check = build_old();
    let new_check = build_new();
    assert_same(&old_check, &new_check);
    assert_eq!(old_check.len(), new_check.len());

    let mut group = c.benchmark_group("mixture/integration/noisy_build");
    group.bench_function("old", |b| b.iter(|| black_box(build_old())));
    group.bench_function("new", |b| b.iter(|| black_box(build_new())));
    group.finish();
}

fn build_old() -> Old {
    let mut state = old(41, 1e-5);
    for layer in 0..3 {
        for qubit in 0..6 {
            state.h(qubit);
            state.t(qubit);
            state.depolarize1(qubit, 0.002);
        }
        for qubit in 0..5 {
            state.cnot(qubit, qubit + 1);
        }
        state.correlated_loss_channel(layer, layer + 1, [0.001, 0.002, 0.003]);
        state.loss_channel(layer + 8, 0.001);
    }
    state
}

fn build_new() -> New {
    let mut state = new(41, 1e-5);
    let mut rng = ppvm_conformance_2::analytic_rng();
    for layer in 0..3 {
        for qubit in 0..6 {
            state.h(qubit);
            state.t(qubit);
            state.depolarize1(qubit, 0.002, &mut rng);
        }
        for qubit in 0..5 {
            state.cnot(qubit, qubit + 1);
        }
        state.correlated_loss_channel(layer, layer + 1, [0.001, 0.002, 0.003], &mut rng);
        state.loss_channel(layer + 8, 0.001, &mut rng);
    }
    state
}
