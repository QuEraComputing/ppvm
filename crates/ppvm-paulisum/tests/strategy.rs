// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_paulisum::{
    prelude::*,
    strategy::{CoefficientThreshold, CombinedStrategy, MaxPauliWeight},
};

#[test]
fn test_coefficient_threshold() {
    let strat = CoefficientThreshold(1e-4);
    let mut state: PauliSum<config::indexmap::ByteFxHashF64<4, CoefficientThreshold>> =
        PauliSum::builder().n_qubits(2).strategy(strat).build();

    state += ("ZZ", 1.0);
    state += ("XX", 1e-3);
    state += ("YY", 1e-5);
    state += ("XZ", -0.5);
    state.truncate();

    let mut target_state: PauliSum<config::indexmap::ByteFxHashF64<4, CoefficientThreshold>> =
        PauliSum::builder().strategy(strat).n_qubits(2).build();
    target_state += ("ZZ", 1.0);
    target_state += ("XX", 1e-3);
    target_state += ("XZ", -0.5);

    assert_eq!(state, target_state);
}

#[test]
fn test_pauli_weight() {
    let strat = MaxPauliWeight(2);

    let mut state: PauliSum<config::indexmap::ByteFxHashF64<4, MaxPauliWeight>> =
        PauliSum::builder().n_qubits(3).strategy(strat).build();

    state += ("XXX", 1.0);
    state += ("ZZI", 1e-4);
    state += ("YZY", -10.0);
    state += ("IXX", -0.1);

    state.truncate();

    let mut target_state: PauliSum<config::indexmap::ByteFxHashF64<4, MaxPauliWeight>> =
        PauliSum::builder().n_qubits(3).strategy(strat).build();

    target_state += ("ZZI", 1e-4);
    target_state += ("IXX", -0.1);

    assert_eq!(state, target_state);
}

#[test]
fn test_combined_strategy() {
    let cutoff_strategy = CoefficientThreshold(1e-4);
    let max_weight_strategy = MaxPauliWeight(2);
    let strat = CombinedStrategy(cutoff_strategy, max_weight_strategy);

    let mut state: PauliSum<
        config::indexmap::ByteFxHashF64<4, CombinedStrategy<CoefficientThreshold, MaxPauliWeight>>,
    > = PauliSum::builder()
        .n_qubits(3)
        .capacity(8)
        .strategy(strat)
        .build();

    state += ("ZZI", 1.0);
    state += ("XIX", 1e-3);
    state += ("IYY", 1e-5);
    state += ("XIZ", -0.5);

    state += ("XXX", 1.0);
    state += ("ZIZ", 1e-4);
    state += ("YZY", -10.0);
    state += ("IXX", -0.1);

    state.truncate();

    let mut target_state: PauliSum<
        config::indexmap::ByteFxHashF64<4, CombinedStrategy<CoefficientThreshold, MaxPauliWeight>>,
    > = PauliSum::builder().n_qubits(3).strategy(strat).build();

    target_state += ("ZZI", 1.0);
    target_state += ("XIX", 1e-3);
    target_state += ("XIZ", -0.5);

    target_state += ("ZIZ", 1e-4);
    target_state += ("IXX", -0.1);

    assert_eq!(state, target_state);
}
