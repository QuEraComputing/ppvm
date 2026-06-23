// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Integration test: `measure_many` (measure a chosen set of qubit indices,
//! reusing one `MeasureScratch`) must agree — in both outcomes and the
//! measurement record — with the established per-qubit `measure` / `measure_all`
//! paths. `measure_all` is exactly the `0..n_qubits` special case.

use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::measure_all::LossyMeasureAll;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<Byte8F64<2>, u128>;

/// Deterministic, entangled, non-Clifford starting state (no measurement yet).
fn build_state(n: usize, seed: u64) -> Tab {
    let mut tab: Tab = GeneralizedTableau::new_with_seed(n, 1e-10, seed);
    for q in 0..n {
        tab.h(q);
    }
    for q in 0..n - 1 {
        tab.cz(q, q + 1);
    }
    // A couple of T gates exercise the case-a (Z-not-a-stabilizer) measurement path.
    tab.t(0);
    tab.t(n / 2);
    tab
}

#[test]
fn measure_many_all_indices_matches_measure_all() {
    let n = 12;
    for seed in 0..25 {
        let batch = build_state(n, seed).measure_many(&(0..n).collect::<Vec<_>>());
        let all = build_state(n, seed).measure_all();
        assert_eq!(batch, all, "seed={seed}");
    }
}

#[test]
fn measure_many_matches_individual_measure_in_order() {
    let n = 10;
    let order = [3usize, 0, 7, 1, 9, 2, 8, 4, 6, 5];
    for seed in 0..25 {
        let batch = build_state(n, seed).measure_many(&order);
        let mut individual_tab = build_state(n, seed);
        let individual: Vec<Option<bool>> =
            order.iter().map(|&q| individual_tab.measure(q)).collect();
        assert_eq!(batch, individual, "seed={seed}");
    }
}

#[test]
fn measure_many_subset_matches_individual() {
    let n = 10;
    let subset = [2usize, 5, 8];
    for seed in 0..25 {
        let batch = build_state(n, seed).measure_many(&subset);
        let mut individual_tab = build_state(n, seed);
        let individual: Vec<Option<bool>> =
            subset.iter().map(|&q| individual_tab.measure(q)).collect();
        assert_eq!(batch, individual, "seed={seed}");
    }
}

#[test]
fn measure_many_empty_returns_empty() {
    let mut tab = build_state(6, 1);
    assert_eq!(tab.measure_many(&[]), Vec::<Option<bool>>::new());
}

#[test]
fn measure_many_lost_qubit_is_none() {
    let mut tab: Tab = GeneralizedTableau::new_with_seed(4, 1e-10, 7);
    tab.h(0);
    tab.cz(0, 1);
    tab.loss_channel(1, 1.0); // p = 1.0 -> qubit 1 is always lost
    let res = tab.measure_many(&[0, 1, 2, 3]);
    assert_eq!(res.len(), 4);
    assert_eq!(res[1], None, "lost qubit must measure as None");
    assert!(
        res[0].is_some() && res[2].is_some() && res[3].is_some(),
        "non-lost qubits must yield a result"
    );
}

/// `measure_many` must leave the measurement record byte-for-byte identical to
/// the per-qubit `measure` loop — including pushing `None` for lost qubits, so
/// downstream `rec[-k]` references (the Stim executor) don't shift. Covers a
/// state with loss, which the outcome-only tests above don't exercise.
#[test]
fn measure_many_record_matches_individual_with_loss() {
    let order = [0usize, 1, 2, 3, 4, 5];

    fn build_lossy(seed: u64) -> Tab {
        let mut tab: Tab = GeneralizedTableau::new_with_seed(6, 1e-10, seed);
        for q in 0..6 {
            tab.h(q);
        }
        for q in 0..5 {
            tab.cz(q, q + 1);
        }
        tab.loss_channel(2, 1.0); // qubit 2 always lost
        tab.loss_channel(4, 1.0); // qubit 4 always lost
        tab
    }

    for seed in 0..25 {
        let mut batch_tab = build_lossy(seed);
        let mut indiv_tab = build_lossy(seed);

        let batch = batch_tab.measure_many(&order);
        let individual: Vec<Option<bool>> = order.iter().map(|&q| indiv_tab.measure(q)).collect();

        assert_eq!(batch, individual, "returned results differ at seed={seed}");
        assert_eq!(
            batch_tab.current_measurement_record(),
            indiv_tab.current_measurement_record(),
            "measurement record differs at seed={seed}"
        );
        // Lost qubits 2 and 4 must each have recorded a `None`.
        assert_eq!(batch_tab.current_measurement_record().len(), order.len());
        assert_eq!(batch_tab.current_measurement_record()[2], None);
        assert_eq!(batch_tab.current_measurement_record()[4], None);
    }
}
