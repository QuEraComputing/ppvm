// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_conformance_2::mixture::{assert_snapshots_close, new, new_snapshot, old, old_snapshot};
use ppvm_traits::traits::{
    Clifford as OldClifford, CliffordBatch as OldCliffordBatch,
    CliffordExtensions as OldExtensions, CliffordExtensionsBatch as OldExtensionsBatch,
    CorrelatedLossChannel as OldCorrelatedLoss, Depolarizing as OldDepolarizing,
    LossChannel as OldLoss, PauliError as OldPauliError, Reset as OldReset,
    ResetLossChannel as OldResetLoss, RotationOne as OldRotation, RotationTwo as OldRotationTwo,
    TGate as OldTGate, TwoQubitPauliError as OldPauli2, U3Gate as OldU3,
};
use ppvm_traits_2::{
    Clifford as NewClifford, CliffordBatch as NewCliffordBatch,
    CliffordExtensions as NewExtensions, CliffordExtensionsBatch as NewExtensionsBatch,
    CorrelatedLossChannel as NewCorrelatedLoss, Depolarizing as NewDepolarizing,
    LossChannel as NewLoss, PauliError as NewPauliError, Reset as NewReset,
    ResetLossChannel as NewResetLoss, RotXY as NewRotXY, RotationOne as NewRotation,
    RotationTwo as NewRotationTwo, TGate as NewTGate, TwoQubitPauliError as NewPauli2,
    U3Gate as NewU3,
};

fn assert_same(old: &ppvm_conformance_2::mixture::Old, new: &ppvm_conformance_2::mixture::New) {
    assert_snapshots_close(old_snapshot(old), new_snapshot(new));
}

#[test]
fn constructor_strict_cutoff_matches_old() {
    for (cutoff, expected_len) in [(0.999, 1), (1.0, 0), (2.0, 0)] {
        let old = old(1, cutoff);
        let new = new(1, cutoff);
        assert_eq!(old.len(), expected_len, "old cutoff {cutoff}");
        assert_eq!(new.len(), expected_len, "new cutoff {cutoff}");
        assert_same(&old, &new);
    }
}

#[test]
fn deterministic_gate_surface_matches_structurally() {
    let (mut old, mut new) = (old(7, 1e-14), new(7, 1e-14));
    old.h(0);
    new.h(0);
    old.x(4);
    new.x(4);
    old.y(5);
    new.y(5);
    old.z(6);
    new.z(6);
    old.cnot(0, 1);
    new.cnot(0, 1);
    old.cz(1, 2);
    new.cz(1, 2);
    old.s(1);
    new.s(1);
    old.sqrt_x(2);
    new.sqrt_x(2);
    old.s_dag(3);
    new.s_dag(3);
    old.sqrt_x_dag(4);
    new.sqrt_x_dag(4);
    old.sqrt_y(5);
    new.sqrt_y(5);
    old.sqrt_y_dag(6);
    new.sqrt_y_dag(6);
    old.cy(2, 3);
    new.cy(2, 3);
    old.h_many(&[7, 8]);
    new.h_many(&[7, 8]);
    old.x_many(&[7]);
    new.x_many(&[7]);
    old.y_many(&[8]);
    new.y_many(&[8]);
    old.z_many(&[9]);
    new.z_many(&[9]);
    old.s_many(&[10]);
    new.s_many(&[10]);
    old.cnot_many(&[(7, 9), (8, 10)]);
    new.cnot_many(&[(7, 9), (8, 10)]);
    old.cz_many(&[(7, 10)]);
    new.cz_many(&[(7, 10)]);
    old.s_dag_many(&[7]);
    new.s_dag_many(&[7]);
    old.sqrt_x_many(&[8]);
    new.sqrt_x_many(&[8]);
    old.sqrt_x_dag_many(&[9]);
    new.sqrt_x_dag_many(&[9]);
    old.sqrt_y_many(&[9, 10]);
    new.sqrt_y_many(&[9, 10]);
    old.sqrt_y_dag_many(&[10]);
    new.sqrt_y_dag_many(&[10]);
    old.cy_many(&[(7, 8)]);
    new.cy_many(&[(7, 8)]);
    old.t(0);
    new.t(0);
    old.t_dag(1);
    new.t_dag(1);
    old.rx(1, 0.37);
    new.rx(1, 0.37);
    old.ry(2, -0.22);
    new.ry(2, -0.22);
    old.rz(3, 0.51);
    new.rz(3, 0.51);
    old.rotate_2([1, 0], [0, 1], 4, 5, 0.19);
    new.rotate_2([1, 0], [0, 1], 4, 5, 0.19);
    old.rz(6, -0.3);
    old.rx(6, -0.4);
    old.rz(6, 0.3);
    new.r(6, 0.3, -0.4);
    old.u3(2, 0.2, -0.4, 0.7);
    new.u3(2, 0.2, -0.4, 0.7);
    assert_same(&old, &new);
}

#[test]
fn analytic_measurement_case_a_and_case_b_match() {
    let (mut old, mut new) = (old(11, 1e-14), new(11, 1e-14));
    old.h(0);
    new.h(0);
    let old_a = old.measure(0);
    let new_a = new.measure(0);
    assert_eq!(old_a, new_a);
    assert_same(&old, &new);
    let old_b = old.measure(0);
    let new_b = new.measure(0);
    assert_eq!(old_b, new_b);
    assert_same(&old, &new);
}

#[test]
fn reset_coalesces_both_outcomes() {
    let (mut old, mut new) = (old(13, 1e-14), new(13, 1e-14));
    old.h(0);
    new.h(0);
    old.reset(0);
    new.reset(0, &mut ppvm_conformance_2::analytic_rng());
    assert_eq!(old.len(), 1);
    assert_eq!(new.len(), 1);
    assert_same(&old, &new);
}

#[test]
fn noise_loss_and_reset_loss_match() {
    let (mut old, mut new) = (old(17, 1e-14), new(17, 1e-14));
    old.h(0);
    new.h(0);
    old.pauli_error(0, [0.11, 0.07, 0.03]);
    new.pauli_error(
        0,
        [0.11, 0.07, 0.03],
        &mut ppvm_conformance_2::analytic_rng(),
    );
    old.depolarize1(1, 0.2);
    new.depolarize1(1, 0.2, &mut ppvm_conformance_2::analytic_rng());
    let p = std::array::from_fn(|i| (i + 1) as f64 / 2000.0);
    old.two_qubit_pauli_error(2, 3, p);
    new.two_qubit_pauli_error(2, 3, p, &mut ppvm_conformance_2::analytic_rng());
    old.loss_channel(4, 0.3);
    new.loss_channel(4, 0.3, &mut ppvm_conformance_2::analytic_rng());
    old.correlated_loss_channel(5, 6, [0.08, 0.12, 0.25]);
    new.correlated_loss_channel(
        5,
        6,
        [0.08, 0.12, 0.25],
        &mut ppvm_conformance_2::analytic_rng(),
    );
    old.reset_loss_channel(4);
    new.reset_loss_channel(4);
    assert_same(&old, &new);
}

#[test]
fn strict_sum_cutoff_boundary_matches() {
    let (mut old, mut new) = (old(19, 0.25), new(19, 0.25));
    old.loss_channel(0, 0.25);
    new.loss_channel(0, 0.25, &mut ppvm_conformance_2::analytic_rng());
    assert_eq!(old.len(), 1);
    assert_eq!(new.len(), 1);
    assert_same(&old, &new);
}

#[test]
fn wide_bitstring_indices_are_supported() {
    type Wide = ppvm_tableau_2::GeneralizedTableauMixture<u128>;
    let mut mixture = Wide::new_with_seed(80, 1e-12, 1e-14, 23);
    mixture.h(79);
    mixture.t(79);
    mixture.loss_channel(70, 0.2, &mut ppvm_conformance_2::analytic_rng());
    assert!(!mixture.is_empty());
}

#[test]
fn seeded_sampler_matches_old_and_is_reproducible() {
    let build = || {
        let (mut old, mut new) = (old(29, 1e-14), new(29, 1e-14));
        old.h(0);
        new.h(0);
        old.cnot(0, 1);
        new.cnot(0, 1);
        old.loss_channel(2, 0.35);
        new.loss_channel(2, 0.35, &mut ppvm_conformance_2::analytic_rng());
        old.two_qubit_pauli_error(3, 4, [0.0; 15]);
        new.two_qubit_pauli_error(3, 4, [0.0; 15], &mut ppvm_conformance_2::analytic_rng());
        old.correlated_loss_channel(5, 6, [0.0; 3]);
        new.correlated_loss_channel(5, 6, [0.0; 3], &mut ppvm_conformance_2::analytic_rng());
        let old_shots = old.sampler().sample_shots_serial(64);
        let new_shots = new.sampler().sample_shots_serial(64);
        assert_eq!(old_shots, new_shots);
        new_shots
    };
    assert_eq!(build(), build());
}
