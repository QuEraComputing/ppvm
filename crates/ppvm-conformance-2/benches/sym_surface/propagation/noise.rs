// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_conformance_2::sym::{NewTerm, OldTerm};
use ppvm_traits::traits::{
    Depolarizing as OldDep, Depolarizing2 as OldDep2, PauliError as OldError,
    TwoQubitPauliError as OldError2,
};
use ppvm_traits_2::{
    Depolarizing as NewDep, Depolarizing2 as NewDep2, PauliError as NewError,
    TwoQubitPauliError as NewError2,
};

use super::paired_args;

pub(super) fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("sym/surface/propagation/noise");
    let targets = &[0, 2, 4][..];
    let pairs = &[(0, 1), (2, 3)][..];

    macro_rules! scalar {
        ($name:literal, $old:expr, $new:expr) => {
            paired_args(
                &mut group,
                $name,
                OldTerm::from(0.01),
                NewTerm::from(0.01),
                $old,
                $new,
            );
        };
    }

    let old_p = [
        OldTerm::from(0.01),
        OldTerm::from(0.02),
        OldTerm::from(0.03),
    ];
    let new_p = [
        NewTerm::from(0.01),
        NewTerm::from(0.02),
        NewTerm::from(0.03),
    ];
    paired_args(
        &mut group,
        "pauli_error",
        old_p.clone(),
        new_p.clone(),
        |s, p| OldError::pauli_error(s, 0, p),
        |s, p| NewError::pauli_error(s, 0, p, &mut ppvm_conformance_2::analytic_rng()),
    );
    paired_args(
        &mut group,
        "batch_pauli_error",
        old_p,
        new_p,
        |s, p| OldError::pauli_error_many(s, targets, p),
        |s, p| NewError::pauli_error_many(s, targets, p, &mut ppvm_conformance_2::analytic_rng()),
    );
    scalar!("x_error", |s, p| OldError::x_error(s, 0, p), |s, p| {
        NewError::x_error(s, 0, p, &mut ppvm_conformance_2::analytic_rng())
    });
    scalar!("y_error", |s, p| OldError::y_error(s, 0, p), |s, p| {
        NewError::y_error(s, 0, p, &mut ppvm_conformance_2::analytic_rng())
    });
    scalar!("z_error", |s, p| OldError::z_error(s, 0, p), |s, p| {
        NewError::z_error(s, 0, p, &mut ppvm_conformance_2::analytic_rng())
    });
    scalar!(
        "batch_x_error",
        |s, p| OldError::x_error_many(s, targets, p),
        |s, p| { NewError::x_error_many(s, targets, p, &mut ppvm_conformance_2::analytic_rng(),) }
    );
    scalar!(
        "batch_y_error",
        |s, p| OldError::y_error_many(s, targets, p),
        |s, p| { NewError::y_error_many(s, targets, p, &mut ppvm_conformance_2::analytic_rng(),) }
    );
    scalar!(
        "batch_z_error",
        |s, p| OldError::z_error_many(s, targets, p),
        |s, p| { NewError::z_error_many(s, targets, p, &mut ppvm_conformance_2::analytic_rng(),) }
    );

    let old_p2 = std::array::from_fn(|i| OldTerm::from((i + 1) as f64 * 1e-4));
    let new_p2 = std::array::from_fn(|i| NewTerm::from((i + 1) as f64 * 1e-4));
    paired_args(
        &mut group,
        "two_qubit_pauli_error",
        old_p2.clone(),
        new_p2.clone(),
        |s, p| OldError2::two_qubit_pauli_error(s, 0, 1, p),
        |s, p| {
            NewError2::two_qubit_pauli_error(s, 0, 1, p, &mut ppvm_conformance_2::analytic_rng())
        },
    );
    paired_args(
        &mut group,
        "batch_two_qubit_pauli_error",
        old_p2,
        new_p2,
        |s, p| OldError2::two_qubit_pauli_error_many(s, pairs, p),
        |s, p| {
            NewError2::two_qubit_pauli_error_many(
                s,
                pairs,
                p,
                &mut ppvm_conformance_2::analytic_rng(),
            )
        },
    );

    scalar!(
        "depolarize1",
        |s, p| OldDep::depolarize1(s, 0, p),
        |s, p| NewDep::depolarize1(s, 0, p, &mut ppvm_conformance_2::analytic_rng())
    );
    scalar!(
        "batch_depolarize1",
        |s, p| OldDep::depolarize1_many(s, targets, p),
        |s, p| { NewDep::depolarize1_many(s, targets, p, &mut ppvm_conformance_2::analytic_rng()) }
    );
    scalar!(
        "depolarize2",
        |s, p| OldDep2::depolarize2(s, 0, 1, p),
        |s, p| NewDep2::depolarize2(s, 0, 1, p, &mut ppvm_conformance_2::analytic_rng())
    );
    scalar!(
        "batch_depolarize2",
        |s, p| OldDep2::depolarize2_many(s, pairs, p),
        |s, p| { NewDep2::depolarize2_many(s, pairs, p, &mut ppvm_conformance_2::analytic_rng(),) }
    );
    group.finish();
}
