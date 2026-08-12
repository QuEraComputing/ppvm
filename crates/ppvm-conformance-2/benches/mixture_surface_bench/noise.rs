// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_conformance_2::mixture::{New, Old};
use ppvm_traits::traits::{
    CorrelatedLossChannel as OldCorrelatedLoss, Depolarizing as OldDepolarizing,
    Depolarizing2 as OldDepolarizing2, LossChannel as OldLoss, PauliError as OldPauliError,
    ResetLossChannel as OldResetLoss, TwoQubitPauliError as OldPauliError2,
};
use ppvm_traits_2::{
    CorrelatedLossChannel as NewCorrelatedLoss, Depolarizing as NewDepolarizing,
    Depolarizing2 as NewDepolarizing2, LossChannel as NewLoss, PauliError as NewPauliError,
    ResetLossChannel as NewResetLoss, TwoQubitPauliError as NewPauliError2,
};

use super::support::{assert_same, bench_mut, branch_pair};

macro_rules! error_one {
    ($c:expr, $old:expr, $new:expr, $method:ident) => {
        bench_mut(
            $c,
            concat!("mixture/noise/", stringify!($method)),
            $old,
            $new,
            |state: &mut Old| state.$method(5, 0.03),
            |state: &mut New| state.$method(5, 0.03, &mut ppvm_conformance_2::analytic_rng()),
        );
    };
}

macro_rules! error_many {
    ($c:expr, $old:expr, $new:expr, $method:ident) => {
        bench_mut(
            $c,
            concat!("mixture/noise/", stringify!($method)),
            $old,
            $new,
            |state: &mut Old| state.$method(&[5, 6, 7], 0.02),
            |state: &mut New| {
                state.$method(&[5, 6, 7], 0.02, &mut ppvm_conformance_2::analytic_rng())
            },
        );
    };
}

pub fn register(c: &mut Criterion) {
    let (old, new) = branch_pair(2);
    let probabilities = [0.011, 0.007, 0.003];
    let probabilities_2 = std::array::from_fn(|i| (i + 1) as f64 / 20_000.0);

    bench_mut(
        c,
        "mixture/noise/pauli_error",
        &old,
        &new,
        move |s: &mut Old| s.pauli_error(5, probabilities),
        move |s: &mut New| s.pauli_error(5, probabilities, &mut ppvm_conformance_2::analytic_rng()),
    );
    error_one!(c, &old, &new, x_error);
    error_one!(c, &old, &new, y_error);
    error_one!(c, &old, &new, z_error);
    bench_mut(
        c,
        "mixture/noise/pauli_error_many",
        &old,
        &new,
        move |s: &mut Old| s.pauli_error_many(&[5, 6, 7], probabilities),
        move |s: &mut New| {
            s.pauli_error_many(
                &[5, 6, 7],
                probabilities,
                &mut ppvm_conformance_2::analytic_rng(),
            )
        },
    );
    error_many!(c, &old, &new, x_error_many);
    error_many!(c, &old, &new, y_error_many);
    error_many!(c, &old, &new, z_error_many);

    bench_mut(
        c,
        "mixture/noise/two_qubit_pauli_error",
        &old,
        &new,
        move |s: &mut Old| s.two_qubit_pauli_error(4, 5, probabilities_2),
        move |s: &mut New| {
            s.two_qubit_pauli_error(
                4,
                5,
                probabilities_2,
                &mut ppvm_conformance_2::analytic_rng(),
            )
        },
    );
    bench_mut(
        c,
        "mixture/noise/two_qubit_pauli_error_many",
        &old,
        &new,
        move |s: &mut Old| s.two_qubit_pauli_error_many(&[(4, 5), (6, 7)], probabilities_2),
        move |s: &mut New| {
            s.two_qubit_pauli_error_many(
                &[(4, 5), (6, 7)],
                probabilities_2,
                &mut ppvm_conformance_2::analytic_rng(),
            )
        },
    );
    bench_mut(
        c,
        "mixture/noise/depolarize1",
        &old,
        &new,
        |s: &mut Old| s.depolarize1(5, 0.03),
        |s: &mut New| s.depolarize1(5, 0.03, &mut ppvm_conformance_2::analytic_rng()),
    );
    bench_mut(
        c,
        "mixture/noise/depolarize1_many",
        &old,
        &new,
        |s: &mut Old| s.depolarize1_many(&[5, 6, 7], 0.02),
        |s: &mut New| s.depolarize1_many(&[5, 6, 7], 0.02, &mut ppvm_conformance_2::analytic_rng()),
    );
    bench_mut(
        c,
        "mixture/noise/depolarize2",
        &old,
        &new,
        |s: &mut Old| s.depolarize2(4, 5, 0.03),
        |s: &mut New| s.depolarize2(4, 5, 0.03, &mut ppvm_conformance_2::analytic_rng()),
    );
    bench_mut(
        c,
        "mixture/noise/depolarize2_many",
        &old,
        &new,
        |s: &mut Old| s.depolarize2_many(&[(4, 5), (6, 7)], 0.02),
        |s: &mut New| {
            s.depolarize2_many(
                &[(4, 5), (6, 7)],
                0.02,
                &mut ppvm_conformance_2::analytic_rng(),
            )
        },
    );

    bench_mut(
        c,
        "mixture/noise/loss_channel",
        &old,
        &new,
        |s: &mut Old| s.loss_channel(5, 0.11),
        |s: &mut New| s.loss_channel(5, 0.11, &mut ppvm_conformance_2::analytic_rng()),
    );
    bench_mut(
        c,
        "mixture/noise/correlated_loss_channel",
        &old,
        &new,
        |s: &mut Old| s.correlated_loss_channel(4, 5, [0.01, 0.02, 0.03]),
        |s: &mut New| {
            s.correlated_loss_channel(
                4,
                5,
                [0.01, 0.02, 0.03],
                &mut ppvm_conformance_2::analytic_rng(),
            )
        },
    );

    let (mut old_lost, mut new_lost) = (old.clone(), new.clone());
    old_lost.loss_channel(5, 0.3);
    new_lost.loss_channel(5, 0.3, &mut ppvm_conformance_2::analytic_rng());
    assert_same(&old_lost, &new_lost);
    bench_mut(
        c,
        "mixture/noise/reset_loss_channel",
        &old_lost,
        &new_lost,
        |s: &mut Old| s.reset_loss_channel(5),
        |s: &mut New| s.reset_loss_channel(5),
    );
}
