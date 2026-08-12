// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_conformance_2::mixture::{New, Old};
use ppvm_traits::traits::{Clifford as OldClifford, LossChannel as OldLoss, Reset as OldReset};
use ppvm_traits_2::{Clifford as NewClifford, LossChannel as NewLoss, Reset as NewReset};

use super::support::{assert_same, bench_mut, bench_output, branch_pair};

macro_rules! reset_one {
    ($c:expr, $old:expr, $new:expr, $method:ident) => {
        bench_mut(
            $c,
            concat!("mixture/reset/", stringify!($method)),
            $old,
            $new,
            |state: &mut Old| state.$method(5),
            |state: &mut New| state.$method(5, &mut ppvm_conformance_2::analytic_rng()),
        );
    };
}

macro_rules! reset_many {
    ($c:expr, $old:expr, $new:expr, $method:ident) => {
        bench_mut(
            $c,
            concat!("mixture/reset/", stringify!($method)),
            $old,
            $new,
            |state: &mut Old| state.$method(&[5, 6, 7]),
            |state: &mut New| state.$method(&[5, 6, 7], &mut ppvm_conformance_2::analytic_rng()),
        );
    };
}

pub fn register(c: &mut Criterion) {
    let (case_b_old, case_b_new) = branch_pair(4);
    bench_output(
        c,
        "mixture/measure/case_b",
        &case_b_old,
        &case_b_new,
        |state: &mut Old| state.measure(5),
        |state: &mut New| state.measure(5),
    );

    let (mut case_a_old, mut case_a_new) = branch_pair(4);
    case_a_old.h(5);
    case_a_new.h(5);
    assert_same(&case_a_old, &case_a_new);
    bench_output(
        c,
        "mixture/measure/case_a",
        &case_a_old,
        &case_a_new,
        |state: &mut Old| state.measure(5),
        |state: &mut New| state.measure(5),
    );

    let (mut lost_old, mut lost_new) = branch_pair(2);
    lost_old.loss_channel(5, 1.0);
    lost_new.loss_channel(5, 1.0, &mut ppvm_conformance_2::analytic_rng());
    assert_same(&lost_old, &lost_new);
    bench_output(
        c,
        "mixture/measure/lost",
        &lost_old,
        &lost_new,
        |state: &mut Old| state.measure(5),
        |state: &mut New| state.measure(5),
    );

    for qubit in 5..8 {
        case_a_old.h(qubit);
        case_a_new.h(qubit);
    }
    assert_same(&case_a_old, &case_a_new);
    reset_one!(c, &case_a_old, &case_a_new, reset);
    reset_one!(c, &case_a_old, &case_a_new, reset_z);
    reset_one!(c, &case_a_old, &case_a_new, reset_x);
    reset_one!(c, &case_a_old, &case_a_new, reset_y);
    reset_many!(c, &case_a_old, &case_a_new, reset_many);
    reset_many!(c, &case_a_old, &case_a_new, reset_z_many);
    reset_many!(c, &case_a_old, &case_a_new, reset_x_many);
    reset_many!(c, &case_a_old, &case_a_new, reset_y_many);
}
