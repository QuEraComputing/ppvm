// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_conformance_2::mixture::{New, Old};
use ppvm_traits::traits::{
    Clifford as OldClifford, CliffordBatch as OldCliffordBatch,
    CliffordExtensions as OldExtensions, CliffordExtensionsBatch as OldExtensionsBatch,
};
use ppvm_traits_2::{
    Clifford as NewClifford, CliffordBatch as NewCliffordBatch,
    CliffordExtensions as NewExtensions, CliffordExtensionsBatch as NewExtensionsBatch,
};

use super::support::{assert_same, bench_mut, branch_pair};

macro_rules! unary {
    ($c:expr, $old:expr, $new:expr, $method:ident) => {
        bench_mut(
            $c,
            concat!("mixture/gate/", stringify!($method)),
            $old,
            $new,
            |state: &mut Old| state.$method(3),
            |state: &mut New| state.$method(3),
        );
    };
}

macro_rules! binary {
    ($c:expr, $old:expr, $new:expr, $method:ident) => {
        bench_mut(
            $c,
            concat!("mixture/gate/", stringify!($method)),
            $old,
            $new,
            |state: &mut Old| state.$method(2, 5),
            |state: &mut New| state.$method(2, 5),
        );
    };
}

macro_rules! unary_many {
    ($c:expr, $old:expr, $new:expr, $method:ident) => {
        bench_mut(
            $c,
            concat!("mixture/gate/", stringify!($method)),
            $old,
            $new,
            |state: &mut Old| state.$method(&[2, 4, 6, 8]),
            |state: &mut New| state.$method(&[2, 4, 6, 8]),
        );
    };
}

macro_rules! binary_many {
    ($c:expr, $old:expr, $new:expr, $method:ident) => {
        bench_mut(
            $c,
            concat!("mixture/gate/", stringify!($method)),
            $old,
            $new,
            |state: &mut Old| state.$method(&[(2, 3), (4, 5), (6, 7)]),
            |state: &mut New| state.$method(&[(2, 3), (4, 5), (6, 7)]),
        );
    };
}

pub fn register(c: &mut Criterion) {
    let (mut old, mut new) = branch_pair(4);
    for qubit in 2..10 {
        old.h(qubit);
        new.h(qubit);
        if qubit % 2 == 0 {
            old.s(qubit);
            new.s(qubit);
        }
    }
    assert_same(&old, &new);

    unary!(c, &old, &new, x);
    unary!(c, &old, &new, y);
    unary!(c, &old, &new, z);
    unary!(c, &old, &new, h);
    unary!(c, &old, &new, s);
    binary!(c, &old, &new, cnot);
    binary!(c, &old, &new, cz);
    binary!(c, &old, &new, cx);
    binary!(c, &old, &new, zcx);
    binary!(c, &old, &new, zcz);

    unary!(c, &old, &new, s_dag);
    unary!(c, &old, &new, sqrt_x);
    unary!(c, &old, &new, sqrt_x_dag);
    unary!(c, &old, &new, sqrt_y);
    unary!(c, &old, &new, sqrt_y_dag);
    binary!(c, &old, &new, cy);
    binary!(c, &old, &new, zcy);

    unary_many!(c, &old, &new, x_many);
    unary_many!(c, &old, &new, y_many);
    unary_many!(c, &old, &new, z_many);
    unary_many!(c, &old, &new, h_many);
    unary_many!(c, &old, &new, s_many);
    binary_many!(c, &old, &new, cnot_many);
    binary_many!(c, &old, &new, cz_many);
    unary_many!(c, &old, &new, s_dag_many);
    unary_many!(c, &old, &new, sqrt_x_many);
    unary_many!(c, &old, &new, sqrt_x_dag_many);
    unary_many!(c, &old, &new, sqrt_y_many);
    unary_many!(c, &old, &new, sqrt_y_dag_many);
    binary_many!(c, &old, &new, cy_many);
}
