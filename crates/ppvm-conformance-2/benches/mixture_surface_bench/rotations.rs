// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_conformance_2::mixture::{New, Old};
use ppvm_traits::char::Pauli as OldPauli;
use ppvm_traits::traits::{
    RotationOne as OldRotationOne, RotationTwo as OldRotationTwo, TGate as OldTGate,
    U3Gate as OldU3,
};
use ppvm_traits_2::{
    Pauli as NewPauli, RotationOne as NewRotationOne, RotationTwo as NewRotationTwo,
    TGate as NewTGate, U3Gate as NewU3,
};

use super::support::{bench_mut, branch_pair};

macro_rules! rotation_two {
    ($c:expr, $old:expr, $new:expr, $method:ident) => {
        bench_mut(
            $c,
            concat!("mixture/rotation/", stringify!($method)),
            $old,
            $new,
            |state: &mut Old| state.$method(4, 5, 0.37),
            |state: &mut New| state.$method(4, 5, 0.37),
        );
    };
}

macro_rules! rotation_two_many {
    ($c:expr, $old:expr, $new:expr, $method:ident) => {
        bench_mut(
            $c,
            concat!("mixture/rotation/", stringify!($method)),
            $old,
            $new,
            |state: &mut Old| state.$method(&[(2, 3), (4, 5), (6, 7)], 0.19),
            |state: &mut New| state.$method(&[(2, 3), (4, 5), (6, 7)], 0.19),
        );
    };
}

pub fn register(c: &mut Criterion) {
    let (old, new) = branch_pair(4);

    bench_mut(
        c,
        "mixture/rotation/t",
        &old,
        &new,
        |s: &mut Old| s.t(3),
        |s: &mut New| s.t(3),
    );
    bench_mut(
        c,
        "mixture/rotation/t_dag",
        &old,
        &new,
        |s: &mut Old| s.t_dag(3),
        |s: &mut New| s.t_dag(3),
    );
    bench_mut(
        c,
        "mixture/rotation/t_many",
        &old,
        &new,
        |s: &mut Old| s.t_many(&[2, 4, 6]),
        |s: &mut New| s.t_many(&[2, 4, 6]),
    );
    bench_mut(
        c,
        "mixture/rotation/t_dag_many",
        &old,
        &new,
        |s: &mut Old| s.t_dag_many(&[2, 4, 6]),
        |s: &mut New| s.t_dag_many(&[2, 4, 6]),
    );

    bench_mut(
        c,
        "mixture/rotation/rotate_1",
        &old,
        &new,
        |s: &mut Old| s.rotate_1(OldPauli::X, 3, 0.31),
        |s: &mut New| s.rotate_1(NewPauli::X, 3, 0.31),
    );
    for (name, old_op, new_op) in [
        (
            "rx",
            Old::rx as fn(&mut Old, usize, f64),
            New::rx as fn(&mut New, usize, f64),
        ),
        ("ry", Old::ry, New::ry),
        ("rz", Old::rz, New::rz),
    ] {
        bench_mut(
            c,
            &format!("mixture/rotation/{name}"),
            &old,
            &new,
            move |s: &mut Old| old_op(s, 3, 0.31),
            move |s: &mut New| new_op(s, 3, 0.31),
        );
    }
    bench_mut(
        c,
        "mixture/rotation/rx_many",
        &old,
        &new,
        |s: &mut Old| s.rx_many(&[2, 4, 6], 0.23),
        |s: &mut New| s.rx_many(&[2, 4, 6], 0.23),
    );
    bench_mut(
        c,
        "mixture/rotation/ry_many",
        &old,
        &new,
        |s: &mut Old| s.ry_many(&[2, 4, 6], 0.23),
        |s: &mut New| s.ry_many(&[2, 4, 6], 0.23),
    );
    bench_mut(
        c,
        "mixture/rotation/rz_many",
        &old,
        &new,
        |s: &mut Old| s.rz_many(&[2, 4, 6], 0.23),
        |s: &mut New| s.rz_many(&[2, 4, 6], 0.23),
    );

    bench_mut(
        c,
        "mixture/rotation/rotate_2",
        &old,
        &new,
        |s: &mut Old| s.rotate_2([1, 0], [0, 1], 4, 5, 0.37),
        |s: &mut New| s.rotate_2([1, 0], [0, 1], 4, 5, 0.37),
    );
    rotation_two!(c, &old, &new, rxx);
    rotation_two!(c, &old, &new, rxy);
    rotation_two!(c, &old, &new, rxz);
    rotation_two!(c, &old, &new, ryx);
    rotation_two!(c, &old, &new, ryy);
    rotation_two!(c, &old, &new, ryz);
    rotation_two!(c, &old, &new, rzx);
    rotation_two!(c, &old, &new, rzy);
    rotation_two!(c, &old, &new, rzz);
    rotation_two_many!(c, &old, &new, rxx_many);
    rotation_two_many!(c, &old, &new, rxy_many);
    rotation_two_many!(c, &old, &new, rxz_many);
    rotation_two_many!(c, &old, &new, ryx_many);
    rotation_two_many!(c, &old, &new, ryy_many);
    rotation_two_many!(c, &old, &new, ryz_many);
    rotation_two_many!(c, &old, &new, rzx_many);
    rotation_two_many!(c, &old, &new, rzy_many);
    rotation_two_many!(c, &old, &new, rzz_many);

    bench_mut(
        c,
        "mixture/rotation/u3",
        &old,
        &new,
        |s: &mut Old| s.u3(3, 0.2, -0.4, 0.7),
        |s: &mut New| s.u3(3, 0.2, -0.4, 0.7),
    );
}
