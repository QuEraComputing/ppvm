// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;

use super::*;

macro_rules! checked {
    ($group:expr, $name:expr, $old:expr, $new:expr, $old_op:expr, $new_op:expr) => {{
        let (mut oc, mut nc) = (($old).clone(), ($new).clone());
        ($old_op)(&mut oc);
        ($new_op)(&mut nc);
        assert_gen_eq(&oc, &nc);
        bench_mut_pair!($group, $name, $old, $new, $old_op, $new_op);
    }};
}

macro_rules! rot2 {
    ($group:expr, $old:expr, $new:expr, $op:ident) => {
        checked!(
            $group,
            concat!("two/", stringify!($op)),
            $old,
            $new,
            |t: &mut OldGen| ppvm_traits::traits::RotationTwo::$op(t, 1, 6, 0.37),
            |t: &mut NewGen| ppvm_traits_2::RotationTwo::<num::complex::Complex64, f64>::$op(
                t, 1, 6, 0.37
            )
        );
    };
}

macro_rules! rot2_many {
    ($group:expr, $old:expr, $new:expr, $op:ident, $pairs:expr) => {
        checked!(
            $group,
            concat!("two-batch/", stringify!($op)),
            $old,
            $new,
            |t: &mut OldGen| ppvm_traits::traits::RotationTwo::$op(t, $pairs, 0.37),
            |t: &mut NewGen| ppvm_traits_2::RotationTwo::<num::complex::Complex64, f64>::$op(
                t, $pairs, 0.37
            )
        );
    };
}

pub fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("tableau-surface/rotation");
    let (old, new) = prepared_gen(16);
    let targets = [0usize, 2, 4, 6];
    let pairs = [(0usize, 1usize), (2, 3), (4, 5), (6, 7)];

    checked!(
        group,
        "t",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::TGate::t(t, 0),
        |t: &mut NewGen| ppvm_traits_2::TGate::t(t, 0)
    );
    checked!(
        group,
        "t_dag",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::TGate::t_dag(t, 0),
        |t: &mut NewGen| ppvm_traits_2::TGate::t_dag(t, 0)
    );
    checked!(
        group,
        "t_many",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::TGate::t_many(t, &targets),
        |t: &mut NewGen| ppvm_traits_2::TGate::t_many(t, &targets)
    );
    checked!(
        group,
        "t_dag_many",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::TGate::t_dag_many(t, &targets),
        |t: &mut NewGen| ppvm_traits_2::TGate::t_dag_many(t, &targets)
    );

    checked!(
        group,
        "one/rotate_1_x",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::RotationOne::rotate_1(
            t,
            ppvm_traits::Pauli::X,
            0,
            0.37
        ),
        |t: &mut NewGen| ppvm_traits_2::RotationOne::<num::complex::Complex64, f64>::rotate_1(
            t,
            ppvm_traits_2::Pauli::X,
            0,
            0.37
        )
    );
    checked!(
        group,
        "one/rx",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::RotationOne::rx(t, 0, 0.37),
        |t: &mut NewGen| ppvm_traits_2::RotationOne::<num::complex::Complex64, f64>::rx(t, 0, 0.37)
    );
    checked!(
        group,
        "one/ry",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::RotationOne::ry(t, 0, 0.37),
        |t: &mut NewGen| ppvm_traits_2::RotationOne::<num::complex::Complex64, f64>::ry(t, 0, 0.37)
    );
    checked!(
        group,
        "one/rz",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::RotationOne::rz(t, 0, 0.37),
        |t: &mut NewGen| ppvm_traits_2::RotationOne::<num::complex::Complex64, f64>::rz(t, 0, 0.37)
    );
    checked!(
        group,
        "one/rx_many",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::RotationOne::rx_many(t, &targets, 0.37),
        |t: &mut NewGen| ppvm_traits_2::RotationOne::<num::complex::Complex64, f64>::rx_many(
            t, &targets, 0.37
        )
    );
    checked!(
        group,
        "one/ry_many",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::RotationOne::ry_many(t, &targets, 0.37),
        |t: &mut NewGen| ppvm_traits_2::RotationOne::<num::complex::Complex64, f64>::ry_many(
            t, &targets, 0.37
        )
    );
    checked!(
        group,
        "one/rz_many",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::RotationOne::rz_many(t, &targets, 0.37),
        |t: &mut NewGen| ppvm_traits_2::RotationOne::<num::complex::Complex64, f64>::rz_many(
            t, &targets, 0.37
        )
    );

    checked!(
        group,
        "one/r_xy",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::RotXY::r(t, 0, 0.2, 0.37),
        |t: &mut NewGen| ppvm_traits_2::RotXY::<num::complex::Complex64, f64>::r(t, 0, 0.2, 0.37)
    );
    checked!(
        group,
        "one/u3",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::U3Gate::u3(t, 0, 0.2, 0.3, 0.4),
        |t: &mut NewGen| ppvm_traits_2::U3Gate::<num::complex::Complex64, f64>::u3(
            t, 0, 0.2, 0.3, 0.4
        )
    );

    checked!(
        group,
        "two/rotate_2_xz",
        old,
        new,
        |t: &mut OldGen| ppvm_traits::traits::RotationTwo::rotate_2(t, [1, 0], [0, 1], 1, 6, 0.37),
        |t: &mut NewGen| ppvm_traits_2::RotationTwo::<num::complex::Complex64, f64>::rotate_2(
            t,
            [1, 0],
            [0, 1],
            1,
            6,
            0.37
        )
    );
    rot2!(group, old, new, rxx);
    rot2!(group, old, new, rxy);
    rot2!(group, old, new, rxz);
    rot2!(group, old, new, ryx);
    rot2!(group, old, new, ryy);
    rot2!(group, old, new, ryz);
    rot2!(group, old, new, rzx);
    rot2!(group, old, new, rzy);
    rot2!(group, old, new, rzz);
    rot2_many!(group, old, new, rxx_many, &pairs);
    rot2_many!(group, old, new, rxy_many, &pairs);
    rot2_many!(group, old, new, rxz_many, &pairs);
    rot2_many!(group, old, new, ryx_many, &pairs);
    rot2_many!(group, old, new, ryy_many, &pairs);
    rot2_many!(group, old, new, ryz_many, &pairs);
    rot2_many!(group, old, new, rzx_many, &pairs);
    rot2_many!(group, old, new, rzy_many, &pairs);
    rot2_many!(group, old, new, rzz_many, &pairs);
    group.finish();
}
