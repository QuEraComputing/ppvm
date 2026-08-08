// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;

macro_rules! single {
    ($group:expr, $old_b:expr, $new_b:expr, $old_g:expr, $new_g:expr, $trait:ident, $op:ident) => {{
        let (mut ob, mut nb) = (($old_b).clone(), ($new_b).clone());
        ppvm_traits::traits::$trait::$op(&mut ob, 0);
        ppvm_traits_2::$trait::$op(&mut nb, 0);
        assert_bare_eq(&ob, &nb);
        let (mut og, mut ng) = (($old_g).clone(), ($new_g).clone());
        ppvm_traits::traits::$trait::$op(&mut og, 0);
        ppvm_traits_2::$trait::$op(&mut ng, 0);
        assert_gen_eq(&og, &ng);
        bench_mut_pair!(
            $group,
            concat!("bare/", stringify!($op)),
            $old_b,
            $new_b,
            |t: &mut OldBare| ppvm_traits::traits::$trait::$op(t, 0),
            |t: &mut NewBare| ppvm_traits_2::$trait::$op(t, 0)
        );
        bench_mut_pair!(
            $group,
            concat!("generalized/", stringify!($op)),
            $old_g,
            $new_g,
            |t: &mut OldGen| ppvm_traits::traits::$trait::$op(t, 0),
            |t: &mut NewGen| ppvm_traits_2::$trait::$op(t, 0)
        );
    }};
}

macro_rules! pair {
    ($group:expr, $old_b:expr, $new_b:expr, $old_g:expr, $new_g:expr, $trait:ident, $op:ident) => {{
        let (mut ob, mut nb) = (($old_b).clone(), ($new_b).clone());
        ppvm_traits::traits::$trait::$op(&mut ob, 1, 70);
        ppvm_traits_2::$trait::$op(&mut nb, 1, 70);
        assert_bare_eq(&ob, &nb);
        let (mut og, mut ng) = (($old_g).clone(), ($new_g).clone());
        ppvm_traits::traits::$trait::$op(&mut og, 1, 70);
        ppvm_traits_2::$trait::$op(&mut ng, 1, 70);
        assert_gen_eq(&og, &ng);
        bench_mut_pair!(
            $group,
            concat!("bare/", stringify!($op)),
            $old_b,
            $new_b,
            |t: &mut OldBare| ppvm_traits::traits::$trait::$op(t, 1, 70),
            |t: &mut NewBare| ppvm_traits_2::$trait::$op(t, 1, 70)
        );
        bench_mut_pair!(
            $group,
            concat!("generalized/", stringify!($op)),
            $old_g,
            $new_g,
            |t: &mut OldGen| ppvm_traits::traits::$trait::$op(t, 1, 70),
            |t: &mut NewGen| ppvm_traits_2::$trait::$op(t, 1, 70)
        );
    }};
}

pub fn bench(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    old_b: &OldBare,
    new_b: &NewBare,
    old_g: &OldGen,
    new_g: &NewGen,
) {
    single!(group, old_b, new_b, old_g, new_g, Clifford, x);
    single!(group, old_b, new_b, old_g, new_g, Clifford, y);
    single!(group, old_b, new_b, old_g, new_g, Clifford, z);
    single!(group, old_b, new_b, old_g, new_g, Clifford, h);
    single!(group, old_b, new_b, old_g, new_g, Clifford, s);
    single!(group, old_b, new_b, old_g, new_g, CliffordExtensions, s_dag);
    single!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordExtensions,
        sqrt_x
    );
    single!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordExtensions,
        sqrt_x_dag
    );
    single!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordExtensions,
        sqrt_y
    );
    single!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordExtensions,
        sqrt_y_dag
    );
    pair!(group, old_b, new_b, old_g, new_g, Clifford, cnot);
    pair!(group, old_b, new_b, old_g, new_g, Clifford, cz);
    pair!(group, old_b, new_b, old_g, new_g, Clifford, cx);
    pair!(group, old_b, new_b, old_g, new_g, Clifford, zcx);
    pair!(group, old_b, new_b, old_g, new_g, Clifford, zcz);
    pair!(group, old_b, new_b, old_g, new_g, CliffordExtensions, cy);
    pair!(group, old_b, new_b, old_g, new_g, CliffordExtensions, zcy);
}
