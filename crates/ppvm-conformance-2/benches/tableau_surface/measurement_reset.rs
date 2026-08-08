// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;

macro_rules! reset_case {
    ($group:expr, $name:expr, $op:ident, $old_b:expr, $new_b:expr,
     $old_g:expr, $new_g:expr, $arg:expr) => {{
        let (mut ob, mut nb) = (($old_b).clone(), ($new_b).clone());
        ppvm_traits::traits::Reset::$op(&mut ob, $arg);
        ppvm_traits_2::Reset::$op(&mut nb, $arg);
        assert_bare_eq(&ob, &nb);
        let (mut og, mut ng) = (($old_g).fork(Some(SEED + 1)), ($new_g).fork(Some(SEED + 1)));
        ppvm_traits::traits::Reset::$op(&mut og, $arg);
        ppvm_traits_2::Reset::$op(&mut ng, $arg);
        assert_gen_eq(&og, &ng);
        bench_mut_pair!(
            $group,
            concat!("bare/", $name),
            $old_b,
            $new_b,
            |t: &mut OldBare| ppvm_traits::traits::Reset::$op(t, $arg),
            |t: &mut NewBare| ppvm_traits_2::Reset::$op(t, $arg)
        );
        bench_mut_pair!(
            $group,
            concat!("generalized/", $name),
            $old_g,
            $new_g,
            |t: &mut OldGen| ppvm_traits::traits::Reset::$op(t, $arg),
            |t: &mut NewGen| ppvm_traits_2::Reset::$op(t, $arg)
        );
    }};
}

pub fn bench(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    old_b: &OldBare,
    new_b: &NewBare,
    old_g: &OldGen,
    new_g: &NewGen,
    qs: &[usize],
) {
    reset_case!(group, "reset", reset, old_b, new_b, old_g, new_g, 0);
    reset_case!(group, "reset_z", reset_z, old_b, new_b, old_g, new_g, 0);
    reset_case!(group, "reset_x", reset_x, old_b, new_b, old_g, new_g, 0);
    reset_case!(group, "reset_y", reset_y, old_b, new_b, old_g, new_g, 0);
    reset_case!(
        group,
        "reset_many",
        reset_many,
        old_b,
        new_b,
        old_g,
        new_g,
        qs
    );
    reset_case!(
        group,
        "reset_z_many",
        reset_z_many,
        old_b,
        new_b,
        old_g,
        new_g,
        qs
    );
    reset_case!(
        group,
        "reset_x_many",
        reset_x_many,
        old_b,
        new_b,
        old_g,
        new_g,
        qs
    );
    reset_case!(
        group,
        "reset_y_many",
        reset_y_many,
        old_b,
        new_b,
        old_g,
        new_g,
        qs
    );
}
