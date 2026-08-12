// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;

macro_rules! batch {
    ($group:expr, $old_b:expr, $new_b:expr, $old_g:expr, $new_g:expr,
     $trait:ident, $op:ident, $args:expr) => {{
        let (mut ob, mut nb) = (($old_b).clone(), ($new_b).clone());
        ppvm_traits::traits::$trait::$op(&mut ob, $args);
        ppvm_traits_2::$trait::$op(&mut nb, $args);
        assert_bare_eq(&ob, &nb);
        let (mut og, mut ng) = (($old_g).clone(), ($new_g).clone());
        ppvm_traits::traits::$trait::$op(&mut og, $args);
        ppvm_traits_2::$trait::$op(&mut ng, $args);
        assert_gen_eq(&og, &ng);
        bench_mut_pair!(
            $group,
            concat!("bare/", stringify!($op)),
            $old_b,
            $new_b,
            |t: &mut OldBare| ppvm_traits::traits::$trait::$op(t, $args),
            |t: &mut NewBare| ppvm_traits_2::$trait::$op(t, $args)
        );
        bench_mut_pair!(
            $group,
            concat!("generalized/", stringify!($op)),
            $old_g,
            $new_g,
            |t: &mut OldGen| ppvm_traits::traits::$trait::$op(t, $args),
            |t: &mut NewGen| ppvm_traits_2::$trait::$op(t, $args)
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
    let qs: Vec<usize> = (0..96).step_by(3).collect();
    let pairs: Vec<(usize, usize)> = (0..24).map(|q| (q, q + 65)).collect();
    batch!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordBatch,
        x_many,
        &qs
    );
    batch!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordBatch,
        y_many,
        &qs
    );
    batch!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordBatch,
        z_many,
        &qs
    );
    batch!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordBatch,
        h_many,
        &qs
    );
    batch!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordBatch,
        s_many,
        &qs
    );
    batch!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordExtensionsBatch,
        s_dag_many,
        &qs
    );
    batch!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordExtensionsBatch,
        sqrt_x_many,
        &qs
    );
    batch!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordExtensionsBatch,
        sqrt_x_dag_many,
        &qs
    );
    batch!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordExtensionsBatch,
        sqrt_y_many,
        &qs
    );
    batch!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordExtensionsBatch,
        sqrt_y_dag_many,
        &qs
    );
    batch!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordBatch,
        cnot_many,
        &pairs
    );
    batch!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordBatch,
        cz_many,
        &pairs
    );
    batch!(
        group,
        old_b,
        new_b,
        old_g,
        new_g,
        CliffordExtensionsBatch,
        cy_many,
        &pairs
    );
}
