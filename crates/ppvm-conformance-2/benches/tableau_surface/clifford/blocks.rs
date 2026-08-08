// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;

pub fn bench(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    old_b: &OldBare,
    new_b: &NewBare,
    old_g: &OldGen,
    new_g: &NewGen,
) {
    let (mut ob, mut nb) = (old_b.clone(), new_b.clone());
    ob.cz_block_pairs(0, 17, 17);
    nb.cz_block_pairs(0, 17, 17);
    assert_bare_eq(&ob, &nb);
    let (mut ob, mut nb) = (old_b.clone(), new_b.clone());
    ob.cz_block_pairs_cross_word(0, 0, 1, 0, 17);
    nb.cz_block_pairs_cross_word(0, 0, 1, 0, 17);
    assert_bare_eq(&ob, &nb);
    let (mut og, mut ng) = (old_g.clone(), new_g.clone());
    og.cz_block(0, 65, 17);
    ng.cz_block(0, 65, 17);
    assert_gen_eq(&og, &ng);

    bench_mut_pair!(
        group,
        "bare/cz_block_pairs",
        old_b,
        new_b,
        |t: &mut OldBare| t.cz_block_pairs(0, 17, 17),
        |t: &mut NewBare| t.cz_block_pairs(0, 17, 17)
    );
    bench_mut_pair!(
        group,
        "bare/cz_block_pairs_cross_word",
        old_b,
        new_b,
        |t: &mut OldBare| t.cz_block_pairs_cross_word(0, 0, 1, 0, 17),
        |t: &mut NewBare| t.cz_block_pairs_cross_word(0, 0, 1, 0, 17)
    );
    bench_mut_pair!(
        group,
        "generalized/cz_block_pairs",
        old_g,
        new_g,
        |t: &mut OldGen| t.cz_block_pairs(0, 17, 17),
        |t: &mut NewGen| t.cz_block_pairs(0, 17, 17)
    );
    bench_mut_pair!(
        group,
        "generalized/cz_block_pairs_cross_word",
        old_g,
        new_g,
        |t: &mut OldGen| t.cz_block_pairs_cross_word(0, 0, 1, 0, 17),
        |t: &mut NewGen| t.cz_block_pairs_cross_word(0, 0, 1, 0, 17)
    );
    bench_mut_pair!(
        group,
        "generalized/cz_block",
        old_g,
        new_g,
        |t: &mut OldGen| t.cz_block(0, 65, 17),
        |t: &mut NewGen| t.cz_block(0, 65, 17)
    );
}

pub fn width_sweep(group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    for n in WIDTHS {
        let (old, new) = prepared_bare(n);
        bench_mut_pair!(
            group,
            format!("width-{n}/bare/h"),
            old,
            new,
            |t: &mut OldBare| ppvm_traits::traits::Clifford::h(t, 0),
            |t: &mut NewBare| ppvm_traits_2::Clifford::h(t, 0)
        );
        bench_mut_pair!(
            group,
            format!("width-{n}/bare/cnot-edge"),
            old,
            new,
            |t: &mut OldBare| ppvm_traits::traits::Clifford::cnot(t, 0, n - 1),
            |t: &mut NewBare| ppvm_traits_2::Clifford::cnot(t, 0, n - 1)
        );
        let (old, new) = prepared_gen(n);
        bench_mut_pair!(
            group,
            format!("width-{n}/generalized/h"),
            old,
            new,
            |t: &mut OldGen| ppvm_traits::traits::Clifford::h(t, 0),
            |t: &mut NewGen| ppvm_traits_2::Clifford::h(t, 0)
        );
        bench_mut_pair!(
            group,
            format!("width-{n}/generalized/cnot-edge"),
            old,
            new,
            |t: &mut OldGen| ppvm_traits::traits::Clifford::cnot(t, 0, n - 1),
            |t: &mut NewGen| ppvm_traits_2::Clifford::cnot(t, 0, n - 1)
        );
    }
}
