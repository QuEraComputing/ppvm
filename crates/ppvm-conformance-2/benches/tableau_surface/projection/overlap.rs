// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use fxhash::FxHashMap;

use super::*;

type RawNew = ppvm_tableau_2::GeneralizedTableau<[usize; 2], u128>;

pub fn bench(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    old: &OldGen,
    new: &NewGen,
) {
    let (old_entries, new_entries) = entries(old, new);
    let b = decomposition(old, new, 7);
    assert_eq!(b.stabilizer, 0);
    let old_b = OldGen::compute_overlap_case_b(&old_entries, b.phase, b.destabilizer);
    let new_b = RawNew::compute_overlap_case_b(&new_entries, b.phase, b.destabilizer);
    assert!((old_b - new_b).abs() <= 1e-12);
    group.bench_function("compute_overlap_case_b/old", |bench| {
        bench.iter(|| {
            std::hint::black_box(OldGen::compute_overlap_case_b(
                &old_entries,
                b.phase,
                b.destabilizer,
            ))
        })
    });
    group.bench_function("compute_overlap_case_b/new", |bench| {
        bench.iter(|| {
            std::hint::black_box(RawNew::compute_overlap_case_b(
                &new_entries,
                b.phase,
                b.destabilizer,
            ))
        })
    });

    let a = decomposition(old, new, 0);
    assert_ne!(a.stabilizer, 0);
    let old_map: FxHashMap<u128, Complex64> = old_entries.iter().map(|&(c, i)| (i, c)).collect();
    let new_map: FxHashMap<u128, Complex64> = new_entries.iter().map(|&(c, i)| (i, c)).collect();
    assert_eq!(old_map, new_map);
    let (old_mask, new_mask) = (
        old.odd_phase_destabilizer_mask(),
        new.odd_phase_destabilizer_mask(),
    );
    assert_eq!(old_mask, new_mask);
    let old_a =
        OldGen::compute_overlap_case_a(&old_map, a.phase, a.destabilizer, a.stabilizer, old_mask);
    let new_a =
        RawNew::compute_overlap_case_a(&new_map, a.phase, a.destabilizer, a.stabilizer, new_mask);
    assert!((old_a - new_a).abs() <= 1e-12);
    group.bench_function("compute_overlap_case_a/old", |bench| {
        bench.iter(|| {
            std::hint::black_box(OldGen::compute_overlap_case_a(
                &old_map,
                a.phase,
                a.destabilizer,
                a.stabilizer,
                old_mask,
            ))
        })
    });
    group.bench_function("compute_overlap_case_a/new", |bench| {
        bench.iter(|| {
            std::hint::black_box(RawNew::compute_overlap_case_a(
                &new_map,
                a.phase,
                a.destabilizer,
                a.stabilizer,
                new_mask,
            ))
        })
    });
}
