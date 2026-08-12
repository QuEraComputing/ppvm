// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;

type OldScratch = ppvm_tableau::measure::MeasureScratch<u128, f64>;
type NewScratch = ppvm_tableau_2::MeasureScratch<u128>;
type OldCaseB = (OldGen, Entries);
type NewCaseB = (NewGen, Entries);

fn case_a_inputs(
    old: &OldGen,
    new: &NewGen,
) -> ((OldGen, OldScratch), (NewGen, NewScratch), Decomposition) {
    let d = decomposition(old, new, 0);
    assert_ne!(d.stabilizer, 0, "qubit 0 must exercise case a");
    let (old_entries, new_entries) = entries(old, new);
    let mut old_tab = old.clone();
    let mut new_tab = new.clone();
    old_tab.coefficients = Vec::new();
    new_tab.coefficients = ppvm_tableau_2::Amplitudes::new();

    let mut old_scratch = OldScratch::new();
    let mut new_scratch = NewScratch::new();
    old_scratch.coeff_map.reserve(old_entries.len());
    new_scratch.coeff_map.reserve(new_entries.len());
    old_scratch
        .coeff_map
        .extend(old_entries.into_iter().map(|(c, i)| (i, c)));
    new_scratch
        .coeff_map
        .extend(new_entries.into_iter().map(|(c, i)| (i, c)));
    ((old_tab, old_scratch), (new_tab, new_scratch), d)
}

fn case_b_inputs(old: &OldGen, new: &NewGen) -> (OldCaseB, NewCaseB, Decomposition) {
    let d = decomposition(old, new, 7);
    assert_eq!(d.stabilizer, 0, "qubit 7 must exercise case b");
    let (old_entries, new_entries) = entries(old, new);
    let mut old_tab = old.clone();
    let mut new_tab = new.clone();
    old_tab.coefficients = Vec::new();
    new_tab.coefficients = ppvm_tableau_2::Amplitudes::new();
    ((old_tab, old_entries), (new_tab, new_entries), d)
}

pub fn bench(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    old: &OldGen,
    new: &NewGen,
) {
    let (old_a, new_a, a) = case_a_inputs(old, new);
    let (mut old_check, mut new_check) = (old_a.clone(), new_a.clone());
    old_check.0.project_case_a(
        false,
        &mut old_check.1,
        a.phase,
        a.stabilizer,
        a.destabilizer,
        0,
    );
    new_check.0.project_case_a(
        false,
        &mut new_check.1,
        a.phase,
        a.stabilizer,
        a.destabilizer,
        0,
    );
    assert_gen_eq(&old_check.0, &new_check.0);
    bench_mut_pair!(
        group,
        "project_case_a",
        old_a,
        new_a,
        |(t, s): &mut (OldGen, OldScratch)| t.project_case_a(
            false,
            s,
            a.phase,
            a.stabilizer,
            a.destabilizer,
            0
        ),
        |(t, s): &mut (NewGen, NewScratch)| t.project_case_a(
            false,
            s,
            a.phase,
            a.stabilizer,
            a.destabilizer,
            0
        )
    );

    let (old_b, new_b, b) = case_b_inputs(old, new);
    let (mut old_check, mut new_check) = (old_b.clone(), new_b.clone());
    old_check
        .0
        .project_case_b(&old_check.1, false, b.phase, b.destabilizer);
    new_check
        .0
        .project_case_b(&new_check.1, false, b.phase, b.destabilizer);
    assert_gen_eq(&old_check.0, &new_check.0);
    bench_mut_pair!(
        group,
        "project_case_b",
        old_b,
        new_b,
        |(t, entries): &mut OldCaseB| { t.project_case_b(entries, false, b.phase, b.destabilizer) },
        |(t, entries): &mut NewCaseB| { t.project_case_b(entries, false, b.phase, b.destabilizer) }
    );
}
