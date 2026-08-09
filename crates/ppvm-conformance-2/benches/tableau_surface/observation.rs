// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;

use super::*;

macro_rules! read_pair {
    ($group:expr, $name:expr, $old:expr, $new:expr, $old_op:expr, $new_op:expr) => {{
        $group.bench_function(concat!($name, "/old"), |b| {
            b.iter(|| std::hint::black_box(($old_op)(&$old)))
        });
        $group.bench_function(concat!($name, "/new"), |b| {
            b.iter(|| std::hint::black_box(($new_op)(&$new)))
        });
    }};
}

pub fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("tableau-surface/observation");
    let (old_b, new_b) = prepared_bare(96);
    let (old_g, new_g) = prepared_gen(96);

    assert_eq!(old_g.n_qubits(), new_g.n_qubits());
    assert_eq!(
        old_g.current_measurement_record(),
        new_g.current_measurement_record()
    );
    read_pair!(
        group,
        "generalized/n_qubits",
        old_g,
        new_g,
        |t: &OldGen| t.n_qubits(),
        |t: &NewGen| t.n_qubits()
    );
    read_pair!(
        group,
        "generalized/current_measurement_record",
        old_g,
        new_g,
        |t: &OldGen| t.current_measurement_record().len(),
        |t: &NewGen| t.current_measurement_record().len()
    );

    assert_eq!(
        old_g.compute_decomposition(65, ppvm_traits::Pauli::Z),
        new_g.compute_decomposition(65, ppvm_traits_2::Pauli::Z)
    );
    read_pair!(
        group,
        "generalized/compute_decomposition",
        old_g,
        new_g,
        |t: &OldGen| t.compute_decomposition(65, ppvm_traits::Pauli::Z),
        |t: &NewGen| t.compute_decomposition(65, ppvm_traits_2::Pauli::Z)
    );
    assert_eq!(
        old_g.odd_phase_destabilizer_mask(),
        new_g.odd_phase_destabilizer_mask()
    );
    read_pair!(
        group,
        "generalized/odd_phase_destabilizer_mask",
        old_g,
        new_g,
        |t: &OldGen| t.odd_phase_destabilizer_mask(),
        |t: &NewGen| t.odd_phase_destabilizer_mask()
    );

    let word = format!("ZZ{}", "I".repeat(94));
    let old_word: ppvm_pauli_word::word::PauliWord<[usize; 2]> = word.as_str().into();
    let new_word: ppvm_pauli_word_2::PauliWord<[usize; 2]> = word.as_str().into();
    let oe = old_g.expectation(&old_word);
    let ne = new_g.expectation(&new_word);
    assert!((oe - ne).abs() <= 1e-10);
    read_pair!(
        group,
        "generalized/expectation",
        old_g,
        new_g,
        |t: &OldGen| t.expectation(&old_word),
        |t: &NewGen| t.expectation(&new_word)
    );
    let (oz, nz) = (old_g.z_expectation(65), new_g.z_expectation(65));
    assert!((oz - nz).abs() <= 1e-10);
    read_pair!(
        group,
        "generalized/z_expectation",
        old_g,
        new_g,
        |t: &OldGen| t.z_expectation(65),
        |t: &NewGen| t.z_expectation(65)
    );

    let old_pattern =
        ppvm_pauli_word::pattern::PauliPattern::parse("Z?{4}").expect("valid old pattern");
    let new_pattern = ppvm_pauli_sum_2::PauliPattern::parse("Z?{4}").expect("valid new pattern");
    let (ot, nt) = (old_g.trace(&old_pattern), new_g.trace(&new_pattern));
    assert!((ot - nt).abs() <= 1e-10);
    read_pair!(
        group,
        "generalized/trace-pattern",
        old_g,
        new_g,
        |t: &OldGen| t.trace(&old_pattern),
        |t: &NewGen| t.trace(&new_pattern)
    );

    record_mutations(&mut group, &old_g, &new_g);

    // Iteration shapes differ (old returns slices, new returns row iterators), so
    // compare equal full traversals but keep the labels explicit rather than
    // pretending an accessor-only ratio is equivalent.
    assert_bare_eq(&old_b, &new_b);
    group.bench_function("bare/stabilizers-traverse/old-api", |b| {
        b.iter(|| std::hint::black_box(old_b.stabilizers().iter().count()))
    });
    group.bench_function("bare/stabilizers-traverse/new-api", |b| {
        b.iter(|| std::hint::black_box(new_b.stabilizer_rows().count()))
    });
    group.bench_function("bare/destabilizers-traverse/old-api", |b| {
        b.iter(|| std::hint::black_box(old_b.destabilizers().iter().count()))
    });
    group.bench_function("bare/destabilizers-traverse/new-api", |b| {
        b.iter(|| std::hint::black_box(new_b.destabilizer_rows().count()))
    });
    group.bench_function("bare/n_qubits/new-only", |b| {
        b.iter(|| std::hint::black_box(new_b.n_qubits()))
    });
    group.bench_function("bare/row_site/new-only", |b| {
        b.iter(|| std::hint::black_box(new_b.row_site(0, 65)))
    });
    group.finish();
}

fn record_mutations(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    old: &OldGen,
    new: &NewGen,
) {
    bench_mut_pair!(
        group,
        "generalized/append_measurement_record",
        old,
        new,
        |t: &mut OldGen| t.append_measurement_record(Some(true)),
        |t: &mut NewGen| t.append_measurement_record(Some(true))
    );
    let mut old_record = old.clone();
    let mut new_record = new.clone();
    old_record.append_measurement_record(Some(true));
    new_record.append_measurement_record(Some(true));
    assert_gen_eq(&old_record, &new_record);
    bench_mut_pair!(
        group,
        "generalized/overwrite_last_measurement_record",
        old_record,
        new_record,
        |t: &mut OldGen| t.overwrite_last_measurement_record(Some(false)),
        |t: &mut NewGen| t.overwrite_last_measurement_record(Some(false))
    );
}
