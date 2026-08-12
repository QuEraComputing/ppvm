// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_conformance_2::tableau::Driver;

use super::*;

macro_rules! gen_checked {
    ($group:expr, $name:expr, $old:expr, $new:expr, $old_op:expr, $new_op:expr) => {{
        let (mut oc, mut nc) = (($old).fork(Some(SEED + 1)), ($new).fork(Some(SEED + 1)));
        let (ro, rn) = (($old_op)(&mut oc), ($new_op)(&mut nc));
        assert_eq!(ro, rn);
        assert_gen_eq(&oc, &nc);
        bench_mut_pair!($group, $name, $old, $new, $old_op, $new_op);
    }};
}

pub fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("tableau-surface/measurement");
    let (old_bd, new_bd) = bare_pair(96);
    let (mut old_br, mut new_br) = bare_pair(96);
    ppvm_traits::traits::Clifford::h(&mut old_br, 0);
    ppvm_traits_2::Clifford::h(&mut new_br, 0);

    assert_eq!(
        ppvm_traits::traits::Measure::measure(&mut old_bd.clone(), 0),
        ppvm_traits_2::Measure::measure(
            &mut new_bd.clone(),
            0,
            &mut ppvm_conformance_2::analytic_rng(),
        )
        .unwrap()
    );
    bench_mut_pair!(
        group,
        "bare/measure-deterministic",
        old_bd,
        new_bd,
        |t: &mut OldBare| ppvm_traits::traits::Measure::measure(t, 0),
        |t: &mut NewBare| {
            ppvm_traits_2::Measure::measure(t, 0, &mut ppvm_conformance_2::analytic_rng())
        }
    );
    assert_eq!(
        ppvm_traits::traits::Measure::measure(&mut old_br.clone(), 0),
        ppvm_traits_2::Measure::measure(
            &mut new_br.clone(),
            0,
            &mut ppvm_conformance_2::analytic_rng(),
        )
        .unwrap()
    );
    bench_mut_pair!(
        group,
        "bare/measure-random",
        old_br,
        new_br,
        |t: &mut OldBare| ppvm_traits::traits::Measure::measure(t, 0),
        |t: &mut NewBare| {
            ppvm_traits_2::Measure::measure(t, 0, &mut ppvm_conformance_2::analytic_rng())
        }
    );

    // NEW-only: old bare `Measure` has no batch method.
    let targets: Vec<usize> = (0..16).collect();
    group.bench_function("bare/measure_many/new-only", |b| {
        b.iter_batched_ref(
            || new_br.clone(),
            |t| {
                std::hint::black_box(ppvm_traits_2::Measure::measure_many(
                    t,
                    &targets,
                    &mut ppvm_conformance_2::analytic_rng(),
                ))
            },
            BatchSize::SmallInput,
        )
    });

    let (old_gd, new_gd) = gen_pair(96);
    let (mut old_gr, mut new_gr) = gen_pair(96);
    for q in 0..16 {
        old_gr.h(q);
        new_gr.h(q);
    }
    gen_checked!(
        group,
        "generalized/measure-deterministic",
        old_gd,
        new_gd,
        |t: &mut OldGen| t.measure(80),
        |t: &mut NewGen| t.measure(80)
    );
    gen_checked!(
        group,
        "generalized/measure-random",
        old_gr,
        new_gr,
        |t: &mut OldGen| t.measure(0),
        |t: &mut NewGen| t.measure(0)
    );
    gen_checked!(
        group,
        "generalized/measure_noisy",
        old_gr,
        new_gr,
        |t: &mut OldGen| t.measure_noisy(0, 0.1),
        |t: &mut NewGen| t.measure_noisy(0, 0.1)
    );
    gen_checked!(
        group,
        "generalized/measure_many",
        old_gr,
        new_gr,
        |t: &mut OldGen| t.measure_many(&targets),
        |t: &mut NewGen| t.measure_many(&targets)
    );
    gen_checked!(
        group,
        "generalized/measure_all",
        old_gr,
        new_gr,
        |t: &mut OldGen| t.measure_all(),
        |t: &mut NewGen| t.measure_all()
    );

    super::measurement_scratch::bench(&mut group, &old_gr, &new_gr, &targets);
    super::measurement_reset::bench(&mut group, &old_br, &new_br, &old_gr, &new_gr, &targets);
    group.finish();
}
