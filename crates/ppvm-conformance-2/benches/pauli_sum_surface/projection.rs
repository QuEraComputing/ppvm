// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hint::black_box;

use criterion::{BatchSize, Criterion};
use ppvm_traits::traits::Projection as OldProjection;
use ppvm_traits_2::Projection as NewProjection;

use super::{NewKey, NewSum, assert_pair, build_new, build_old};

pub fn bench(c: &mut Criterion) {
    comparable(c, false);
    comparable(c, true);
    corrected_mixed(c, false);
    corrected_mixed(c, true);
}

fn comparable(c: &mut Criterion, one: bool) {
    let data = vec![("IIIIIIII".to_string(), 1.0), ("IIIZIIII".to_string(), 1.0)];
    let (old, new) = (build_old(&data), build_new(&data));
    let (mut op, mut np) = (old.clone(), new.clone());
    if one {
        op.p1(3);
        np.p1(3);
    } else {
        op.p0(3);
        np.p0(3);
    }
    assert_pair(&op, &np);
    let name = if one { "p1" } else { "p0" };
    let mut group = c.benchmark_group(format!("pauli_sum_surface/projection/{name}_iz_unit"));
    group.bench_function("old", |b| {
        b.iter_batched_ref(
            || old.clone(),
            |s| {
                if one { s.p1(3) } else { s.p0(3) }
            },
            BatchSize::LargeInput,
        )
    });
    group.bench_function("new", |b| {
        b.iter_batched_ref(
            || new.clone(),
            |s| {
                if one { s.p1(3) } else { s.p0(3) }
            },
            BatchSize::LargeInput,
        )
    });
    group.finish();
}

fn corrected_mixed(c: &mut Criterion, one: bool) {
    let data = vec![
        ("IIIIIIII".to_string(), 2.0),
        ("IIIXIIII".to_string(), -0.75),
        ("IIIYIIII".to_string(), 1.25),
        ("IIIZIIII".to_string(), -0.5),
    ];
    let seed: NewSum = build_new(&data);
    let mut probe = seed.clone();
    if one {
        probe.p1(3);
    } else {
        probe.p0(3);
    }
    assert_eq!(probe.get(&NewKey::from("IIIXIIII")), Some(0.0));
    assert_eq!(probe.get(&NewKey::from("IIIYIIII")), Some(0.0));
    let name = if one { "p1" } else { "p0" };
    let mut group = c.benchmark_group(format!(
        "pauli_sum_surface/projection/{name}_corrected_mixed_new_only"
    ));
    group.bench_function("new", |b| {
        b.iter_batched_ref(
            || seed.clone(),
            |s| {
                if one {
                    s.p1(black_box(3))
                } else {
                    s.p0(black_box(3))
                }
            },
            BatchSize::LargeInput,
        )
    });
    group.finish();
}
