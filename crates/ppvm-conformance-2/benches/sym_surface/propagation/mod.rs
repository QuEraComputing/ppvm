// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod blockers;
mod clifford;
mod noise;
mod rotation_one;
mod rotation_two;

use criterion::{BatchSize, BenchmarkGroup, measurement::WallTime};
use ppvm_conformance_2::sym::{
    NewSymKey, NewSymSum, NewTerm, OldSymSum, OldTerm, new_sym_support, new_view, old_sym_support,
    old_view,
};

use super::{assert_real, new_sum, old_sum};

pub(super) fn bench(c: &mut criterion::Criterion) {
    clifford::bench(c);
    rotation_one::bench(c);
    rotation_two::bench(c);
    noise::bench(c);
    blockers::bench(c);
}

fn fixture() -> (OldSymSum, NewSymSum) {
    let mut old = old_sum(8);
    old += ("XYZIXYZI", OldTerm::from(1.0));
    old += ("ZYXIZYXI", OldTerm::var(0).sin() + OldTerm::var(1).cos());
    let mut new = new_sum(8);
    new += (NewSymKey::from("XYZIXYZI"), NewTerm::from(1.0));
    new += (
        NewSymKey::from("ZYXIZYXI"),
        NewTerm::var(0).sin() + NewTerm::var(1).cos(),
    );
    assert_eq!(old.capacity(), new.capacity());
    (old, new)
}

fn paired<O, N>(group: &mut BenchmarkGroup<'_, WallTime>, name: &str, old_op: O, new_op: N)
where
    O: Fn(&mut OldSymSum) + Copy,
    N: Fn(&mut NewSymSum) + Copy,
{
    let (old, new) = fixture();
    let mut old_expected = old.clone();
    let mut new_expected = new.clone();
    old_op(&mut old_expected);
    new_op(&mut new_expected);
    assert_sums(&old_expected, &new_expected);

    group.bench_function(format!("new/{name}"), |b| {
        b.iter_batched(
            || new.clone(),
            |mut sum| {
                new_op(&mut sum);
                sum
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function(format!("old/{name}"), |b| {
        b.iter_batched(
            || old.clone(),
            |mut sum| {
                old_op(&mut sum);
                sum
            },
            BatchSize::SmallInput,
        )
    });
}

fn paired_args<OA, NA, O, N>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    name: &str,
    old_arg: OA,
    new_arg: NA,
    old_op: O,
    new_op: N,
) where
    OA: Clone,
    NA: Clone,
    O: Fn(&mut OldSymSum, OA) + Copy,
    N: Fn(&mut NewSymSum, NA) + Copy,
{
    let (old, new) = fixture();
    let mut old_expected = old.clone();
    let mut new_expected = new.clone();
    old_op(&mut old_expected, old_arg.clone());
    new_op(&mut new_expected, new_arg.clone());
    assert_sums(&old_expected, &new_expected);

    group.bench_function(format!("new/{name}"), |b| {
        b.iter_batched(
            || (new.clone(), new_arg.clone()),
            |(mut sum, arg)| {
                new_op(&mut sum, arg);
                sum
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function(format!("old/{name}"), |b| {
        b.iter_batched(
            || (old.clone(), old_arg.clone()),
            |(mut sum, arg)| {
                old_op(&mut sum, arg);
                sum
            },
            BatchSize::SmallInput,
        )
    });
}

fn assert_sums(old: &OldSymSum, new: &NewSymSum) {
    let old = old_sym_support(old);
    let new = new_sym_support(new);
    assert_eq!(
        old.iter().map(|(key, _)| key).collect::<Vec<_>>(),
        new.iter().map(|(key, _)| key).collect::<Vec<_>>()
    );
    for ((_, old), (_, new)) in old.iter().zip(&new) {
        assert_eq!(old_view(old).monomials, new_view(new).monomials);
        assert_real(
            old.eval(&[0.3, -0.7]).unwrap(),
            new.eval(&[0.3, -0.7]).unwrap(),
        );
    }
}
