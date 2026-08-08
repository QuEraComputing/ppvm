// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_traits::traits::RotationTwo as OldRotationTwo;
use ppvm_traits_2::RotationTwo as NewRotationTwo;

use super::{NewSum, OldSum, bench_mut};

type RotOld = fn(&mut OldSum, usize, usize, f64);
type RotNew = fn(&mut NewSum, usize, usize, f64);
type ManyOld = fn(&mut OldSum, &[(usize, usize)], f64);
type ManyNew = fn(&mut NewSum, &[(usize, usize)], f64);

const THETA: f64 = 0.29;
const PAIRS: &[(usize, usize)] = &[(0, 1), (2, 3), (4, 5), (6, 7)];

pub fn bench(c: &mut Criterion) {
    bench_mut(
        c,
        "rotation_two/rotate_2_generic_xz",
        |s| s.rotate_2([1, 0], [0, 1], 2, 5, THETA),
        |s| s.rotate_2([1, 0], [0, 1], 2, 5, THETA),
    );
    named(c);
    batch(c);
}

fn named(c: &mut Criterion) {
    let gates: [(&str, RotOld, RotNew); 9] = [
        (
            "rxx",
            |s, a, b, t| s.rxx(a, b, t),
            |s, a, b, t| s.rxx(a, b, t),
        ),
        (
            "rxy",
            |s, a, b, t| s.rxy(a, b, t),
            |s, a, b, t| s.rxy(a, b, t),
        ),
        (
            "rxz",
            |s, a, b, t| s.rxz(a, b, t),
            |s, a, b, t| s.rxz(a, b, t),
        ),
        (
            "ryx",
            |s, a, b, t| s.ryx(a, b, t),
            |s, a, b, t| s.ryx(a, b, t),
        ),
        (
            "ryy",
            |s, a, b, t| s.ryy(a, b, t),
            |s, a, b, t| s.ryy(a, b, t),
        ),
        (
            "ryz",
            |s, a, b, t| s.ryz(a, b, t),
            |s, a, b, t| s.ryz(a, b, t),
        ),
        (
            "rzx",
            |s, a, b, t| s.rzx(a, b, t),
            |s, a, b, t| s.rzx(a, b, t),
        ),
        (
            "rzy",
            |s, a, b, t| s.rzy(a, b, t),
            |s, a, b, t| s.rzy(a, b, t),
        ),
        (
            "rzz",
            |s, a, b, t| s.rzz(a, b, t),
            |s, a, b, t| s.rzz(a, b, t),
        ),
    ];
    for (name, old, new) in gates {
        bench_mut(
            c,
            &format!("rotation_two/{name}"),
            move |s| old(s, 2, 5, THETA),
            move |s| new(s, 2, 5, THETA),
        );
    }
}

fn batch(c: &mut Criterion) {
    let gates: [(&str, ManyOld, ManyNew); 9] = [
        (
            "rxx",
            |s, p, t| s.rxx_many(p, t),
            |s, p, t| s.rxx_many(p, t),
        ),
        (
            "rxy",
            |s, p, t| s.rxy_many(p, t),
            |s, p, t| s.rxy_many(p, t),
        ),
        (
            "rxz",
            |s, p, t| s.rxz_many(p, t),
            |s, p, t| s.rxz_many(p, t),
        ),
        (
            "ryx",
            |s, p, t| s.ryx_many(p, t),
            |s, p, t| s.ryx_many(p, t),
        ),
        (
            "ryy",
            |s, p, t| s.ryy_many(p, t),
            |s, p, t| s.ryy_many(p, t),
        ),
        (
            "ryz",
            |s, p, t| s.ryz_many(p, t),
            |s, p, t| s.ryz_many(p, t),
        ),
        (
            "rzx",
            |s, p, t| s.rzx_many(p, t),
            |s, p, t| s.rzx_many(p, t),
        ),
        (
            "rzy",
            |s, p, t| s.rzy_many(p, t),
            |s, p, t| s.rzy_many(p, t),
        ),
        (
            "rzz",
            |s, p, t| s.rzz_many(p, t),
            |s, p, t| s.rzz_many(p, t),
        ),
    ];
    for (name, old, new) in gates {
        bench_mut(
            c,
            &format!("rotation_two_batch/{name}"),
            move |s| old(s, PAIRS, THETA),
            move |s| new(s, PAIRS, THETA),
        );
    }
}
