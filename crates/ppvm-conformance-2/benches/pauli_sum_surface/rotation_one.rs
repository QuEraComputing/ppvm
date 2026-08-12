// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_traits::char::Pauli as OldPauli;
use ppvm_traits::traits::{RotXY as OldRotXY, RotationOne as OldRotationOne};
use ppvm_traits_2::{Pauli as NewPauli, RotXY as NewRotXY, RotationOne as NewRotationOne};

use super::{NewSum, OldSum, bench_mut};

type RotOld = fn(&mut OldSum, usize, f64);
type RotNew = fn(&mut NewSum, usize, f64);
type ManyOld = fn(&mut OldSum, &[usize], f64);
type ManyNew = fn(&mut NewSum, &[usize], f64);

const THETA: f64 = 0.37;
const TARGETS: &[usize] = &[0, 2, 4, 6];

pub fn bench(c: &mut Criterion) {
    generic(c);
    named(c);
    batch(c);
    bench_mut(
        c,
        "rotation_one/rot_xy_r",
        |s| s.r(3, 0.61, THETA),
        |s| s.r(3, 0.61, THETA),
    );
}

fn generic(c: &mut Criterion) {
    let axes = [
        ("i", OldPauli::I, NewPauli::I),
        ("x", OldPauli::X, NewPauli::X),
        ("y", OldPauli::Y, NewPauli::Y),
        ("z", OldPauli::Z, NewPauli::Z),
    ];
    for (name, old_axis, new_axis) in axes {
        bench_mut(
            c,
            &format!("rotation_one/rotate_1_{name}"),
            move |s| s.rotate_1(old_axis, 3, THETA),
            move |s| s.rotate_1(new_axis, 3, THETA),
        );
    }
}

fn named(c: &mut Criterion) {
    let gates: [(&str, RotOld, RotNew); 3] = [
        (
            "rotation_one/rx",
            |s, q, t| s.rx(q, t),
            |s, q, t| s.rx(q, t),
        ),
        (
            "rotation_one/ry",
            |s, q, t| s.ry(q, t),
            |s, q, t| s.ry(q, t),
        ),
        (
            "rotation_one/rz",
            |s, q, t| s.rz(q, t),
            |s, q, t| s.rz(q, t),
        ),
    ];
    for (name, old, new) in gates {
        bench_mut(
            c,
            name,
            move |s| old(s, 3, THETA),
            move |s| new(s, 3, THETA),
        );
    }
}

fn batch(c: &mut Criterion) {
    let gates: [(&str, ManyOld, ManyNew); 3] = [
        (
            "rotation_one_batch/rx",
            |s, q, t| s.rx_many(q, t),
            |s, q, t| s.rx_many(q, t),
        ),
        (
            "rotation_one_batch/ry",
            |s, q, t| s.ry_many(q, t),
            |s, q, t| s.ry_many(q, t),
        ),
        (
            "rotation_one_batch/rz",
            |s, q, t| s.rz_many(q, t),
            |s, q, t| s.rz_many(q, t),
        ),
    ];
    for (name, old, new) in gates {
        bench_mut(
            c,
            name,
            move |s| old(s, TARGETS, THETA),
            move |s| new(s, TARGETS, THETA),
        );
    }
}
