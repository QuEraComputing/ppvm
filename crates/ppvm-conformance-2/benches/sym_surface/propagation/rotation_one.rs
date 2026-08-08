// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_conformance_2::sym::{NewTerm, OldTerm};
use ppvm_traits::{
    char::Pauli as OldPauli,
    traits::{RotXY as OldXY, RotationOne as Old},
};
use ppvm_traits_2::{Pauli as NewPauli, RotXY as NewXY, RotationOne as New};

use super::paired_args;

pub(super) fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("sym/surface/propagation/rotation_one");
    macro_rules! pair {
        ($name:literal, $old:expr, $new:expr) => {
            paired_args(
                &mut group,
                $name,
                OldTerm::var(0),
                NewTerm::var(0),
                $old,
                $new,
            );
        };
    }

    pair!(
        "generic_rotate_1",
        |s, a| Old::rotate_1(s, OldPauli::X, 0, a),
        |s, a| New::rotate_1(s, NewPauli::X, 0, a)
    );
    pair!("rx", |s, a| Old::rx(s, 0, a), |s, a| New::rx(s, 0, a));
    pair!("ry", |s, a| Old::ry(s, 0, a), |s, a| New::ry(s, 0, a));
    pair!("rz", |s, a| Old::rz(s, 0, a), |s, a| New::rz(s, 0, a));

    let targets = &[0, 2, 4][..];
    pair!("batch_rx", |s, a| Old::rx_many(s, targets, a), |s, a| {
        New::rx_many(s, targets, a)
    });
    pair!("batch_ry", |s, a| Old::ry_many(s, targets, a), |s, a| {
        New::ry_many(s, targets, a)
    });
    pair!("batch_rz", |s, a| Old::rz_many(s, targets, a), |s, a| {
        New::rz_many(s, targets, a)
    });

    paired_args(
        &mut group,
        "rot_xy_r",
        (OldTerm::from(0.4), OldTerm::var(0)),
        (NewTerm::from(0.4), NewTerm::var(0)),
        |s, (axis, theta)| OldXY::r(s, 0, axis, theta),
        |s, (axis, theta)| NewXY::r(s, 0, axis, theta),
    );
    group.finish();
}
