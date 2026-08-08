// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_conformance_2::sym::{NewTerm, OldTerm};
use ppvm_traits::traits::RotationTwo as Old;
use ppvm_traits_2::RotationTwo as New;

use super::paired_args;

pub(super) fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("sym/surface/propagation/rotation_two");
    let pairs = &[(0, 1), (2, 3)][..];

    paired_args(
        &mut group,
        "generic_rotate_2",
        OldTerm::var(0),
        NewTerm::var(0),
        |s, a| Old::rotate_2(s, [1, 0], [0, 1], 0, 1, a),
        |s, a| New::rotate_2(s, [1, 0], [0, 1], 0, 1, a),
    );

    macro_rules! named {
        ($single:ident, $batch:ident) => {
            paired_args(
                &mut group,
                stringify!($single),
                OldTerm::var(0),
                NewTerm::var(0),
                |s, a| Old::$single(s, 0, 1, a),
                |s, a| New::$single(s, 0, 1, a),
            );
            paired_args(
                &mut group,
                concat!("batch_", stringify!($single)),
                OldTerm::var(0),
                NewTerm::var(0),
                |s, a| Old::$batch(s, pairs, a),
                |s, a| New::$batch(s, pairs, a),
            );
        };
    }

    named!(rxx, rxx_many);
    named!(rxy, rxy_many);
    named!(rxz, rxz_many);
    named!(ryx, ryx_many);
    named!(ryy, ryy_many);
    named!(ryz, ryz_many);
    named!(rzx, rzx_many);
    named!(rzy, rzy_many);
    named!(rzz, rzz_many);
    group.finish();
}
