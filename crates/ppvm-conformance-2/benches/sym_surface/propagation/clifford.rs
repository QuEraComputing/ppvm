// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_traits::traits::{
    Clifford as Old, CliffordBatch as OldBatch, CliffordExtensions as OldExt,
    CliffordExtensionsBatch as OldExtBatch,
};
use ppvm_traits_2::{
    Clifford as New, CliffordBatch as NewBatch, CliffordExtensions as NewExt,
    CliffordExtensionsBatch as NewExtBatch,
};

use super::paired;

pub(super) fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("sym/surface/propagation/clifford");
    macro_rules! pair {
        ($name:literal, $old:expr, $new:expr) => {
            paired(&mut group, $name, $old, $new);
        };
    }

    pair!("x", |s| Old::x(s, 0), |s| New::x(s, 0));
    pair!("y", |s| Old::y(s, 0), |s| New::y(s, 0));
    pair!("z", |s| Old::z(s, 0), |s| New::z(s, 0));
    pair!("h", |s| Old::h(s, 0), |s| New::h(s, 0));
    pair!("s", |s| Old::s(s, 0), |s| New::s(s, 0));
    pair!("cnot", |s| Old::cnot(s, 0, 1), |s| New::cnot(s, 0, 1));
    pair!("cz", |s| Old::cz(s, 0, 1), |s| New::cz(s, 0, 1));

    pair!("alias_cx", |s| Old::cx(s, 0, 1), |s| New::cx(s, 0, 1));
    pair!("alias_zcx", |s| Old::zcx(s, 0, 1), |s| New::zcx(s, 0, 1));
    pair!("alias_zcz", |s| Old::zcz(s, 0, 1), |s| New::zcz(s, 0, 1));

    pair!("s_dag", |s| OldExt::s_dag(s, 0), |s| NewExt::s_dag(s, 0));
    pair!("sqrt_x", |s| OldExt::sqrt_x(s, 0), |s| NewExt::sqrt_x(s, 0));
    pair!("sqrt_x_dag", |s| OldExt::sqrt_x_dag(s, 0), |s| {
        NewExt::sqrt_x_dag(s, 0)
    });
    pair!("sqrt_y", |s| OldExt::sqrt_y(s, 0), |s| NewExt::sqrt_y(s, 0));
    pair!("sqrt_y_dag", |s| OldExt::sqrt_y_dag(s, 0), |s| {
        NewExt::sqrt_y_dag(s, 0)
    });
    pair!("cy", |s| OldExt::cy(s, 0, 1), |s| NewExt::cy(s, 0, 1));
    pair!("alias_zcy", |s| OldExt::zcy(s, 0, 1), |s| NewExt::zcy(
        s, 0, 1
    ));

    let targets = &[0, 2, 4][..];
    let pairs = &[(0, 1), (2, 3)][..];
    pair!("batch_x", |s| OldBatch::x_many(s, targets), |s| {
        NewBatch::x_many(s, targets)
    });
    pair!("batch_y", |s| OldBatch::y_many(s, targets), |s| {
        NewBatch::y_many(s, targets)
    });
    pair!("batch_z", |s| OldBatch::z_many(s, targets), |s| {
        NewBatch::z_many(s, targets)
    });
    pair!("batch_h", |s| OldBatch::h_many(s, targets), |s| {
        NewBatch::h_many(s, targets)
    });
    pair!("batch_s", |s| OldBatch::s_many(s, targets), |s| {
        NewBatch::s_many(s, targets)
    });
    pair!("batch_cnot", |s| OldBatch::cnot_many(s, pairs), |s| {
        NewBatch::cnot_many(s, pairs)
    });
    pair!("batch_cz", |s| OldBatch::cz_many(s, pairs), |s| {
        NewBatch::cz_many(s, pairs)
    });
    pair!(
        "batch_s_dag",
        |s| OldExtBatch::s_dag_many(s, targets),
        |s| NewExtBatch::s_dag_many(s, targets)
    );
    pair!(
        "batch_sqrt_x",
        |s| OldExtBatch::sqrt_x_many(s, targets),
        |s| NewExtBatch::sqrt_x_many(s, targets)
    );
    pair!(
        "batch_sqrt_x_dag",
        |s| OldExtBatch::sqrt_x_dag_many(s, targets),
        |s| NewExtBatch::sqrt_x_dag_many(s, targets)
    );
    pair!(
        "batch_sqrt_y",
        |s| OldExtBatch::sqrt_y_many(s, targets),
        |s| NewExtBatch::sqrt_y_many(s, targets)
    );
    pair!(
        "batch_sqrt_y_dag",
        |s| OldExtBatch::sqrt_y_dag_many(s, targets),
        |s| NewExtBatch::sqrt_y_dag_many(s, targets)
    );
    pair!("batch_cy", |s| OldExtBatch::cy_many(s, pairs), |s| {
        NewExtBatch::cy_many(s, pairs)
    });
    group.finish();
}
