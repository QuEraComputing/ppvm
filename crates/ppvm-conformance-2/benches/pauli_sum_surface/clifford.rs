// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_traits::traits::{
    Clifford as OldClifford, CliffordBatch as OldCliffordBatch,
    CliffordExtensions as OldExtensions, CliffordExtensionsBatch as OldExtensionsBatch,
};
use ppvm_traits_2::{
    Clifford as NewClifford, CliffordBatch as NewCliffordBatch,
    CliffordExtensions as NewExtensions, CliffordExtensionsBatch as NewExtensionsBatch,
};

use super::{NewSum, OldSum, bench_mut};

type OneOld = fn(&mut OldSum, usize);
type OneNew = fn(&mut NewSum, usize);
type TwoOld = fn(&mut OldSum, usize, usize);
type TwoNew = fn(&mut NewSum, usize, usize);
type ManyOld = fn(&mut OldSum, &[usize]);
type ManyNew = fn(&mut NewSum, &[usize]);
type PairsOld = fn(&mut OldSum, &[(usize, usize)]);
type PairsNew = fn(&mut NewSum, &[(usize, usize)]);

pub fn bench(c: &mut Criterion) {
    singles(c);
    aliases(c);
    batches(c);
}

fn singles(c: &mut Criterion) {
    let one: [(&str, OneOld, OneNew); 10] = [
        ("clifford/x", |s, q| s.x(q), |s, q| s.x(q)),
        ("clifford/y", |s, q| s.y(q), |s, q| s.y(q)),
        ("clifford/z", |s, q| s.z(q), |s, q| s.z(q)),
        ("clifford/h", |s, q| s.h(q), |s, q| s.h(q)),
        ("clifford/s", |s, q| s.s(q), |s, q| s.s(q)),
        ("clifford/s_dag", |s, q| s.s_dag(q), |s, q| s.s_dag(q)),
        ("clifford/sqrt_x", |s, q| s.sqrt_x(q), |s, q| s.sqrt_x(q)),
        (
            "clifford/sqrt_x_dag",
            |s, q| s.sqrt_x_dag(q),
            |s, q| s.sqrt_x_dag(q),
        ),
        ("clifford/sqrt_y", |s, q| s.sqrt_y(q), |s, q| s.sqrt_y(q)),
        (
            "clifford/sqrt_y_dag",
            |s, q| s.sqrt_y_dag(q),
            |s, q| s.sqrt_y_dag(q),
        ),
    ];
    for (name, old, new) in one {
        bench_mut(c, name, move |s| old(s, 3), move |s| new(s, 3));
    }
    let two: [(&str, TwoOld, TwoNew); 3] = [
        (
            "clifford/cnot",
            |s, a, b| s.cnot(a, b),
            |s, a, b| s.cnot(a, b),
        ),
        ("clifford/cz", |s, a, b| s.cz(a, b), |s, a, b| s.cz(a, b)),
        ("clifford/cy", |s, a, b| s.cy(a, b), |s, a, b| s.cy(a, b)),
    ];
    for (name, old, new) in two {
        bench_mut(c, name, move |s| old(s, 2, 5), move |s| new(s, 2, 5));
    }
}

fn aliases(c: &mut Criterion) {
    let aliases: [(&str, TwoOld, TwoNew); 4] = [
        (
            "clifford/cx_alias",
            |s, a, b| s.cx(a, b),
            |s, a, b| s.cx(a, b),
        ),
        (
            "clifford/zcx_alias",
            |s, a, b| s.zcx(a, b),
            |s, a, b| s.zcx(a, b),
        ),
        (
            "clifford/zcz_alias",
            |s, a, b| s.zcz(a, b),
            |s, a, b| s.zcz(a, b),
        ),
        (
            "clifford/zcy_alias",
            |s, a, b| s.zcy(a, b),
            |s, a, b| s.zcy(a, b),
        ),
    ];
    for (name, old, new) in aliases {
        bench_mut(c, name, move |s| old(s, 1, 6), move |s| new(s, 1, 6));
    }
}

fn batches(c: &mut Criterion) {
    const TARGETS: &[usize] = &[0, 2, 4, 6];
    const PAIRS: &[(usize, usize)] = &[(0, 1), (2, 3), (4, 5), (6, 7)];
    let one: [(&str, ManyOld, ManyNew); 10] = [
        ("clifford_batch/x", |s, q| s.x_many(q), |s, q| s.x_many(q)),
        ("clifford_batch/y", |s, q| s.y_many(q), |s, q| s.y_many(q)),
        ("clifford_batch/z", |s, q| s.z_many(q), |s, q| s.z_many(q)),
        ("clifford_batch/h", |s, q| s.h_many(q), |s, q| s.h_many(q)),
        ("clifford_batch/s", |s, q| s.s_many(q), |s, q| s.s_many(q)),
        (
            "clifford_batch/s_dag",
            |s, q| s.s_dag_many(q),
            |s, q| s.s_dag_many(q),
        ),
        (
            "clifford_batch/sqrt_x",
            |s, q| s.sqrt_x_many(q),
            |s, q| s.sqrt_x_many(q),
        ),
        (
            "clifford_batch/sqrt_x_dag",
            |s, q| s.sqrt_x_dag_many(q),
            |s, q| s.sqrt_x_dag_many(q),
        ),
        (
            "clifford_batch/sqrt_y",
            |s, q| s.sqrt_y_many(q),
            |s, q| s.sqrt_y_many(q),
        ),
        (
            "clifford_batch/sqrt_y_dag",
            |s, q| s.sqrt_y_dag_many(q),
            |s, q| s.sqrt_y_dag_many(q),
        ),
    ];
    for (name, old, new) in one {
        bench_mut(c, name, move |s| old(s, TARGETS), move |s| new(s, TARGETS));
    }
    let pairs: [(&str, PairsOld, PairsNew); 3] = [
        (
            "clifford_batch/cnot",
            |s, p| s.cnot_many(p),
            |s, p| s.cnot_many(p),
        ),
        (
            "clifford_batch/cz",
            |s, p| s.cz_many(p),
            |s, p| s.cz_many(p),
        ),
        (
            "clifford_batch/cy",
            |s, p| s.cy_many(p),
            |s, p| s.cy_many(p),
        ),
    ];
    for (name, old, new) in pairs {
        bench_mut(c, name, move |s| old(s, PAIRS), move |s| new(s, PAIRS));
    }
}
