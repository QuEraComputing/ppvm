// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_traits::traits::{Clifford as _, CliffordExtensions as _};
use ppvm_traits_2::{Clifford as _, CliffordExtensions as _, Indexable};
use std::hint::black_box;

use super::common::*;

macro_rules! gate1 {
    ($g:expr, $new:expr, $old:expr, $name:literal, $method:ident, $q:expr) => {{
        let mut nc = $new.clone();
        let mut oc = $old.clone();
        nc.$method($q);
        oc.$method($q);
        assert_eq!(nc.to_string(), oc.to_string(), $name);
        $g.bench_function(concat!("new/", $name), |b| {
            b.iter_batched(
                || $new.clone(),
                |mut w| {
                    w.$method(black_box($q));
                    black_box(w)
                },
                BatchSize::SmallInput,
            )
        });
        $g.bench_function(concat!("old/", $name), |b| {
            b.iter_batched(
                || $old.clone(),
                |mut w| {
                    w.$method(black_box($q));
                    black_box(w)
                },
                BatchSize::SmallInput,
            )
        });
    }};
}

macro_rules! gate2 {
    ($g:expr, $new:expr, $old:expr, $name:literal, $method:ident, $a:expr, $b:expr) => {{
        let mut nc = $new.clone();
        let mut oc = $old.clone();
        nc.$method($a, $b);
        oc.$method($a, $b);
        assert_eq!(nc.to_string(), oc.to_string(), $name);
        $g.bench_function(concat!("new/", $name), |bench| {
            bench.iter_batched(
                || $new.clone(),
                |mut w| {
                    w.$method(black_box($a), black_box($b));
                    black_box(w)
                },
                BatchSize::SmallInput,
            )
        });
        $g.bench_function(concat!("old/", $name), |bench| {
            bench.iter_batched(
                || $old.clone(),
                |mut w| {
                    w.$method(black_box($a), black_box($b));
                    black_box(w)
                },
                BatchSize::SmallInput,
            )
        });
    }};
}

pub fn bench(c: &mut Criterion) {
    let mut g = c.benchmark_group("word_surface/ordinary/clifford/256");
    let s = ordinary_string(WIDTH);
    let new = NewWord::from(s.as_str());
    let old = OldWord::from(s.as_str());
    black_box(new.key_hash());
    const Q: usize = 126;
    const R: usize = 127;
    gate1!(g, new, old, "x", x, Q);
    gate1!(g, new, old, "y", y, Q);
    gate1!(g, new, old, "z", z, Q);
    gate1!(g, new, old, "h", h, Q);
    gate1!(g, new, old, "s", s, Q);
    gate2!(g, new, old, "cnot", cnot, Q, R);
    gate2!(g, new, old, "cz", cz, Q, R);
    gate2!(g, new, old, "cx_alias", cx, Q, R);
    gate2!(g, new, old, "zcx_alias", zcx, Q, R);
    gate2!(g, new, old, "zcz_alias", zcz, Q, R);
    gate1!(g, new, old, "s_dag", s_dag, Q);
    gate1!(g, new, old, "sqrt_x", sqrt_x, Q);
    gate1!(g, new, old, "sqrt_x_dag", sqrt_x_dag, Q);
    gate1!(g, new, old, "sqrt_y", sqrt_y, Q);
    gate1!(g, new, old, "sqrt_y_dag", sqrt_y_dag, Q);
    gate2!(g, new, old, "cy", cy, Q, R);
    gate2!(g, new, old, "zcy_alias", zcy, Q, R);
    g.finish();
}
