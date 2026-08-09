// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The Phase-6 perf gate: the **columnar** `ColumnStore` (structure-of-arrays)
//! backend against the shipped `HashMapStore` (hash-join) backend, per *op
//! class*, plus one end-to-end workload.
//!
//! Both backends are the same `-2` engine (`Sum<S, P>`) over the same key type,
//! the same coefficient, the same policy and the same capacity hint — the ONLY
//! variable is `S`. Storage is pinned to `[u8; 8]` on both sides for the same
//! reason `pauli_sum_bench` pins it: a
//! storage-width codegen delta would otherwise fold into the layout ratio.
//!
//! The op classes are chosen to separate where SoA can win from where it must
//! lose (design: §"Backends are containers; columnar is expressible from day
//! one"):
//!
//! * **coefficient-only passes** — `scale` (one contiguous `*=`), `reduce` (the
//!   prefix-sum compaction), the `CoefficientThreshold` retain scan. The column
//!   store reads `len` contiguous `f64`s; the hash map walks *buckets*, key and
//!   all. This is where the layout is supposed to pay.
//! * **key-reading in-place passes** — `x` (sign flip), `pauli_error`
//!   (`scale_by_key`), the `MaxPauliWeight` retain scan. The callbacks take
//!   `&K`, so the columnar backend must **materialize** each key from its planes;
//!   the hash map hands out a borrow of the stored key.
//! * **hash-join paths** — `h`/`cnot` (re-key), `rx` (rotation branch merge),
//!   `overlap`, `from_terms`. Both backends probe by key here; SoA is expected to
//!   lose or draw.
//! * **end-to-end** — the TFIM Trotter propagation, which is a weighted mix of
//!   all three.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use ppvm_conformance_2::{random_terms, seeded_rng};
use ppvm_pauli_sum_2::{
    CoefficientThreshold, ColumnStore, CombinedPolicy, HashMapStore, MaxPauliWeight, NoPolicy,
    PauliWord, Policy, Sum,
};
use ppvm_traits_2::{Clifford, PauliError, RotationOne, RotationTwo};

type Key = PauliWord<[u8; 8]>;
type HashSum<P = NoPolicy> = Sum<HashMapStore<Key, f64>, P>;
type ColSum<P = NoPolicy> = Sum<ColumnStore<Key, f64>, P>;

/// Qubit width for the moderate-support op-class targets.
const N: usize = 16;
/// Support size for the op-class targets.
const TERMS: usize = 1000;

fn terms(seed: u64, count: usize, n: usize) -> Vec<(String, f64)> {
    let mut rng = seeded_rng(seed);
    random_terms(&mut rng, n, count)
}

fn hash_sum<P: Policy<Key, f64>>(n: usize, policy: P, t: &[(String, f64)]) -> HashSum<P> {
    let mut s = HashSum::with_capacity(n, policy, n * 10);
    for (w, c) in t {
        s += (Key::from(w.as_str()), *c);
    }
    s
}

fn col_sum<P: Policy<Key, f64>>(n: usize, policy: P, t: &[(String, f64)]) -> ColSum<P> {
    let mut s = ColSum::with_capacity(n, policy, n * 10);
    for (w, c) in t {
        s += (Key::from(w.as_str()), *c);
    }
    s
}

/// Emit one `hash`/`column` pair of timings for an in-place, support-preserving
/// op (an involution or a pure coefficient pass), so the two are measured on the
/// same support built the same way.
macro_rules! ab {
    ($group:expr, $name:literal, $t:expr, |$s:ident| $body:expr) => {{
        let mut $s = hash_sum(N, NoPolicy, $t);
        $group.bench_function(concat!("hash/", $name), |b| b.iter(|| $body));
        let mut $s = col_sum(N, NoPolicy, $t);
        $group.bench_function(concat!("column/", $name), |b| b.iter(|| $body));
    }};
}

// ---------------------------------------------------------------------------
// 1. Coefficient-only passes — where SoA is supposed to pay.
// ---------------------------------------------------------------------------

fn bench_scale(c: &mut Criterion) {
    let mut g = c.benchmark_group("column_store/scale");
    let t = terms(3, TERMS, N);
    ab!(g, "scale", &t, |s| s.scale(black_box(&1.0000001)));
    g.finish();
}

fn bench_reduce(c: &mut Criterion) {
    let mut g = c.benchmark_group("column_store/reduce");
    let t = terms(4, TERMS, N);
    // Nothing is zero, so `reduce` is a pure scan with no removal on either side
    // (the compaction's write index tracks the read index).
    ab!(g, "reduce", &t, |s| s.reduce());
    g.finish();
}

fn bench_truncate_threshold(c: &mut Criterion) {
    let mut g = c.benchmark_group("column_store/truncate_threshold");
    let t = terms(5, TERMS, N);
    // A threshold below every coefficient: the retain scan runs in full and keeps
    // everything, so this times the *scan*, not the removal.
    let policy = CoefficientThreshold { threshold: 1e-9 };

    let mut h = hash_sum(N, policy, &t);
    g.bench_function("hash/truncate", |b| b.iter(|| h.truncate()));
    let mut c2 = col_sum(N, policy, &t);
    g.bench_function("column/truncate", |b| b.iter(|| c2.truncate()));
    g.finish();
}

fn bench_truncate_weight(c: &mut Criterion) {
    let mut g = c.benchmark_group("column_store/truncate_weight");
    let t = terms(6, TERMS, N);
    // The weight predicate reads the KEY, so the columnar side pays a key
    // materialization per surviving term that the hash map does not.
    let policy = MaxPauliWeight(N);

    let mut h = hash_sum(N, policy, &t);
    g.bench_function("hash/truncate", |b| b.iter(|| h.truncate()));
    let mut c2 = col_sum(N, policy, &t);
    g.bench_function("column/truncate", |b| b.iter(|| c2.truncate()));

    // The `usize::MAX` disable sentinel must be ~free on BOTH backends (the
    // policy skips the retain pass entirely).
    // `MaxPauliWeight::default()` IS the `usize::MAX` sentinel, alone: the policy
    // must return without touching the support on either backend.
    let disabled = MaxPauliWeight::default();
    let mut h = hash_sum(N, disabled, &t);
    g.bench_function("hash/truncate_disabled", |b| b.iter(|| h.truncate()));
    let mut c2 = col_sum(N, disabled, &t);
    g.bench_function("column/truncate_disabled", |b| b.iter(|| c2.truncate()));
    g.finish();
}

// ---------------------------------------------------------------------------
// 2. Key-reading in-place passes — the `&K` callbacks that force a
//    materialization on the columnar side.
// ---------------------------------------------------------------------------

fn bench_sign_flip(c: &mut Criterion) {
    let mut g = c.benchmark_group("column_store/clifford_x");
    let t = terms(1, TERMS, N);
    ab!(g, "x", &t, |s| s.x(black_box(0)));
    g.finish();
}

fn bench_pauli_error(c: &mut Criterion) {
    let mut g = c.benchmark_group("column_store/pauli_error");
    let t = terms(7, TERMS, N);
    let hash_seed = hash_sum(N, NoPolicy, &t);
    g.bench_function("hash/pauli_error", |b| {
        b.iter_batched_ref(
            || hash_seed.clone(),
            |sum| {
                sum.pauli_error(
                    black_box(0),
                    [1e-4, 1e-4, 1e-4],
                    &mut ppvm_conformance_2::analytic_rng(),
                )
            },
            criterion::BatchSize::LargeInput,
        )
    });
    let col_seed = col_sum(N, NoPolicy, &t);
    g.bench_function("column/pauli_error", |b| {
        b.iter_batched_ref(
            || col_seed.clone(),
            |sum| {
                sum.pauli_error(
                    black_box(0),
                    [1e-4, 1e-4, 1e-4],
                    &mut ppvm_conformance_2::analytic_rng(),
                )
            },
            criterion::BatchSize::LargeInput,
        )
    });
    g.finish();
}

// ---------------------------------------------------------------------------
// 3. Hash-join paths — where SoA must still probe by key.
// ---------------------------------------------------------------------------

fn bench_clifford_h(c: &mut Criterion) {
    let mut g = c.benchmark_group("column_store/clifford_h");
    let t = terms(1, TERMS, N);
    ab!(g, "h", &t, |s| s.h(black_box(0)));
    g.finish();
}

fn bench_clifford_cnot(c: &mut Criterion) {
    let mut g = c.benchmark_group("column_store/clifford_cnot");
    let t = terms(1, TERMS, N);
    ab!(g, "cnot", &t, |s| s.cnot(black_box(0), black_box(1)));
    g.finish();
}

fn bench_rotation(c: &mut Criterion) {
    let mut g = c.benchmark_group("column_store/rx");
    let t = terms(8, TERMS, N);
    let mut closure_probe = col_sum(N, NoPolicy, &t);
    closure_probe.rx(0, 0.1);
    let anticommuting: Vec<_> = closure_probe
        .iter()
        .map(|(key, _)| key)
        .filter(|key| ppvm_traits_2::PauliBits::z_bit(key, 0))
        .collect();
    assert!(anticommuting.iter().all(|key| {
        let target = ppvm_traits_2::PauliBits::toggled_bits(key, 0, true, false);
        closure_probe.contains_key(&target)
    }));
    // The probe verifies that every generated branch key is already present.
    // Repeated `rx(θ)` therefore changes coefficients but does not grow support.
    ab!(g, "rx", &t, |s| s.rx(black_box(0), black_box(0.1)));
    g.finish();
}

fn bench_rotation_growth(c: &mut Criterion) {
    const SITES: usize = 8;
    let mut g = c.benchmark_group("column_store/rx_growth");
    let t = terms(8, TERMS, N);
    let hseed = hash_sum(N, NoPolicy, &t);
    let cseed = col_sum(N, NoPolicy, &t);

    g.bench_function("hash/rx_sweep", |b| {
        b.iter_batched_ref(
            || hseed.clone(),
            |sum| {
                for i in 0..SITES {
                    sum.rx(i, 0.1);
                }
            },
            criterion::BatchSize::LargeInput,
        )
    });
    g.bench_function("column/rx_sweep", |b| {
        b.iter_batched_ref(
            || cseed.clone(),
            |sum| {
                for i in 0..SITES {
                    sum.rx(i, 0.1);
                }
            },
            criterion::BatchSize::LargeInput,
        )
    });
    g.finish();
}

fn bench_overlap(c: &mut Criterion) {
    let mut g = c.benchmark_group("column_store/overlap");
    let a = terms(9, TERMS, N);
    let b2 = terms(10, TERMS, N);

    let ha = hash_sum(N, NoPolicy, &a);
    let hb = hash_sum(N, NoPolicy, &b2);
    g.bench_function("hash/overlap", |b| b.iter(|| black_box(ha.overlap(&hb))));

    let ca = col_sum(N, NoPolicy, &a);
    let cb = col_sum(N, NoPolicy, &b2);
    g.bench_function("column/overlap", |b| b.iter(|| black_box(ca.overlap(&cb))));
    g.finish();
}

fn bench_build(c: &mut Criterion) {
    let mut g = c.benchmark_group("column_store/build");
    let t = terms(2, TERMS, N);
    let keyed: Vec<(Key, f64)> = t.iter().map(|(w, c)| (Key::from(w.as_str()), *c)).collect();

    g.bench_function("hash/from_terms", |b| {
        b.iter(|| {
            let s: HashSum = HashSum::from_terms(N, keyed.iter().cloned());
            black_box(s.len())
        })
    });
    g.bench_function("column/from_terms", |b| {
        b.iter(|| {
            let s: ColSum = ColSum::from_terms(N, keyed.iter().cloned());
            black_box(s.len())
        })
    });
    g.finish();
}

// ---------------------------------------------------------------------------
// 4. End-to-end: the TFIM Trotter propagation (the headline integration
//    workload's shape), driven on both backends.
// ---------------------------------------------------------------------------

const TROTTER_N: usize = 12;
const TROTTER_STEPS: usize = 10;
const THETA_X: f64 = 0.1;
const THETA_ZZ: f64 = 0.0125;
const NOISE: [f64; 3] = [2.5e-5, 2.5e-5, 2.5e-5];

macro_rules! trotter_body {
    ($sum:ident) => {{
        for _ in 0..TROTTER_STEPS {
            for i in 0..TROTTER_N {
                $sum.pauli_error(i, NOISE, &mut ppvm_conformance_2::analytic_rng());
                $sum.truncate();
                $sum.rx(i, THETA_X);
                $sum.truncate();
            }
            for i in 0..TROTTER_N - 1 {
                $sum.pauli_error(i + 1, NOISE, &mut ppvm_conformance_2::analytic_rng());
                $sum.truncate();
                $sum.pauli_error(i, NOISE, &mut ppvm_conformance_2::analytic_rng());
                $sum.truncate();
                $sum.rzz(i, i + 1, THETA_ZZ);
                $sum.truncate();
            }
        }
        black_box($sum.len())
    }};
}

fn bench_trotter(c: &mut Criterion) {
    let mut g = c.benchmark_group("column_store/trotter_tfim_n12");
    let policy = CombinedPolicy(
        CoefficientThreshold { threshold: 1e-6 },
        MaxPauliWeight::default(),
    );
    let cap = TROTTER_N * TROTTER_N;

    let seed: Vec<(Key, f64)> = (0..TROTTER_N)
        .map(|i| {
            let mut w = Key::new(TROTTER_N);
            ppvm_traits_2::PauliBits::set_z_bit(&mut w, i, true);
            (w, 1.0)
        })
        .collect();

    let hash_seed = {
        let mut s = HashSum::with_capacity(TROTTER_N, policy, cap);
        for (w, c) in &seed {
            s += (*w, *c);
        }
        s
    };
    g.bench_function("hash/trotter", |b| {
        b.iter_batched_ref(
            || hash_seed.clone(),
            |sum| trotter_body!(sum),
            criterion::BatchSize::SmallInput,
        )
    });

    let col_seed = {
        let mut s = ColSum::with_capacity(TROTTER_N, policy, cap);
        for (w, c) in &seed {
            s += (*w, *c);
        }
        s
    };
    g.bench_function("column/trotter", |b| {
        b.iter_batched_ref(
            || col_seed.clone(),
            |sum| trotter_body!(sum),
            criterion::BatchSize::SmallInput,
        )
    });
    g.finish();
}

criterion_group!(
    benches,
    bench_scale,
    bench_reduce,
    bench_truncate_threshold,
    bench_truncate_weight,
    bench_sign_flip,
    bench_pauli_error,
    bench_clifford_h,
    bench_clifford_cnot,
    bench_rotation,
    bench_rotation_growth,
    bench_overlap,
    bench_build,
    bench_trotter,
);
criterion_main!(benches);
