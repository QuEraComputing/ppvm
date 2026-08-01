// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Comparative benchmarks: new `ppvm-pauli-sum-2::PauliSum<f64>` vs the old
//! `ppvm-pauli-sum::PauliSum<f64>` on the hot sparse-sum paths, so the refactor's
//! Phase-3 performance gate reads off as a new/old ratio per target.
//!
//! Targets (design: `traits-2-implementation-plan.md` Phase 3 perf gate):
//! * **Clifford gate on a moderate-support sum** — `h(q)` and `cnot(c, t)` on a
//!   ~1000-term `PauliSum`. This is the Phase-2 watch item come due: the new word
//!   dropped `Copy` (lazy `OnceLock` hash cache), so the per-term key clones on the
//!   re-key path may cost against the old `Copy` word. The ratio is the gate here.
//! * **`accumulate_batch` + `reduce`** of a produced batch — `from_terms` (new) vs
//!   the old `+=` build.
//! * **`scale`** by a constant.
//! * **`Pair` overlap** of two sums.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use ppvm_conformance_2::{
    NewKey, build_new_sum, build_old_sum, random_terms, reduce_old, seeded_rng,
};
use ppvm_pauli_sum_2::PauliWord as NewPauliWord;

use ppvm_traits::traits::Clifford as OldClifford;
use ppvm_traits::traits::{PauliError as OldPauliError, RotationOne as OldRotationOne};
use ppvm_traits_2::Clifford as NewClifford;
use ppvm_traits_2::{PauliError as NewPauliError, RotationOne as NewRotationOne};

/// Qubit width for the moderate-support sums (fits both `u64` and `[u8; 8]`).
const N: usize = 16;
/// Support size for the "moderate support" Clifford/scale/overlap targets.
const TERMS: usize = 1000;

/// A shared ~`TERMS`-term term list, seeded so new and old build identical sums.
fn terms(seed: u64, count: usize) -> Vec<(String, f64)> {
    let mut rng = seeded_rng(seed);
    random_terms(&mut rng, N, count)
}

// ---------------------------------------------------------------------------
// 1. Clifford gate on a ~1000-term sum — the Copy-drop watch item.
// ---------------------------------------------------------------------------

fn bench_clifford_h(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_sum/clifford_h");
    let t = terms(1, TERMS);

    // `h` is an involution, so the support size is invariant across iterations —
    // each timed `h(0)` is a full re-key of the whole ~1000-term support (clone
    // every key, drain the sign, accumulate, reduce).
    let mut new = build_new_sum(N, &t);
    g.bench_function("new/h", |b| b.iter(|| new.h(black_box(0))));

    let mut old = build_old_sum(N, &t);
    g.bench_function("old/h", |b| b.iter(|| old.h(black_box(0))));

    g.finish();
}

fn bench_clifford_x(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_sum/clifford_x");
    let t = terms(1, TERMS);

    // `X` is a pure-sign gate: the Pauli word is fixed and only each term's
    // coefficient flips `±1` (`XPX = (−1)^z P`). It is an involution, so every
    // timed `x(0)` is a full in-place sign pass over the whole ~1000-term support.
    // This is the pure-sign fast path (`sign_flip_by_key`, no map rebuild) vs the
    // old crate's in-place `scale`.
    let mut new = build_new_sum(N, &t);
    g.bench_function("new/x", |b| b.iter(|| new.x(black_box(0))));

    let mut old = build_old_sum(N, &t);
    g.bench_function("old/x", |b| b.iter(|| old.x(black_box(0))));

    g.finish();
}

fn bench_clifford_cnot(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_sum/clifford_cnot");
    let t = terms(1, TERMS);

    let mut new = build_new_sum(N, &t);
    g.bench_function("new/cnot", |b| {
        b.iter(|| new.cnot(black_box(0), black_box(1)))
    });

    let mut old = build_old_sum(N, &t);
    g.bench_function("old/cnot", |b| {
        b.iter(|| old.cnot(black_box(0), black_box(1)))
    });

    g.finish();
}

// ---------------------------------------------------------------------------
// 2. accumulate_batch + reduce of a produced batch.
// ---------------------------------------------------------------------------

fn bench_build_batch(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_sum/build_batch");
    let t = terms(2, TERMS);

    // New: `from_terms` = accumulate_batch + reduce over the pre-keyed batch.
    let new_terms: Vec<(NewPauliWord, f64)> = t
        .iter()
        .map(|(w, c)| (NewKey::from(w.as_str()), *c))
        .collect();
    g.bench_function("new/from_terms", |b| {
        b.iter(|| {
            let s = ppvm_conformance_2::NewSum::from_terms(N, new_terms.iter().cloned());
            black_box(s.len());
        })
    });

    // Old: build the map by `+=` then `truncate()` (the reduce equivalent).
    g.bench_function("old/build_add_assign", |b| {
        b.iter(|| {
            let mut s = build_old_sum(N, &t);
            reduce_old(&mut s);
            black_box(s.len());
        })
    });

    g.finish();
}

// ---------------------------------------------------------------------------
// 3. scale by a constant.
// ---------------------------------------------------------------------------

fn bench_scale(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_sum/scale");
    let t = terms(3, TERMS);

    let mut new = build_new_sum(N, &t);
    g.bench_function("new/scale", |b| b.iter(|| new.scale(black_box(&1.0000001))));

    let mut old = build_old_sum(N, &t);
    g.bench_function("old/scale", |b| b.iter(|| old *= black_box(1.0000001)));

    g.finish();
}

// ---------------------------------------------------------------------------
// 4. Pair overlap of two sums.
// ---------------------------------------------------------------------------

fn bench_overlap(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_sum/overlap");
    let a = terms(4, TERMS);
    let b_terms = terms(5, TERMS);

    let new_a = build_new_sum(N, &a);
    let new_b = build_new_sum(N, &b_terms);
    g.bench_function("new/overlap", |b| {
        b.iter(|| black_box(new_a.overlap(black_box(&new_b))))
    });

    let old_a = build_old_sum(N, &a);
    let old_b = build_old_sum(N, &b_terms);
    g.bench_function("old/overlap", |b| {
        b.iter(|| black_box(old_a.overlap(black_box(&old_b))))
    });

    g.finish();
}

// ---------------------------------------------------------------------------
// 5. rx(q, θ) on a ~1000-term sum — the branching producer (the non-Clifford hot
//    path: fan-out + accumulate_batch + reduce), new vs old.
// ---------------------------------------------------------------------------

fn bench_rotation_rx(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_sum/rotation_rx");
    let t = terms(6, TERMS);

    // A generic angle (both sin and cos nonzero) so every anticommuting term
    // fans out to two, and the branch keys merge back through accumulate_batch +
    // reduce. rx is NOT an involution, so we rebuild the sum each iteration to
    // keep the timed workload's support size stable (otherwise the support would
    // grow/rotate across iterations and the ratio would drift).
    const THETA: f64 = 0.7;

    let new0 = build_new_sum(N, &t);
    g.bench_function("new/rx", |b| {
        b.iter_batched_ref(
            || new0.clone(),
            |s| s.rx(black_box(0), black_box(THETA)),
            criterion::BatchSize::SmallInput,
        )
    });

    let old0 = build_old_sum(N, &t);
    g.bench_function("old/rx", |b| {
        b.iter_batched_ref(
            || old0.clone(),
            |s| s.rx(black_box(0), black_box(THETA)),
            criterion::BatchSize::SmallInput,
        )
    });

    g.finish();
}

// ---------------------------------------------------------------------------
// 6. pauli_error(q, [..]) on a ~1000-term sum — the diagonal in-place scale,
//    new vs old.
// ---------------------------------------------------------------------------

fn bench_pauli_error(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_sum/pauli_error");
    let t = terms(7, TERMS);

    // A diagonal channel is idempotent in shape (keys are fixed; only the ±λ
    // scale changes), so repeated application over the same map is a stable,
    // support-invariant in-place pass — no rebuild needed.
    const P: [f64; 3] = [0.05, 0.1, 0.15];

    let mut new = build_new_sum(N, &t);
    g.bench_function("new/pauli_error", |b| {
        b.iter(|| new.pauli_error(black_box(0), black_box(P)))
    });

    let mut old = build_old_sum(N, &t);
    g.bench_function("old/pauli_error", |b| {
        b.iter(|| old.pauli_error(black_box(0), black_box(P)))
    });

    g.finish();
}

criterion_group!(
    benches,
    bench_clifford_h,
    bench_clifford_x,
    bench_clifford_cnot,
    bench_build_batch,
    bench_scale,
    bench_overlap,
    bench_rotation_rx,
    bench_pauli_error,
);
criterion_main!(benches);
