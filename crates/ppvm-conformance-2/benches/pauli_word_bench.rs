// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Comparative benchmarks: new `ppvm-pauli-word-2::PauliWord` vs the old
//! `ppvm-pauli-word` reference on the hot Pauli-word paths, so the refactor's
//! performance gate can be read off as a new/old ratio per target.
//!
//! Targets (design: `traits-2-implementation-plan.md` Phase 2 perf gate):
//! * Pauli×Pauli product — new `KeyProduct::key_mul` vs the old phased
//!   `MulAssign` (the kernel the new one was ported from).
//! * one Clifford conjugation — `cnot` on the bare bit map, new vs old.
//! * `key_hash()` cold (first compute) vs warm (cached `OnceLock`) — the lazy
//!   hashing the design trades `Copy` for.
//! * `weight()` — the fused popcount, new vs old.

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use ppvm_pauli_word_2::PauliWord as NewWord;
use ppvm_traits_2::{Clifford as NewClifford, Indexable, KeyProduct};

use ppvm_pauli_word::phase::PhasedPauliWord;
use ppvm_pauli_word::word::PauliWord as OldBareWord;
use ppvm_traits::traits::{Clifford as OldClifford, PauliWordTrait};

type New = NewWord<u64>;
type OldBare = OldBareWord<u64>;
type OldPhased = PhasedPauliWord<u64>;

const LHS: &str = "XYZIXYZIXYZIXYZI"; // 16 qubits
const RHS: &str = "ZZXYIIYXZXYZIXYZ";

fn bench_product(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_word/product");

    let na: New = LHS.into();
    let nb: New = RHS.into();
    g.bench_function("new/key_mul", |b| {
        b.iter(|| black_box(black_box(&na).key_mul(black_box(&nb))))
    });

    let oa: OldPhased = format!("+{LHS}").as_str().into();
    let ob: OldPhased = format!("+{RHS}").as_str().into();
    g.bench_function("old/phased_mul", |b| {
        b.iter(|| black_box(black_box(oa) * black_box(ob)))
    });

    g.finish();
}

fn bench_clifford_cnot(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_word/cnot");

    // One `cnot` conjugation per timed iteration, mutating a single persistent
    // word in place — no per-iteration clone/copy to confound the measurement.
    // The word stays a valid 16-qubit state (the CNOT bit map keeps it bounded).
    // This measures each crate's *true* per-op cost, including its hashing
    // policy: the old bare word eagerly `rehash()`es on every gate, the new word
    // only invalidates its lazy `OnceLock`.
    let mut nw: New = LHS.into();
    g.bench_function("new/cnot", |b| {
        b.iter(|| nw.cnot(black_box(0), black_box(1)))
    });

    let mut ow: OldBare = OldBare::from(LHS);
    g.bench_function("old/cnot", |b| {
        b.iter(|| ow.cnot(black_box(0), black_box(1)))
    });

    g.finish();
}

fn bench_key_hash(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_word/key_hash");

    // Cold: a never-hashed template, cloned in the untimed setup (clone of an
    // empty `OnceLock` stays empty) so the timed closure measures only the
    // first-compute of the digest.
    let cold_template: New = LHS.into();
    g.bench_function("new/cold", |b| {
        b.iter_batched(
            || cold_template.clone(),
            |w| {
                black_box(w.key_hash());
                w
            },
            BatchSize::SmallInput,
        )
    });

    // Warm: hash once up front, then re-read the cached value.
    let warm: New = LHS.into();
    black_box(warm.key_hash());
    g.bench_function("new/warm", |b| {
        b.iter(|| black_box(black_box(&warm).key_hash()))
    });

    g.finish();
}

fn bench_weight(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_word/weight");

    let nw: New = LHS.into();
    g.bench_function("new/weight", |b| {
        b.iter(|| black_box(ppvm_traits_2::Word::weight(black_box(&nw))))
    });

    let ow: OldBare = OldBare::from(LHS);
    g.bench_function("old/weight", |b| {
        b.iter(|| black_box(PauliWordTrait::weight(black_box(&ow))))
    });

    g.finish();
}

criterion_group!(
    benches,
    bench_product,
    bench_clifford_cnot,
    bench_key_hash,
    bench_weight
);
criterion_main!(benches);
