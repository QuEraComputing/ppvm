// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Comparative benchmarks: new `ppvm-phased-pauli-word-2::PhasedPauliWord`
//! (`Phased<PauliWord>`) vs the old `ppvm-pauli-word::PhasedPauliWord` on the hot
//! phased-word paths, so the refactor's performance gate reads off as a new/old
//! ratio per target.
//!
//! Targets (design: `traits-2-implementation-plan.md` Phase 2 perf gate; there is
//! **no** `key_hash` bench — `Phased` is deliberately non-`Indexable`):
//! * the full ℤ₄ phased product (`key_mul` + phase accumulation) — new phased
//!   `MulAssign` vs old phased `MulAssign`.
//! * one Clifford conjugation *with phase tracking* — `cnot`, new vs old (the
//!   phased word tracks the sign delta; the delta the bare word dropped).

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use ppvm_phased_pauli_word_2::PhasedPauliWord as NewPhased;
use ppvm_traits_2::Clifford as NewClifford;

use ppvm_pauli_word::phase::PhasedPauliWord as OldPhasedTy;
use ppvm_traits::traits::Clifford as OldClifford;

type OldPhased = OldPhasedTy<u64>;

const LHS: &str = "+XYZIXYZIXYZIXYZI"; // 16 qubits
const RHS: &str = "+ZZXYIIYXZXYZIXYZ";

fn bench_product(c: &mut Criterion) {
    let mut g = c.benchmark_group("phased_pauli_word/product");

    // New: `MulAssign` (the full ℤ₄ product = key_mul bits + phase accumulation),
    // mutating a persistent word in place so the measurement is the per-op cost.
    let mut na: NewPhased = LHS.into();
    let nb: NewPhased = RHS.into();
    g.bench_function("new/phased_mul_assign", |b| {
        b.iter(|| {
            na *= black_box(&nb);
            black_box(&na);
        })
    });

    // Old: `MulAssign` (old phased word is `Copy`; mutate in place likewise).
    let mut oa: OldPhased = LHS.into();
    let ob: OldPhased = RHS.into();
    g.bench_function("old/phased_mul_assign", |b| {
        b.iter(|| {
            oa *= black_box(&ob);
            black_box(&oa);
        })
    });

    g.finish();
}

fn bench_clifford_cnot(c: &mut Criterion) {
    let mut g = c.benchmark_group("phased_pauli_word/cnot");

    // One phase-tracking `cnot` conjugation per timed iteration, mutating a single
    // persistent word in place — no per-iteration clone/copy to confound it. Both
    // crates track the ℤ₄ sign delta here (this is the phased word's job).
    let mut nw: NewPhased = LHS.into();
    g.bench_function("new/cnot", |b| {
        b.iter(|| nw.cnot(black_box(0), black_box(1)))
    });

    let mut ow: OldPhased = LHS.into();
    g.bench_function("old/cnot", |b| {
        b.iter(|| ow.cnot(black_box(0), black_box(1)))
    });

    g.finish();
}

criterion_group!(benches, bench_product, bench_clifford_cnot);
criterion_main!(benches);
