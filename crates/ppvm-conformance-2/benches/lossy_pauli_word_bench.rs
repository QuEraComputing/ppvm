// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Comparative benchmarks: new `ppvm-lossy-pauli-word-2::LossyPauliWord` vs the
//! old `ppvm-pauli-word::loss::LossyPauliWord` on the hot lossy-word paths, so the
//! refactor's performance gate can be read off as a new/old ratio per target.
//!
//! Targets (design: task brief; `traits-2-implementation-plan.md` Phase 2 perf
//! gate):
//! * **product** — the Pauli twisted product on the present-site projection of a
//!   lossy word (new `KeyProduct::key_mul` vs the old phased `MulAssign`). The
//!   projection to ordinary words is done once in setup (untimed); neither lossy
//!   type carries a native product, so the loss-agnostic product on the present
//!   sites is the operation lossy propagation calls.
//! * **cnot** — one Clifford conjugation on a persistent lossy word. This exposes
//!   the hashing policy: the old lossy word eagerly `rehash()`es all three planes
//!   on every gate, the new one only clears its lazy `OnceLock` components.
//! * **key_hash** — new cold (first compute) vs warm (cached), plus the old
//!   word's warm read, so the design's lazy/`Copy`-drop trade-off is visible.
//! * **weight** / **loss_weight** — the fused popcounts, new vs old.

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use ppvm_lossy_pauli_word_2::LossyPauliWord as NewLossyWord;
use ppvm_pauli_word_2::PauliWord as NewWord;
use ppvm_traits_2::{
    Clifford as NewClifford, IdentityBuildHasher, Indexable, KeyProduct, LossySite, Pauli, Word,
};

use ppvm_pauli_word::loss::LossyPauliWord as OldLossyWord;
use ppvm_pauli_word::phase::PhasedPauliWord;
use ppvm_traits::traits::{Clifford as OldClifford, PauliWordTrait};

use std::hash::BuildHasher;

type NewLossy = NewLossyWord<u64>;
type OldLossy = OldLossyWord<u64>;
type New = NewWord<u64>;
type OldPhased = PhasedPauliWord<u64>;

// 16-qubit lossy words with a mix of present Paulis and loss.
const LHS: &str = "XYZLXYZIXLZIXYZL";
const RHS: &str = "ZLXYIILXZXYZLXYZ";

/// Present-site projection of a *new* lossy word as a plain Pauli string
/// (`Lost ↦ I`).
fn new_projection(w: &NewLossy) -> String {
    (0..w.n_sites())
        .map(|i| match w.get(i) {
            LossySite::Present(Pauli::X) => 'X',
            LossySite::Present(Pauli::Y) => 'Y',
            LossySite::Present(Pauli::Z) => 'Z',
            LossySite::Present(Pauli::I) | LossySite::Lost => 'I',
        })
        .collect()
}

/// Present-site projection of an *old* lossy word (`L ↦ I`).
fn old_projection(w: &OldLossy) -> String {
    use ppvm_traits::char::Pauli as OldPauli;
    (0..w.n_qubits())
        .map(|i| match w.get(i) {
            OldPauli::L => 'I',
            p => p.to_string().chars().next().unwrap(),
        })
        .collect()
}

fn bench_product(c: &mut Criterion) {
    let mut g = c.benchmark_group("lossy_pauli_word/product");

    // Project both lossy operands once (untimed); time only the product kernel.
    let nl: NewLossy = LHS.into();
    let nr: NewLossy = RHS.into();
    let na: New = new_projection(&nl).as_str().into();
    let nb: New = new_projection(&nr).as_str().into();
    g.bench_function("new/key_mul", |b| {
        b.iter(|| black_box(black_box(&na).key_mul(black_box(&nb))))
    });

    let ol: OldLossy = OldLossy::from(LHS);
    let or: OldLossy = OldLossy::from(RHS);
    let oa: OldPhased = format!("+{}", old_projection(&ol)).as_str().into();
    let ob: OldPhased = format!("+{}", old_projection(&or)).as_str().into();
    g.bench_function("old/phased_mul", |b| {
        b.iter(|| black_box(black_box(oa) * black_box(ob)))
    });

    g.finish();
}

fn bench_clifford_cnot(c: &mut Criterion) {
    let mut g = c.benchmark_group("lossy_pauli_word/cnot");

    // One `cnot` per timed iteration, mutating a single persistent word in place.
    // Qubits 0 and 1 are present in `LHS` (`X`,`Y`), so the gate does real work.
    // The old lossy word eagerly `rehash()`es all three planes each gate; the new
    // one only invalidates its lazy `OnceLock` — this is what the ratio exposes.
    let mut nw: NewLossy = LHS.into();
    g.bench_function("new/cnot", |b| {
        b.iter(|| nw.cnot(black_box(0), black_box(1)))
    });

    let mut ow: OldLossy = OldLossy::from(LHS);
    g.bench_function("old/cnot", |b| {
        b.iter(|| ow.cnot(black_box(0), black_box(1)))
    });

    g.finish();
}

fn bench_key_hash(c: &mut Criterion) {
    let mut g = c.benchmark_group("lossy_pauli_word/key_hash");

    // Cold: a never-hashed template, cloned in untimed setup (clone of an empty
    // `OnceLock` stays empty), so the timed closure measures only first-compute.
    let cold_template: NewLossy = LHS.into();
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

    // Warm: hash once up front, then re-read both cached components.
    let warm: NewLossy = LHS.into();
    black_box(warm.key_hash());
    g.bench_function("new/warm", |b| {
        b.iter(|| black_box(black_box(&warm).key_hash()))
    });

    // Old: the hash is computed eagerly at construction and stored; reading it
    // back (what hashbrown does) is a field read. `IdentityBuildHasher` captures
    // exactly the `u64` the old `Hash` impl writes.
    let ow: OldLossy = OldLossy::from(LHS);
    let bh = IdentityBuildHasher;
    g.bench_function("old/warm", |b| {
        b.iter(|| black_box(bh.hash_one(black_box(&ow))))
    });

    g.finish();
}

fn bench_weight(c: &mut Criterion) {
    let mut g = c.benchmark_group("lossy_pauli_word/weight");

    let nw: NewLossy = LHS.into();
    g.bench_function("new/weight", |b| {
        b.iter(|| black_box(Word::weight(black_box(&nw))))
    });

    let ow: OldLossy = OldLossy::from(LHS);
    g.bench_function("old/weight", |b| {
        b.iter(|| black_box(PauliWordTrait::weight(black_box(&ow))))
    });

    g.finish();
}

fn bench_loss_weight(c: &mut Criterion) {
    let mut g = c.benchmark_group("lossy_pauli_word/loss_weight");

    let nw: NewLossy = LHS.into();
    g.bench_function("new/loss_weight", |b| {
        b.iter(|| black_box(black_box(&nw).loss_weight()))
    });

    let ow: OldLossy = OldLossy::from(LHS);
    g.bench_function("old/loss_weight", |b| {
        b.iter(|| black_box(PauliWordTrait::loss_weight(black_box(&ow))))
    });

    g.finish();
}

criterion_group!(
    benches,
    bench_product,
    bench_clifford_cnot,
    bench_key_hash,
    bench_weight,
    bench_loss_weight
);
criterion_main!(benches);
