// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! End-to-end **circuit-propagation** benchmark: the whole Trotter workload
//! (ported from `ppvm-pauli-sum::benches::trotter`) propagated through BOTH the
//! new `ppvm-pauli-sum-2::Sum` and the old `ppvm-pauli-sum::PauliSum`, reporting
//! the same-build new/old TOTAL-circuit wall-clock ratio.
//!
//! **Why this exists (what the single-gate microbench misses).** The other
//! `pauli_sum_bench` targets are single-gate MICRObenches — a tight `b.iter(||
//! new.h(0))` loop over one warm sum. A per-gate fresh allocation looks nearly
//! free there: the allocator hands back the same warm page every iteration, so
//! the cost of *not* reusing the storage double-buffer (the `HashMapStore`
//! `aux`/`scratch`) barely registers. A DEEP circuit is different — each gate
//! sees a differently-shaped support, so a per-gate `HashMap::with_capacity`
//! (instead of the persistent aux ping-pong) pays real allocation churn that
//! compounds over thousands of gates. This bench is the workload that surfaces
//! it.
//!
//! **Fair-comparison note (storage width).** Both sides are pinned to `[u8; 8]`
//! (64-qubit-capacity) storage so the ratio is engine-to-engine, matching the
//! `pauli_sum_bench` module's rationale: the shipped `PauliSum` default is
//! `u64`-backed, but `BitArray<u64>` single-bit ops differ from `[u8; 8]` by a
//! few percent in both directions, which would fold a storage-codegen delta into
//! the engine ratio. Correctness on the shipped `u64` default is covered by the
//! differential suite.
//!
//! Params mirror the old `trotter` bench: `n = 12`, `h = 1`, `dt = 0.1`,
//! `time = 1.0` (→ 10 Trotter steps), `J = 1/8`, `CombinedStrategy(
//! CoefficientThreshold(1e-6), MaxPauliWeight(usize::MAX))` so allocation churn
//! over many gates is real (weight-cap disabled, coefficient floor only — the old
//! bench's config).

use criterion::{Criterion, criterion_group, criterion_main};

// --- Old crate: CombinedStrategy config + gate traits. -------------------------
use ppvm_pauli_sum::config::fxhash::ByteF64 as OldByteF64;
use ppvm_pauli_sum::strategy::{
    CoefficientThreshold as OldCoeffThreshold, CombinedStrategy, MaxPauliWeight as OldMaxWeight,
};
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_traits::traits::{
    Clifford as OldClifford, PauliError as OldPauliError, RotationOne as OldRotationOne,
};

// --- New crate: storage-matched `[u8; 8]` sum + CombinedPolicy + gate traits. --
use ppvm_pauli_sum_2::{
    CoefficientThreshold as NewCoeffThreshold, CombinedPolicy, HashMapStore,
    MaxPauliWeight as NewMaxWeight, PauliWord as NewPauliWord, Sum,
};
use ppvm_traits_2::{
    Clifford as NewClifford, PauliError as NewPauliError, RotationOne as NewRotationOne,
};

/// Shared truncation floor (both crates).
const THRESHOLD: f64 = 1e-6;

// Old side: `[u8; 8]` + CombinedStrategy(1e-6, MaxPauliWeight).
type OldStrat = CombinedStrategy<OldCoeffThreshold, OldMaxWeight>;
type OldCfg = OldByteF64<8, OldStrat>;
type OldSum = OldPauliSum<OldCfg>;

// New side: storage-matched `[u8; 8]` key + CombinedPolicy(1e-6, MaxPauliWeight).
type NewKey = NewPauliWord<[u8; 8]>;
type NewPolicy = CombinedPolicy<NewCoeffThreshold, NewMaxWeight>;
type NewSum = Sum<HashMapStore<NewKey, f64>, NewPolicy>;

fn old_strat() -> OldStrat {
    CombinedStrategy(OldCoeffThreshold(THRESHOLD), OldMaxWeight(usize::MAX))
}

fn new_policy() -> NewPolicy {
    CombinedPolicy(
        NewCoeffThreshold {
            threshold: THRESHOLD,
        },
        NewMaxWeight(usize::MAX),
    )
}

/// The initial observable `Σ_i Z_i` as `(pauli_string, coeff)` terms.
fn sum_z_terms(n: usize) -> Vec<(String, f64)> {
    (0..n)
        .map(|i| {
            let s: String = (0..n).map(|j| if j == i { 'Z' } else { 'I' }).collect();
            (s, 1.0)
        })
        .collect()
}

fn build_old(n: usize) -> OldSum {
    let mut s: OldSum = OldPauliSum::builder()
        .n_qubits(n)
        .strategy(old_strat())
        .capacity(n.pow(2))
        .build();
    for (word, c) in sum_z_terms(n) {
        s += (word.as_str(), c);
    }
    s
}

fn build_new(n: usize) -> NewSum {
    NewSum::from_terms_with_policy(
        n,
        new_policy(),
        sum_z_terms(n)
            .into_iter()
            .map(|(word, c)| (NewKey::from(word.as_str()), c)),
    )
}

/// One `rzz(a, b, θ)` decomposed as `cnot(a, b); rz(b, θ); cnot(a, b)` — the same
/// decomposition the differential test validates, applied identically on both
/// sides so the benchmarked circuits are the same.
#[inline]
fn old_rzz(state: &mut OldSum, a: usize, b: usize, theta: f64) {
    state.cnot(a, b);
    state.rz(b, theta);
    state.cnot(a, b);
}

#[inline]
fn new_rzz(state: &mut NewSum, a: usize, b: usize, theta: f64) {
    state.cnot(a, b);
    state.rz(b, theta);
    state.cnot(a, b);
}

fn trotter_old(
    state: &mut OldSum,
    n: usize,
    steps: usize,
    theta_x: f64,
    theta_zz: f64,
    noise: [f64; 3],
) {
    for _ in 0..steps {
        for i in 0..n {
            state.pauli_error(i, noise);
            state.truncate();
            state.rx(i, theta_x);
            state.truncate();
        }
        for i in 0..n - 1 {
            state.pauli_error(i + 1, noise);
            state.truncate();
            state.pauli_error(i, noise);
            state.truncate();
            old_rzz(state, i, i + 1, theta_zz);
            state.truncate();
        }
    }
}

fn trotter_new(
    state: &mut NewSum,
    n: usize,
    steps: usize,
    theta_x: f64,
    theta_zz: f64,
    noise: [f64; 3],
) {
    for _ in 0..steps {
        for i in 0..n {
            state.pauli_error(i, noise);
            state.truncate();
            state.rx(i, theta_x);
            state.truncate();
        }
        for i in 0..n - 1 {
            state.pauli_error(i + 1, noise);
            state.truncate();
            state.pauli_error(i, noise);
            state.truncate();
            new_rzz(state, i, i + 1, theta_zz);
            state.truncate();
        }
    }
}

fn bench_trotter(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_sum/integration_trotter");

    // Fuller params (the old `trotter` bench): n = 12, 10 Trotter steps.
    let n = 12usize;
    let h = 1.0_f64;
    let dt = 0.1 / h;
    let time = 1.0 / h;
    let j = 1.0 / 8.0 * h;
    let steps = (time / dt) as usize;
    let theta_x = dt * h;
    let theta_zz = dt * j;
    let noise = [1e-4 / 4.0; 3];

    // Both seeds are the `Σ Z_i` observable, built once; each timed iteration
    // clones the seed and propagates the whole circuit (`iter_batched_ref`,
    // matching the old trotter bench's `SmallInput` clone-per-iter).
    let new_seed = build_new(n);
    g.bench_function("new/trotter", |b| {
        b.iter_batched_ref(
            || new_seed.clone(),
            |state| trotter_new(state, n, steps, theta_x, theta_zz, noise),
            criterion::BatchSize::SmallInput,
        )
    });

    let old_seed = build_old(n);
    g.bench_function("old/trotter", |b| {
        b.iter_batched_ref(
            || old_seed.clone(),
            |state| trotter_old(state, n, steps, theta_x, theta_zz, noise),
            criterion::BatchSize::SmallInput,
        )
    });

    g.finish();
}

criterion_group!(benches, bench_trotter);
criterion_main!(benches);
