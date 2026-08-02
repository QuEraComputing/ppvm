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
    // Same capacity override as `build_old` (`n²`) and the same accumulating
    // `+=` seeding path, so the ratio is engine-to-engine rather than an
    // artifact of the two maps being sized differently.
    let mut s = NewSum::with_capacity(n, new_policy(), n.pow(2));
    for (word, c) in sum_z_terms(n) {
        s += (NewKey::from(word.as_str()), c);
    }
    s
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

/// The Trotter circuit parameters shared by every benchmark in this file
/// (the old `trotter` bench's): `n = 12`, `h = 1`, `dt = 0.1`, `time = 1.0`
/// (→ 10 steps), `J = 1/8`.
struct Params {
    n: usize,
    steps: usize,
    theta_x: f64,
    theta_zz: f64,
    noise: [f64; 3],
}

fn params() -> Params {
    let h = 1.0_f64;
    let dt = 0.1 / h;
    let time = 1.0 / h;
    let j = 1.0 / 8.0 * h;
    Params {
        n: 12,
        steps: (time / dt) as usize,
        theta_x: dt * h,
        theta_zz: dt * j,
        noise: [1e-4 / 4.0; 3],
    }
}

// ---------------------------------------------------------------------------
// 2. Per-op benchmarks on a REALISTIC support (`cnot` re-key and `truncate`).
// ---------------------------------------------------------------------------
//
// These are the two coverage holes the end-to-end regression exposed: `cnot`
// (the `RekeyBijective` path) and `truncate` (the policy's retain scans, run
// after every gate in the Trotter loop) had **no** benchmark at all, so a
// regression in either was invisible.
//
// They are *not* the usual tight one-gate microbench. The two problems with
// that shape are (a) the support is a toy, so per-term costs (hashing, probe
// locality) are unrepresentative, and (b) the allocator recycles one warm page.
// Here the state is a support grown by actually propagating the Trotter circuit,
// and each timed iteration does a full sweep over it.
//
// Neither benchmark clones per iteration — both timed bodies are **state
// preserving**, so `b.iter` can hammer one realistic support directly and the
// measurement is the operation, not a `HashMap` clone:
//   * `cnot(a, b)` is an involution: applying it twice restores the support,
//     so a forward+backward pair per bond leaves the state exactly as it was.
//   * `truncate()` is idempotent: after the first call nothing more is
//     droppable, which is precisely the Trotter loop's steady state (truncate
//     runs after every gate and usually removes little).

/// Grow a realistic support by propagating `steps` Trotter steps.
fn grown_new(p: &Params, steps: usize) -> NewSum {
    let mut s = build_new(p.n);
    trotter_new(&mut s, p.n, steps, p.theta_x, p.theta_zz, p.noise);
    s
}

fn grown_old(p: &Params, steps: usize) -> OldSum {
    let mut s = build_old(p.n);
    trotter_old(&mut s, p.n, steps, p.theta_x, p.theta_zz, p.noise);
    s
}

fn bench_rekey_and_truncate(c: &mut Criterion) {
    let p = params();
    // Grow with the FULL circuit, so the benchmarked support is the one the deep
    // circuit actually ends up propagating rather than an early-and-small one.
    let mut new_state = grown_new(&p, p.steps);
    let mut old_state = grown_old(&p, p.steps);
    assert_eq!(
        new_state.len(),
        old_state.data().len(),
        "the two engines must be benchmarked on the same-size support"
    );
    // Printed so the reported per-op times can be read per term, and so a future
    // change to the growth circuit that silently shrinks the support is visible.
    println!(
        "[pauli_sum/rekey_cnot + truncate] benchmarked support: {} terms",
        new_state.len()
    );

    // --- `cnot` re-key: one forward+backward sweep over every bond. ----------
    let mut g = c.benchmark_group("pauli_sum/rekey_cnot");
    let before = new_state.len();
    g.bench_function("new/cnot_sweep", |b| {
        b.iter(|| {
            for i in 0..p.n - 1 {
                new_state.cnot(i, i + 1);
                new_state.cnot(i, i + 1);
            }
        })
    });
    assert_eq!(
        new_state.len(),
        before,
        "cnot sweep must be state-preserving (cnot is an involution)"
    );

    let before = old_state.data().len();
    g.bench_function("old/cnot_sweep", |b| {
        b.iter(|| {
            for i in 0..p.n - 1 {
                old_state.cnot(i, i + 1);
                old_state.cnot(i, i + 1);
            }
        })
    });
    assert_eq!(
        old_state.data().len(),
        before,
        "cnot sweep must be state-preserving (cnot is an involution)"
    );
    g.finish();

    // --- `truncate`: the idempotent steady-state retain scan. ----------------
    let mut g = c.benchmark_group("pauli_sum/truncate");
    new_state.truncate();
    old_state.truncate();
    g.bench_function("new/truncate", |b| b.iter(|| new_state.truncate()));
    g.bench_function("old/truncate", |b| b.iter(|| old_state.truncate()));
    g.finish();

    // --- The remaining per-op classes: `pauli_error` (diagonal channel) and
    //     `rx` (branching rotation), each a full sweep over all n qubits.
    //
    // Unlike `cnot`/`truncate` these are NOT state-preserving — the channel
    // contracts every coefficient and the rotation mixes them — so a plain
    // `b.iter` would drive the support toward denormals over millions of
    // iterations and silently corrupt the measurement. `iter_batched_ref`
    // restores a clean clone per iteration, and criterion does NOT time the
    // setup closure, so the clone does not enter the reported number.
    let noise = p.noise;
    let theta = p.theta_x;

    let mut g = c.benchmark_group("pauli_sum/pauli_error");
    g.bench_function("new/pauli_error_sweep", |b| {
        b.iter_batched_ref(
            || new_state.clone(),
            |s| {
                for i in 0..p.n {
                    s.pauli_error(i, noise);
                }
            },
            criterion::BatchSize::LargeInput,
        )
    });
    g.bench_function("old/pauli_error_sweep", |b| {
        b.iter_batched_ref(
            || old_state.clone(),
            |s| {
                for i in 0..p.n {
                    s.pauli_error(i, noise);
                }
            },
            criterion::BatchSize::LargeInput,
        )
    });
    g.finish();

    let mut g = c.benchmark_group("pauli_sum/rx");
    g.bench_function("new/rx_sweep", |b| {
        b.iter_batched_ref(
            || new_state.clone(),
            |s| {
                for i in 0..p.n {
                    s.rx(i, theta);
                }
            },
            criterion::BatchSize::LargeInput,
        )
    });
    g.bench_function("old/rx_sweep", |b| {
        b.iter_batched_ref(
            || old_state.clone(),
            |s| {
                for i in 0..p.n {
                    s.rx(i, theta);
                }
            },
            criterion::BatchSize::LargeInput,
        )
    });
    g.finish();
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

criterion_group!(benches, bench_trotter, bench_rekey_and_truncate);
criterion_main!(benches);
