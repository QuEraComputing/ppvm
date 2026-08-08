// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The integration-baseline workloads that `pauli_sum_integration` does not
//! cover, each as a same-build new-vs-old **total wall-clock** comparison:
//!
//! * **workload 2 — qubit-scaling sweep**: the noisy-TFIM Trotter circuit under
//!   `CoefficientThreshold(1e-6)` alone, swept over `n`. Reported as a ratio
//!   *curve*: this is the regime that exposes hash quality / bucket clustering at
//!   high fill, which a single width cannot show.
//! * **workload 3 — untruncated deep circuit**: rotation layers + a `cnot` ring
//!   with **no** truncation, so the support grows by pure fan-out. The opposite
//!   regime from the truncation-bounded Trotter loop: it stresses map resizing,
//!   capacity hints and the aux ping-pong under monotone growth.
//! * **workload 4 — truncation cost grid**: `truncate()` over
//!   `(weight profile) × (policy cell)`, including the `usize::MAX` disable
//!   sentinel, which must be ~free on both sides.
//!
//! # Separate benchmark target
//!
//! These longer workload sweeps are kept out of `pauli_sum_integration` so that
//! target remains focused and each benchmark process has a manageable runtime.
//!
//! # Fair-config note
//!
//! Storage is pinned on both sides — `[u8; 8]` for the circuit workloads (64-qubit
//! capacity), `[u8; 16]` for the 128-qubit truncation grid — with `f64`
//! coefficients and the algebraically identical policy/strategy pair, so every
//! ratio is engine-to-engine rather than a storage-codegen artifact.

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use std::collections::BTreeMap;

use ppvm_pauli_sum::config::fxhash::ByteF64 as OldByteF64;
use ppvm_pauli_sum::strategy::{
    CoefficientThreshold as OldCoeffThreshold, CombinedStrategy, MaxPauliWeight as OldMaxWeight,
};
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_traits::traits::{
    Clifford as OldClifford, PauliError as OldPauliError, RotationOne as OldRotationOne,
};

use ppvm_pauli_sum_2::{
    CoefficientThreshold as NewCoeffThreshold, CombinedPolicy, HashMapStore,
    MaxPauliWeight as NewMaxWeight, NoPolicy, PauliWord as NewPauliWord, Sum,
};
use ppvm_traits_2::{
    Clifford as NewClifford, PauliError as NewPauliError, RotationOne as NewRotationOne,
};

const THRESHOLD: f64 = 1e-6;

type NewKey = NewPauliWord<[u8; 8]>;
type NewWideKey = NewPauliWord<[u8; 16]>;

type OldThreshSum = OldPauliSum<OldByteF64<8, OldCoeffThreshold>>;
type NewThreshSum = Sum<HashMapStore<NewKey, f64>, NewCoeffThreshold>;

type OldExactSum = OldPauliSum<OldByteF64<8>>;
type NewExactSum = Sum<HashMapStore<NewKey, f64>, NoPolicy>;

macro_rules! old_support {
    ($sum:expr) => {
        $sum.data()
            .iter()
            .map(|(key, coeff)| (key.to_string(), *coeff))
            .collect::<BTreeMap<_, _>>()
    };
}

macro_rules! new_support {
    ($sum:expr) => {
        $sum.iter()
            .map(|(key, coeff)| (key.to_string(), coeff))
            .collect::<BTreeMap<_, _>>()
    };
}

#[track_caller]
fn assert_supports_match(
    old: BTreeMap<String, f64>,
    new: BTreeMap<String, f64>,
    tolerance: f64,
    label: &str,
) {
    assert_eq!(
        old.keys().collect::<Vec<_>>(),
        new.keys().collect::<Vec<_>>(),
        "[{label}] canonical support keys differ"
    );
    for (key, old_coeff) in old {
        let new_coeff = new[&key];
        let bound = tolerance.max(old_coeff.abs() * tolerance);
        assert!(
            (old_coeff - new_coeff).abs() <= bound,
            "[{label}] coefficient at {key} differs: old={old_coeff}, new={new_coeff}, tol={bound}"
        );
    }
}

#[track_caller]
fn assert_threshold_supports_match(
    old: BTreeMap<String, f64>,
    new: BTreeMap<String, f64>,
    label: &str,
) {
    // Merge order may put a coefficient a few ulps to either side of the
    // inclusive threshold. Keep the integration suite's intentional 1% band.
    const BOUNDARY_BAND: f64 = THRESHOLD * 1.01;
    let solid = |support: BTreeMap<String, f64>| {
        support
            .into_iter()
            .filter(|(_, coeff)| coeff.abs() >= BOUNDARY_BAND)
            .collect()
    };
    assert_supports_match(solid(old), solid(new), 1e-9, label);
}

/// `Σ_i Z_i` as `(string, coeff)` terms.
fn sum_z_terms(n: usize) -> Vec<(String, f64)> {
    (0..n)
        .map(|i| {
            let s: String = (0..n).map(|j| if j == i { 'Z' } else { 'I' }).collect();
            (s, 1.0)
        })
        .collect()
}

/// The Trotter body using the same explicit `cnot; rz; cnot` decomposition on
/// both sides, with caller-driven truncation after every single operation.
macro_rules! trotter_evolve {
    ($state:expr, $n:expr, $steps:expr, $theta_x:expr, $theta_zz:expr, $noise:expr) => {{
        for _ in 0..$steps {
            for i in 0..$n {
                $state.pauli_error(i, $noise);
                $state.truncate();
                $state.rx(i, $theta_x);
                $state.truncate();
            }
            for i in 0..$n - 1 {
                $state.pauli_error(i + 1, $noise);
                $state.truncate();
                $state.pauli_error(i, $noise);
                $state.truncate();
                $state.cnot(i, i + 1);
                $state.rz(i + 1, $theta_zz);
                $state.cnot(i, i + 1);
                $state.truncate();
            }
        }
    }};
}

// ---------------------------------------------------------------------------
// Workload 2 — qubit-scaling sweep.
// ---------------------------------------------------------------------------

/// Sweep widths. Old's example sweeps `(2..65).step_by(4)`; the top of that range
/// at `J = 1.0` runs minutes per sample, so this takes the same shape over the
/// widths that still fit a criterion sample budget while reaching a ~10³–10⁴-term
/// support at the wide end.
const SWEEP_WIDTHS: [usize; 6] = [4, 12, 20, 28, 36, 44];
/// Trotter steps per sweep point.
const SWEEP_STEPS: usize = 10;

fn bench_qubit_sweep(c: &mut Criterion) {
    let dt = 0.1_f64;
    let theta_x = dt;
    // `J = 1.0` (not the headline `1/8`): drives the support large, which is the
    // point of the sweep — high map fill is where bucket clustering shows.
    let theta_zz = dt;
    let noise = [1e-4 / 4.0; 3];

    let mut g = c.benchmark_group("pauli_sum/workload_qubit_sweep");
    for n in SWEEP_WIDTHS {
        let mut old_seed: OldThreshSum = OldPauliSum::builder()
            .n_qubits(n)
            .strategy(OldCoeffThreshold(THRESHOLD))
            .capacity(n.pow(2))
            .build();
        for (w, v) in sum_z_terms(n) {
            old_seed += (w.as_str(), v);
        }
        let mut new_seed = NewThreshSum::with_capacity(
            n,
            NewCoeffThreshold {
                threshold: THRESHOLD,
            },
            n.pow(2),
        );
        for (w, v) in sum_z_terms(n) {
            new_seed += (NewKey::from(w.as_str()), v);
        }

        // Validate canonical keys and coefficients outside the timed path, then
        // report support so workload-size movement remains visible in the log.
        {
            let mut new_probe = new_seed.clone();
            let mut old_probe = old_seed.clone();
            trotter_evolve!(new_probe, n, SWEEP_STEPS, theta_x, theta_zz, noise);
            trotter_evolve!(old_probe, n, SWEEP_STEPS, theta_x, theta_zz, noise);
            assert_threshold_supports_match(
                old_support!(old_probe),
                new_support!(new_probe),
                &format!("qubit sweep n={n}"),
            );
            println!(
                "[workload_qubit_sweep] n={n}: final support {} terms",
                new_probe.len()
            );
        }

        g.bench_function(format!("new/n{n}"), |b| {
            b.iter_batched_ref(
                || new_seed.clone(),
                |s| trotter_evolve!(s, n, SWEEP_STEPS, theta_x, theta_zz, noise),
                BatchSize::SmallInput,
            )
        });
        g.bench_function(format!("old/n{n}"), |b| {
            b.iter_batched_ref(
                || old_seed.clone(),
                |s| trotter_evolve!(s, n, SWEEP_STEPS, theta_x, theta_zz, noise),
                BatchSize::SmallInput,
            )
        });
    }
    g.finish();
}

// ---------------------------------------------------------------------------
// Attribution — the decomposed Trotter control with the Clifford re-key ABLATED.
// ---------------------------------------------------------------------------
//
// `full` uses the explicit `cnot; rz; cnot` decomposition (`n = 12`, 10 steps,
// `Combined(1e-6, MAX)`, capacity `n²`); `no_rekey` removes the two `cnot`s from
// each decomposition while retaining the other operations, truncate placement,
// policy, seed and storage. This is an attribution control for the
// `RekeyBijective` path, not the native-`rzz` integration workload.
//
// (The ablated circuit is not physically the same evolution — that is fine and
// deliberate: it is a timing control, never a numeric one.)

/// The Trotter body with the `rzz` decomposition's `cnot` pair removed.
macro_rules! trotter_evolve_no_rekey {
    ($state:expr, $n:expr, $steps:expr, $theta_x:expr, $theta_zz:expr, $noise:expr) => {{
        for _ in 0..$steps {
            for i in 0..$n {
                $state.pauli_error(i, $noise);
                $state.truncate();
                $state.rx(i, $theta_x);
                $state.truncate();
            }
            for i in 0..$n - 1 {
                $state.pauli_error(i + 1, $noise);
                $state.truncate();
                $state.pauli_error(i, $noise);
                $state.truncate();
                $state.rz(i + 1, $theta_zz);
                $state.truncate();
            }
        }
    }};
}

fn bench_trotter_ablation(c: &mut Criterion) {
    let n = 12usize;
    let dt = 0.1_f64;
    let steps = 10usize;
    let theta_x = dt;
    let theta_zz = dt / 8.0;
    let noise = [1e-4 / 4.0; 3];

    let combined_old = CombinedStrategy(OldCoeffThreshold(THRESHOLD), OldMaxWeight(usize::MAX));
    let combined_new = CombinedPolicy(
        NewCoeffThreshold {
            threshold: THRESHOLD,
        },
        NewMaxWeight(usize::MAX),
    );

    let mut old_seed: OldPauliSum<
        OldByteF64<8, CombinedStrategy<OldCoeffThreshold, OldMaxWeight>>,
    > = OldPauliSum::builder()
        .n_qubits(n)
        .strategy(combined_old)
        .capacity(n.pow(2))
        .build();
    for (w, v) in sum_z_terms(n) {
        old_seed += (w.as_str(), v);
    }
    let mut new_seed: Sum<
        HashMapStore<NewKey, f64>,
        CombinedPolicy<NewCoeffThreshold, NewMaxWeight>,
    > = Sum::with_capacity(n, combined_new, n.pow(2));
    for (w, v) in sum_z_terms(n) {
        new_seed += (NewKey::from(w.as_str()), v);
    }

    for (label, no_rekey) in [("full", false), ("no_rekey", true)] {
        let mut old_probe = old_seed.clone();
        let mut new_probe = new_seed.clone();
        if no_rekey {
            trotter_evolve_no_rekey!(old_probe, n, steps, theta_x, theta_zz, noise);
            trotter_evolve_no_rekey!(new_probe, n, steps, theta_x, theta_zz, noise);
        } else {
            trotter_evolve!(old_probe, n, steps, theta_x, theta_zz, noise);
            trotter_evolve!(new_probe, n, steps, theta_x, theta_zz, noise);
        }
        assert_threshold_supports_match(
            old_support!(old_probe),
            new_support!(new_probe),
            &format!("trotter ablation {label}"),
        );
    }

    let mut g = c.benchmark_group("pauli_sum/workload_trotter_ablation");
    g.bench_function("new/full", |b| {
        b.iter_batched_ref(
            || new_seed.clone(),
            |s| trotter_evolve!(s, n, steps, theta_x, theta_zz, noise),
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/full", |b| {
        b.iter_batched_ref(
            || old_seed.clone(),
            |s| trotter_evolve!(s, n, steps, theta_x, theta_zz, noise),
            BatchSize::SmallInput,
        )
    });
    g.bench_function("new/no_rekey", |b| {
        b.iter_batched_ref(
            || new_seed.clone(),
            |s| trotter_evolve_no_rekey!(s, n, steps, theta_x, theta_zz, noise),
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/no_rekey", |b| {
        b.iter_batched_ref(
            || old_seed.clone(),
            |s| trotter_evolve_no_rekey!(s, n, steps, theta_x, theta_zz, noise),
            BatchSize::SmallInput,
        )
    });
    g.finish();
}

// ---------------------------------------------------------------------------
// Workload 3 — untruncated deep random circuit (pure fan-out growth).
// ---------------------------------------------------------------------------

/// `depth` × (rotation layer + `cnot` ring), closed by a final rotation layer —
/// old's `benches/random-circuit.rs` shape. No truncation anywhere.
macro_rules! random_circuit_evolve {
    ($state:expr, $n:expr, $depth:expr) => {{
        for _ in 0..$depth {
            for i in 0..$n {
                $state.rz(i, 1.1);
                $state.ry(i, 2.1);
                $state.rz(i, 1.1);
            }
            for i in 0..$n {
                $state.cnot(i, (i + 1) % $n);
            }
        }
        for i in 0..$n {
            $state.rz(i, 1.1);
            $state.ry(i, 2.1);
            $state.rz(i, 1.1);
        }
    }};
}

fn bench_untruncated_circuit(c: &mut Criterion) {
    // `n = 8, depth = 2` reaches a large fraction of the 4⁸ = 65 536 key space
    // without making the benchmark impractically large. This exercises map
    // growth and the auxiliary-store ping-pong without truncation.
    let n = 8usize;
    let depth = 2usize;
    let capacity = n.pow(2);
    let zz: String = (0..n).map(|i| if i < 2 { 'Z' } else { 'I' }).collect();

    let mut old_seed: OldExactSum = OldPauliSum::builder()
        .n_qubits(n)
        .capacity(capacity)
        .build();
    old_seed += (zz.as_str(), 1.0);
    let mut new_seed = NewExactSum::with_capacity(n, NoPolicy, capacity);
    new_seed += (NewKey::from(zz.as_str()), 1.0);

    {
        let mut new_probe = new_seed.clone();
        let mut old_probe = old_seed.clone();
        random_circuit_evolve!(new_probe, n, depth);
        random_circuit_evolve!(old_probe, n, depth);
        assert_supports_match(
            old_support!(old_probe),
            new_support!(new_probe),
            1e-9,
            "untruncated circuit",
        );
        println!(
            "[workload_random_circuit] n={n} depth={depth}: final support {} terms",
            new_probe.len()
        );
    }

    let mut g = c.benchmark_group("pauli_sum/workload_random_circuit");
    g.bench_function("new/circuit", |b| {
        b.iter_batched_ref(
            || new_seed.clone(),
            |s| random_circuit_evolve!(s, n, depth),
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/circuit", |b| {
        b.iter_batched_ref(
            || old_seed.clone(),
            |s| random_circuit_evolve!(s, n, depth),
            BatchSize::SmallInput,
        )
    });
    g.finish();
}

// ---------------------------------------------------------------------------
// Workload 4 — truncation cost grid (n = 128, 1000 terms).
// ---------------------------------------------------------------------------

const TRUNC_N: usize = 128;
const TRUNC_TERMS: usize = 1000;

/// The bench's term list for one weight profile (same construction as
/// `tests/pauli_sum_truncation_boundary_diff.rs`, whose cells this times).
fn profile_terms(target_weight: usize) -> Vec<(String, f64)> {
    let stride = (TRUNC_N / target_weight).max(1);
    (0..TRUNC_TERMS)
        .map(|k| {
            let mut w = vec!['I'; TRUNC_N];
            for j in 0..target_weight {
                let pos = (j * stride + k) % TRUNC_N;
                w[pos] = ['X', 'Y', 'Z'][(k + j) % 3];
            }
            let extra = (k * 7 + 3) % TRUNC_N;
            if w[extra] == 'I' {
                w[extra] = ['X', 'Y', 'Z'][k % 3];
            }
            (w.into_iter().collect::<String>(), 1.0 / (k as f64 + 1.0))
        })
        .collect()
}

/// Time one `(profile, policy)` cell on both engines. `truncate()` is timed from
/// a fresh clone each iteration (criterion does not time the setup closure), so
/// the measurement is the *dropping* pass, not the idempotent re-scan.
macro_rules! truncate_cell {
    ($group:expr, $label:expr, $terms:expr, $old_cfg:ty, $old_strat:expr, $new_pol:ty, $new_pol_val:expr) => {{
        let mut old: OldPauliSum<$old_cfg> = OldPauliSum::builder()
            .n_qubits(TRUNC_N)
            .strategy($old_strat)
            .capacity(TRUNC_TERMS * 2)
            .build();
        for (w, v) in $terms {
            old += (w.as_str(), *v);
        }
        let mut new: Sum<HashMapStore<NewWideKey, f64>, $new_pol> =
            Sum::with_capacity(TRUNC_N, $new_pol_val, TRUNC_TERMS * 2);
        for (w, v) in $terms {
            new += (NewWideKey::from(w.as_str()), *v);
        }

        let mut old_probe = old.clone();
        let mut new_probe = new.clone();
        old_probe.truncate();
        new_probe.truncate();
        assert_supports_match(
            old_support!(old_probe),
            new_support!(new_probe),
            0.0,
            &$label,
        );

        $group.bench_function(format!("new/{}", $label), |b| {
            b.iter_batched_ref(|| new.clone(), |s| s.truncate(), BatchSize::LargeInput)
        });
        $group.bench_function(format!("old/{}", $label), |b| {
            b.iter_batched_ref(|| old.clone(), |s| s.truncate(), BatchSize::LargeInput)
        });
    }};
}

fn bench_truncation_grid(c: &mut Criterion) {
    let mut g = c.benchmark_group("pauli_sum/workload_truncate");
    for target in [3usize, 50, 120] {
        let terms = profile_terms(target);
        for w in [10usize, 1000, usize::MAX] {
            let label = if w == usize::MAX {
                format!("w{target}/max_sentinel")
            } else {
                format!("w{target}/cut{w}")
            };
            truncate_cell!(
                g,
                label,
                &terms,
                OldByteF64<16, OldMaxWeight>,
                OldMaxWeight(w),
                NewMaxWeight,
                NewMaxWeight(w)
            );
        }
        truncate_cell!(
            g,
            format!("w{target}/threshold"),
            &terms,
            OldByteF64<16, OldCoeffThreshold>,
            OldCoeffThreshold(1e-12),
            NewCoeffThreshold,
            NewCoeffThreshold { threshold: 1e-12 }
        );
        truncate_cell!(
            g,
            format!("w{target}/combined"),
            &terms,
            OldByteF64<16, CombinedStrategy<OldCoeffThreshold, OldMaxWeight>>,
            CombinedStrategy(OldCoeffThreshold(1e-12), OldMaxWeight(10)),
            CombinedPolicy<NewCoeffThreshold, NewMaxWeight>,
            CombinedPolicy(NewCoeffThreshold { threshold: 1e-12 }, NewMaxWeight(10))
        );
    }
    g.finish();
}

criterion_group!(
    benches,
    bench_qubit_sweep,
    bench_trotter_ablation,
    bench_untruncated_circuit,
    bench_truncation_grid
);
criterion_main!(benches);
