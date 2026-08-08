// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Per-op-class time-attribution harness for the deep Trotter workload, new vs
//! old, SAME build (both engines in one binary — the sound instrument; cross-build
//! absolutes swing from code-alignment/Mytkowicz bias). Wraps each op class
//! (pauli_error / rx / native rzz / the *explicit* truncate) with `Instant`
//! accumulators so the end-to-end new-vs-old gap can be split by operation. This is
//! the tool that decomposed the `ps2.trotter.perf` regression; kept for the
//! continuation of that investigation (see docs/log.md → "▶ Continue here").
//!
//! Note: the per-op timings *include* each gate's INTERNAL auto-truncate (a
//! user-facing-behaviour gap being fixed — see the log). To isolate the internal
//! truncate, A/B by removing `self.policy.truncate(&mut self.storage)` from the
//! `Sum` gate drivers (`apply`/`rekey_bijective`/`rotate_in_place`) and re-run.
//! Run:  cargo run --release -p ppvm-conformance-2 --example trotter_attrib

use std::time::{Duration, Instant};

use ppvm_pauli_sum::config::fxhash::ByteF64 as OldByteF64;
use ppvm_pauli_sum::strategy::{
    CoefficientThreshold as OldCoeffThreshold, CombinedStrategy, MaxPauliWeight as OldMaxWeight,
};
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_traits::traits::{
    PauliError as OldPauliError, RotationOne as OldRotationOne, RotationTwo as OldRotationTwo,
};

use ppvm_pauli_sum_2::{
    CoefficientThreshold as NewCoeffThreshold, CombinedPolicy, HashMapStore,
    MaxPauliWeight as NewMaxWeight, PauliWord as NewPauliWord, Sum,
};
use ppvm_traits_2::{
    PauliError as NewPauliError, RotationOne as NewRotationOne, RotationTwo as NewRotationTwo,
};

const THRESHOLD: f64 = 1e-6;

type OldStrat = CombinedStrategy<OldCoeffThreshold, OldMaxWeight>;
type OldSum = OldPauliSum<OldByteF64<8, OldStrat>>;
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
fn sum_z_terms(n: usize) -> Vec<(String, f64)> {
    (0..n)
        .map(|i| {
            (
                (0..n)
                    .map(|j| if j == i { 'Z' } else { 'I' })
                    .collect::<String>(),
                1.0,
            )
        })
        .collect()
}
fn build_old(n: usize) -> OldSum {
    let mut s: OldSum = OldPauliSum::builder()
        .n_qubits(n)
        .strategy(old_strat())
        .capacity(n.pow(2))
        .build();
    for (w, c) in sum_z_terms(n) {
        s += (w.as_str(), c);
    }
    s
}
fn build_new(n: usize) -> NewSum {
    NewSum::from_terms_with_policy(
        n,
        new_policy(),
        sum_z_terms(n)
            .into_iter()
            .map(|(w, c)| (NewKey::from(w.as_str()), c)),
    )
}

#[derive(Default, Clone, Copy)]
struct T {
    pauli_error: Duration,
    rx: Duration,
    rz: Duration,
    cnot: Duration,
    truncate: Duration,
    total: Duration,
}
macro_rules! timed {
    ($acc:expr, $body:expr) => {{
        let s = Instant::now();
        $body;
        $acc += s.elapsed();
    }};
}

fn trotter_old(
    state: &mut OldSum,
    n: usize,
    steps: usize,
    tx: f64,
    tzz: f64,
    noise: [f64; 3],
    t: &mut T,
) {
    for _ in 0..steps {
        for i in 0..n {
            timed!(t.pauli_error, state.pauli_error(i, noise));
            timed!(t.truncate, state.truncate());
            timed!(t.rx, state.rx(i, tx));
            timed!(t.truncate, state.truncate());
        }
        for i in 0..n - 1 {
            timed!(t.pauli_error, state.pauli_error(i + 1, noise));
            timed!(t.truncate, state.truncate());
            timed!(t.pauli_error, state.pauli_error(i, noise));
            timed!(t.truncate, state.truncate());
            timed!(t.rz, state.rzz(i, i + 1, tzz));
            timed!(t.truncate, state.truncate());
        }
    }
}
fn trotter_new(
    state: &mut NewSum,
    n: usize,
    steps: usize,
    tx: f64,
    tzz: f64,
    noise: [f64; 3],
    t: &mut T,
) {
    for _ in 0..steps {
        for i in 0..n {
            timed!(t.pauli_error, state.pauli_error(i, noise));
            timed!(t.truncate, state.truncate());
            timed!(t.rx, state.rx(i, tx));
            timed!(t.truncate, state.truncate());
        }
        for i in 0..n - 1 {
            timed!(t.pauli_error, state.pauli_error(i + 1, noise));
            timed!(t.truncate, state.truncate());
            timed!(t.pauli_error, state.pauli_error(i, noise));
            timed!(t.truncate, state.truncate());
            timed!(t.rz, state.rzz(i, i + 1, tzz));
            timed!(t.truncate, state.truncate());
        }
    }
}

fn report(label: &str, t: &T, iters: u32) {
    let ms = |d: Duration| d.as_secs_f64() * 1e3 / iters as f64;
    println!(
        "{label:>4}: total={:.4}ms  pauli_error={:.4}  rx={:.4}  rz={:.4}  cnot={:.4}  truncate={:.4}",
        ms(t.total),
        ms(t.pauli_error),
        ms(t.rx),
        ms(t.rz),
        ms(t.cnot),
        ms(t.truncate)
    );
}

fn main() {
    let n = 12usize;
    let h = 1.0_f64;
    let dt = 0.1 / h;
    let time = 1.0 / h;
    let j = 1.0 / 8.0 * h;
    let steps = (time / dt) as usize;
    let tx = dt * h;
    let tzz = dt * j;
    let noise = [1e-4 / 4.0; 3];
    let iters: u32 = 300;

    let new_seed = build_new(n);
    let old_seed = build_old(n);

    // Interleave to average out any drift; accumulate per-class.
    let mut tn = T::default();
    let mut to = T::default();
    // warm up
    for _ in 0..20 {
        trotter_new(
            &mut new_seed.clone(),
            n,
            steps,
            tx,
            tzz,
            noise,
            &mut T::default(),
        );
        trotter_old(
            &mut old_seed.clone(),
            n,
            steps,
            tx,
            tzz,
            noise,
            &mut T::default(),
        );
    }
    for _ in 0..iters {
        let mut s = new_seed.clone();
        let g = Instant::now();
        trotter_new(&mut s, n, steps, tx, tzz, noise, &mut tn);
        tn.total += g.elapsed();
        let mut s = old_seed.clone();
        let g = Instant::now();
        trotter_old(&mut s, n, steps, tx, tzz, noise, &mut to);
        to.total += g.elapsed();
    }
    report("new", &tn, iters);
    report("old", &to, iters);
    let r = |a: Duration, b: Duration| a.as_secs_f64() / b.as_secs_f64();
    println!(
        "ratio new/old: total={:.3}  pauli_error={:.3}  rx={:.3}  rz={:.3}  cnot={:.3}  truncate={:.3}",
        r(tn.total, to.total),
        r(tn.pauli_error, to.pauli_error),
        r(tn.rx, to.rx),
        r(tn.rz, to.rz),
        r(tn.cnot, to.cnot),
        r(tn.truncate, to.truncate)
    );
    // Attribution: each new op class minus its old counterpart = contribution to gap.
    let gap = |a: Duration, b: Duration| (a.as_secs_f64() - b.as_secs_f64()) * 1e3 / iters as f64;
    println!(
        "gap contribution (new-old, ms/iter): pauli_error={:.4}  rx={:.4}  rz={:.4}  cnot={:.4}  truncate={:.4}  SUM_ops={:.4}  total={:.4}",
        gap(tn.pauli_error, to.pauli_error),
        gap(tn.rx, to.rx),
        gap(tn.rz, to.rz),
        gap(tn.cnot, to.cnot),
        gap(tn.truncate, to.truncate),
        gap(tn.pauli_error, to.pauli_error)
            + gap(tn.rx, to.rx)
            + gap(tn.rz, to.rz)
            + gap(tn.cnot, to.cnot)
            + gap(tn.truncate, to.truncate),
        gap(tn.total, to.total)
    );
}
