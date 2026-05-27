//! Decides whether parallelizing the noise-channel "Phase 1" (per-entry
//! `fork` + Pauli, before the merge) is worth implementing.
//!
//! Group A ("phase-split"): on a realistic populated sum, times Phase 1
//! (serial vs parallel) against Phase 2 (`insert_or_merge_batch`). The
//! Phase1/Phase2 ratio is the Amdahl ceiling — if the merge dominates,
//! parallel forks can't help much.
//!
//! Group B ("phase1-scaling"): controlled microbench. Clones one
//! representative tableau N times and times the serial vs parallel fork-loop
//! across entry counts and qubit counts, to find the crossover N.
//!
//! The parallel variant is prototyped inline here (no library change) — the
//! bench owns the concrete `Vec`, so it can call `par_iter` directly.

use std::hint::black_box;
use std::time::Duration;

use criterion::{BatchSize, BenchmarkId, Criterion, criterion_group, criterion_main};
use num::complex::Complex64;
use rayon::prelude::*;
use smallvec::SmallVec;

use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_runtime::traits::{Clifford, CliffordExtensions, Depolarizing, LossChannel, TGate};
use ppvm_tableau::data::GeneralizedTableau;
use ppvm_tableau_sum::{
    data::GeneralizedTableauSum,
    storage::{EntryStore, vec::VecStorage, word_fingerprint},
};

type Cfg = Byte8F64<2>;
type Idx = u128;
type Coef = Vec<(Complex64, Idx)>;
type Tableau = GeneralizedTableau<Cfg, Idx, Coef>;
type VecBackedSum = GeneralizedTableauSum<Cfg, Idx, Coef, VecStorage<Cfg, Idx, Coef>>;

const N_QUBITS: usize = 17;
const P_LOSS: f64 = 1e-4;
const P_DEPOL: f64 = 1e-4;
const SUM_CUTOFF: f64 = 1e-7;
const COEFF_THRESHOLD: f64 = 1e-10;
const SEED: u64 = 42;
const ADDR: usize = 0;

/// T-gates used to grow each tableau's *internal* coefficient count (K) in
/// Group B. Kept deliberately low: we want realistic, cheap forks while we
/// scale the number of *sum* terms (entries, N) up to ~10k. K ≈ 2^T_GATES.
const T_GATES: usize = 3;

/// Group C builds a realistic high-entry-count sum (genuinely distinct
/// entries from a noisy circuit, not clones) to time the *full* noise call —
/// fork-loop + merge — and see whether the serial merge caps parallel gains.
const HIGH_NQ: usize = 30;
const HIGH_P: f64 = 0.05;
const HIGH_CUTOFF: f64 = 1e-5;
const TARGET_N: usize = 4000;

/// Depolarize-shaped Phase 1: 3 forks + X/Y/Z per non-lost entry. Serial.
fn gen_serial(entries: &[(Tableau, f64)], addr0: usize) -> Vec<(Tableau, f64, u64)> {
    entries
        .iter()
        .enumerate()
        .flat_map(|(i, (tab, p_sum))| branch_entry(i as u64, tab, *p_sum, addr0))
        .collect()
}

/// Same work as `gen_serial`, but the per-entry fork loop runs on rayon.
fn gen_parallel(entries: &[(Tableau, f64)], addr0: usize) -> Vec<(Tableau, f64, u64)> {
    entries
        .par_iter()
        .enumerate()
        .flat_map_iter(|(i, (tab, p_sum))| branch_entry(i as u64, tab, *p_sum, addr0))
        .collect()
}

/// The independent per-entry work. SmallVec<[_; 3]> keeps the 3 branches
/// inline (no per-entry heap alloc), so serial and parallel are comparable.
/// Each branch carries the parent's word-fingerprint (X/Y/Z leave words
/// unchanged), matching what the noise channels now feed the merge.
#[inline]
fn branch_entry(
    seed: u64,
    tab: &Tableau,
    p_sum: f64,
    addr0: usize,
) -> SmallVec<[(Tableau, f64, u64); 3]> {
    let mut out: SmallVec<[(Tableau, f64, u64); 3]> = SmallVec::new();
    if tab.is_lost[addr0] {
        return out;
    }
    let word_fp = word_fingerprint(tab);
    let mut bx = tab.fork(Some(seed));
    let mut by = tab.fork(Some(seed ^ 1));
    let mut bz = tab.fork(Some(seed ^ 2));
    bx.x(addr0);
    by.y(addr0);
    bz.z(addr0);
    let w = p_sum / 3.0;
    out.push((bx, w, word_fp));
    out.push((by, w, word_fp));
    out.push((bz, w, word_fp));
    out
}

/// msd-noisy-shaped circuit (mirrors storage_compare) to grow a realistic
/// entry count and tableau variety.
fn apply_circuit(tab: &mut VecBackedSum) {
    for q in 0..N_QUBITS {
        tab.sqrt_y(q);
        tab.loss_channel(q, P_LOSS);
        tab.depolarize(q, P_DEPOL);
    }
    for q in [0, 7, 12] {
        tab.t(q);
        tab.loss_channel(q, P_LOSS);
        tab.depolarize(q, P_DEPOL);
    }
    for [i, j] in [
        [1, 3],
        [7, 10],
        [12, 14],
        [13, 16],
        [4, 7],
        [8, 10],
        [11, 14],
        [15, 16],
    ] {
        tab.cz(i, j);
        tab.loss_channel(i, P_LOSS);
        tab.loss_channel(j, P_LOSS);
        tab.depolarize(i, P_DEPOL);
        tab.depolarize(j, P_DEPOL);
    }
    for q in [0, 2, 5, 6, 8, 10, 12] {
        tab.sqrt_y_adj(q);
        tab.loss_channel(q, P_LOSS);
    }
}

fn populated_sum() -> VecBackedSum {
    let mut tab = GeneralizedTableauSum::new_with_seed(N_QUBITS, COEFF_THRESHOLD, SUM_CUTOFF, SEED);
    apply_circuit(&mut tab);
    tab
}

/// One representative tableau at `n_qubits`, with internal coefficient count
/// (K) grown by `t_gates` T gates so fork's clone cost is realistic.
fn base_tableau(n_qubits: usize, t_gates: usize) -> Tableau {
    let mut sum: VecBackedSum =
        GeneralizedTableauSum::new_with_seed(n_qubits, COEFF_THRESHOLD, SUM_CUTOFF, SEED);
    for q in 0..n_qubits {
        sum.h(q);
    }
    for k in 0..t_gates {
        sum.t(k % n_qubits);
    }
    sum.entries.entries[0].0.clone()
}

fn phase_split(c: &mut Criterion) {
    let sum = populated_sum();
    let entries = sum.entries.entries.clone();
    let storage = sum.entries.clone();
    let n = entries.len();
    let k_sample = entries.first().map(|e| e.0.coefficients.len()).unwrap_or(0);
    let threads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(0);
    let branches = gen_serial(&entries, ADDR);
    eprintln!(
        "[phase-split] entries N={n}, sample K={k_sample}, branches={}, threads={threads}",
        branches.len()
    );

    let mut group = c.benchmark_group("phase-split");
    group.bench_function("phase1-gen-serial", |b| {
        b.iter(|| black_box(gen_serial(black_box(&entries), ADDR)))
    });
    group.bench_function("phase1-gen-parallel", |b| {
        b.iter(|| black_box(gen_parallel(black_box(&entries), ADDR)))
    });
    group.bench_function("phase2-merge", |b| {
        b.iter_batched(
            || (storage.clone(), branches.clone()),
            |(mut s, br)| black_box(s.insert_or_merge_batch(br, &SUM_CUTOFF)),
            BatchSize::SmallInput,
        )
    });
    group.finish();
}

fn phase1_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("phase1-scaling");
    for &nq in &[17usize, 64usize] {
        let base = base_tableau(nq, T_GATES);
        let k = base.coefficients.len();
        eprintln!("[phase1-scaling] nq={nq}, K={k}");
        for &n in &[64usize, 512, 2048, 10000] {
            let entries: Vec<(Tableau, f64)> =
                (0..n).map(|_| (base.clone(), 1.0 / n as f64)).collect();
            group.bench_with_input(
                BenchmarkId::new(format!("serial/nq{nq}/k{k}"), n),
                &entries,
                |b, e| b.iter(|| black_box(gen_serial(black_box(e), ADDR))),
            );
            group.bench_with_input(
                BenchmarkId::new(format!("parallel/nq{nq}/k{k}"), n),
                &entries,
                |b, e| b.iter(|| black_box(gen_parallel(black_box(e), ADDR))),
            );
        }
    }
    group.finish();
}

/// Grow a sum with genuinely distinct entries until it reaches ~TARGET_N,
/// using a noisy Clifford+T circuit. Built once in setup (untimed).
fn high_n_sum() -> VecBackedSum {
    let mut tab = GeneralizedTableauSum::new_with_seed(HIGH_NQ, COEFF_THRESHOLD, HIGH_CUTOFF, SEED);
    'outer: for layer in 0..100 {
        for q in 0..HIGH_NQ {
            tab.sqrt_y(q);
            tab.depolarize(q, HIGH_P);
            if tab.len() >= TARGET_N {
                break 'outer;
            }
        }
        for q in (0..HIGH_NQ - 1).step_by(2) {
            tab.cz(q, q + 1);
        }
        tab.t(layer % HIGH_NQ);
        tab.loss_channel(layer % HIGH_NQ, HIGH_P);
        if tab.len() >= TARGET_N {
            break;
        }
    }
    tab
}

/// Full noise call (depolarize-shaped) at high N: fork-loop + merge, serial
/// vs parallel fork-loop (merge stays serial in both). The serial-vs-parallel
/// ratio here is the *real* end-to-end gain, with the merge included.
fn phase_full(c: &mut Criterion) {
    let sum = high_n_sum();
    let storage = sum.entries.clone();
    let n = storage.entries.len();
    let branches = gen_serial(&storage.entries, ADDR);
    eprintln!(
        "[phase-full] N={n}, branches={}, threads={}",
        branches.len(),
        std::thread::available_parallelism()
            .map(|x| x.get())
            .unwrap_or(0)
    );

    let mut group = c.benchmark_group("phase-full");
    group.bench_function("gen-serial", |b| {
        b.iter(|| black_box(gen_serial(black_box(&storage.entries), ADDR)))
    });
    group.bench_function("gen-parallel", |b| {
        b.iter(|| black_box(gen_parallel(black_box(&storage.entries), ADDR)))
    });
    group.bench_function("merge", |b| {
        b.iter_batched(
            || (storage.clone(), branches.clone()),
            |(mut s, br)| black_box(s.insert_or_merge_batch(br, &HIGH_CUTOFF)),
            BatchSize::LargeInput,
        )
    });
    group.bench_function("full-serial", |b| {
        b.iter_batched(
            || storage.clone(),
            |mut s| {
                let br = gen_serial(&s.entries, ADDR);
                black_box(s.insert_or_merge_batch(br, &HIGH_CUTOFF))
            },
            BatchSize::LargeInput,
        )
    });
    group.bench_function("full-parallel", |b| {
        b.iter_batched(
            || storage.clone(),
            |mut s| {
                let br = gen_parallel(&s.entries, ADDR);
                black_box(s.insert_or_merge_batch(br, &HIGH_CUTOFF))
            },
            BatchSize::LargeInput,
        )
    });
    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(2))
        .sample_size(20);
    targets = phase_split, phase1_scaling, phase_full
}
criterion_main!(benches);
