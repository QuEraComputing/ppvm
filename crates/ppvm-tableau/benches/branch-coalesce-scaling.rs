// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Branch-coalesce scaling: **sort-merge vs. hash-map**, head to head.
//!
//! PR #154 replaced the `FxHashMap` coalesce in the T-gate hot path
//! (`GeneralizedTableau::branch_with_coefficients`, `data.rs`) with a
//! sort-merge, measuring ~10× on `cultivation_d5`. That win was found on
//! one circuit; this benchmark asks the follow-up question directly:
//!
//! > Does the sort-merge advantage **persist as the branch count `m`
//! > grows**, and is there a regime where the hash coalesce wins again?
//!
//! Because PR #154 *deleted* the hash path from the default build (it now
//! survives only behind the `rayon` feature), there is no way to A/B the
//! two strategies through the public gate API. So both coalesce routines
//! are reimplemented here as free functions, faithful to their sources:
//!
//! * [`coalesce_sortmerge`] — a verbatim port of the sequential sort-merge
//!   in `branch_with_coefficients` (both the `u64`-packed fast path and the
//!   generic `(I, u32)` fallback), specialised to `IndexType = u128`.
//! * [`coalesce_hashmap`] — the pre-#154 `FxHashMap` coalesce, matching
//!   `branch_coefficients_seq` in `data.rs`.
//!
//! Both consume the **same real input**: a coefficient vector grown to an
//! exact size by applying H+T gates to a fresh 80-qubit tableau, plus the
//! genuine decomposition (`compute_decomposition` /
//! `odd_phase_destabilizer_mask`) of the next T gate. `verify_equivalence`
//! asserts the two produce identical coefficient sets before any timing, so
//! a drifted port fails loudly rather than benchmarking a lie.
//!
//! ## Mapping T gates ↔ branches
//!
//! T gates touch only the coefficient vector, never the tableau, so `k`
//! branching T gates on distinct qubits produce exactly `m = 2^k` branches
//! (no truncation here — the threshold is 0). The benchmark therefore
//! sweeps `m = 2^j` directly; that *is* the "number of T gates" axis. Real
//! circuits with truncation reach high T-gate counts at a bounded `m`, and
//! that bounded `m` is exactly what the sweep covers. 40 *untruncated*
//! branching T gates would be 2^40 ≈ 10^12 branches — out of reach for any
//! coalesce — so the honest variable is `m`, not the raw T count.
//!
//! ## Two collision regimes
//!
//! * **doubling** — the benched T gate flips a *fresh* index bit, so every
//!   branch lands on a new index: output `= 2·m`, zero merges. This is the
//!   canonical per-T-gate cost in a growing circuit.
//! * **merge** — the benched T gate flips an index bit the set is already
//!   closed under, so every branch coalesces onto an existing entry: output
//!   `= m`, all merges. This is the collision-heavy regime (the flavour of
//!   the measurement case-a path) and the most likely place for the hash
//!   coalesce to claw back ground.
//!
//! ## Reading the results
//!
//! Within each group (`branch-coalesce-doubling`, `branch-coalesce-merge`)
//! the `hashmap` and `sortmerge` lines are parameterised by `m`. Compare
//! them at each `m` and watch for a crossover. The `m = 32768 → 65536` step
//! also straddles the packed-path cutoff (`m ≤ 65535`): above it the
//! sort-merge drops to its generic `(I, u32)` fallback, so any change in the
//! gap there is the packing's contribution.
//!
//! Run:
//! ```bash
//! cargo bench -p ppvm-tableau --bench branch-coalesce-scaling
//! ```
//! Push the ceiling (default `m ≤ 2^20`) with e.g. `PPVM_BRANCH_MAX_EXP=22`.

use std::cmp::Ordering;
use std::hint::black_box;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use fxhash::FxHashMap;
use num::complex::{Complex, Complex64};
use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;

/// 128-bit index / 2×u64 storage — the 80-qubit regime used elsewhere
/// (`measure-scaling.rs`, `profile_scaling.rs`).
type Tab = GeneralizedTableau<Byte8F64<2>, u128>;

const N_QUBITS: usize = 80;

/// `exp(iπ/8)·cos(π/8)` — the non-branch (`coefficient_factor`) weight a
/// `T` gate applies. Copied from `gates/tgate.rs`.
const COS: Complex64 = Complex {
    re: 0.8535533905932737,
    im: 0.3535533905932738,
};
/// `-i·exp(iπ/8)·sin(π/8)` — the branch (`branch_factor`) weight. From
/// `gates/tgate.rs`.
const SIN: Complex64 = Complex {
    re: 0.14644660940672624,
    im: -0.3535533905932738,
};

/// Phase index → unit complex, matching `COMPLEX_PHASE_CONVERSION` in `data.rs`.
const PHASE: [Complex64; 4] = [
    Complex { re: 1.0, im: 0.0 },
    Complex { re: 0.0, im: 1.0 },
    Complex { re: -1.0, im: 0.0 },
    Complex { re: 0.0, im: -1.0 },
];

/// The decomposition + scaling factors a single branching T gate applies to
/// every coefficient. Mirrors the arguments threaded through
/// `branch_with_coefficients`.
#[derive(Clone, Copy)]
struct Params {
    stab_bits: u128,
    destab_bits: u128,
    odd_mask: u128,
    phase_decomp: u8,
    coefficient_factor: Complex64,
    branch_factor: Complex64,
    cutoff_sq: f64,
}

/// Verbatim port of `compute_phase_with_mask_static` (`data.rs`), specialised
/// to `u128`.
#[inline]
fn compute_phase_with_mask(destab_bits: u128, basis: u128, stab_bits: u128, odd_mask: u128) -> u8 {
    let mut phase = (2 * ((destab_bits & basis).count_ones() as u8)) % 4;
    let active = basis & stab_bits;
    let parity = (active & odd_mask).count_ones() % 2;
    phase = (phase + 2 * parity as u8) % 4;
    phase
}

/// Pre-#154 hash coalesce: one `FxHashMap` probe per branch + non-branch
/// contribution, then a magnitude-cutoff sweep into the output vector.
/// Matches `branch_coefficients_seq` (`data.rs`).
fn coalesce_hashmap(input: &[(Complex64, u128)], p: &Params) -> Vec<(Complex64, u128)> {
    let mut map: FxHashMap<u128, Complex64> =
        FxHashMap::with_capacity_and_hasher(2 * input.len(), Default::default());
    for &(coeff, idx) in input {
        let branch_index = idx ^ p.stab_bits;
        let bpc = compute_phase_with_mask(p.destab_bits, idx, p.stab_bits, p.odd_mask);
        let branch_phase = (bpc + p.phase_decomp) % 4;
        let pf = PHASE[branch_phase as usize];
        let branch_coefficient = pf * coeff * p.branch_factor;
        let nonbranch_coefficient = coeff * p.coefficient_factor;
        *map.entry(branch_index).or_insert(Complex64::new(0.0, 0.0)) += branch_coefficient;
        *map.entry(idx).or_insert(Complex64::new(0.0, 0.0)) += nonbranch_coefficient;
    }
    let mut out = Vec::with_capacity(map.len());
    for (idx, coeff) in map {
        if coeff.norm_sqr() > p.cutoff_sq {
            out.push((coeff, idx));
        }
    }
    out
}

/// Verbatim port of the sequential sort-merge in `branch_with_coefficients`
/// (`data.rs`), specialised to `u128`. Keeps both the `u64`-packed fast path
/// (engaged when `m ≤ 65535` and every branch key fits in 47 bits) and the
/// generic `(u128, u32)` fallback, so the benchmark exercises whichever path
/// the real code would take at a given `m`.
fn coalesce_sortmerge(input: &[(Complex64, u128)], p: &Params) -> Vec<(Complex64, u128)> {
    let n = input.len();
    let cutoff_sq = p.cutoff_sq;

    let mut nb: Vec<(u128, Complex64)> = Vec::with_capacity(n);
    let mut brv: Vec<Complex64> = Vec::with_capacity(n);
    let mut packed: Vec<u64> = Vec::with_capacity(n);
    let mut packable = n <= 0xFFFF;
    let mut nb_sorted = true;
    let mut prev: Option<u128> = None;

    for (pos, &(coeff, idx)) in (0_u32..).zip(input) {
        let branch_index = idx ^ p.stab_bits;
        let bpc = compute_phase_with_mask(p.destab_bits, idx, p.stab_bits, p.odd_mask);
        let branch_phase = (bpc + p.phase_decomp) % 4;
        let pf = PHASE[branch_phase as usize];
        brv.push(pf * coeff * p.branch_factor);
        if branch_index < (1u128 << 47) {
            packed.push(((branch_index as u64) << 16) | (pos as u64));
        } else {
            packable = false;
            packed.push(pos as u64);
        }
        nb.push((idx, coeff * p.coefficient_factor));
        if let Some(pp) = prev
            && idx < pp
        {
            nb_sorted = false;
        }
        prev = Some(idx);
    }

    let mut out: Vec<(Complex64, u128)> = Vec::with_capacity(nb.len() + brv.len());
    let mut i = 0;

    if packable {
        if !nb_sorted {
            nb.sort_unstable_by_key(|a| a.0);
        }
        packed.sort_unstable();
        let mut j = 0;
        while i < nb.len() && j < packed.len() {
            let bp = (packed[j] & 0xFFFF) as usize;
            let bk = (packed[j] >> 16) as u128;
            match nb[i].0.cmp(&bk) {
                Ordering::Less => {
                    if nb[i].1.norm_sqr() > cutoff_sq {
                        out.push((nb[i].1, nb[i].0));
                    }
                    i += 1;
                }
                Ordering::Greater => {
                    let v = brv[bp];
                    if v.norm_sqr() > cutoff_sq {
                        out.push((v, bk));
                    }
                    j += 1;
                }
                Ordering::Equal => {
                    let mut sv = nb[i].1;
                    sv += brv[bp];
                    if sv.norm_sqr() > cutoff_sq {
                        out.push((sv, nb[i].0));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        while j < packed.len() {
            let bp = (packed[j] & 0xFFFF) as usize;
            let bk = (packed[j] >> 16) as u128;
            let v = brv[bp];
            if v.norm_sqr() > cutoff_sq {
                out.push((v, bk));
            }
            j += 1;
        }
    } else {
        let mut brk: Vec<(u128, u32)> = (0_u32..)
            .zip(nb.iter())
            .map(|(pp, &(idx, _))| (idx ^ p.stab_bits, pp))
            .collect();
        if !nb_sorted {
            nb.sort_unstable_by_key(|a| a.0);
        }
        brk.sort_unstable_by_key(|a| a.0);
        let mut j = 0;
        while i < nb.len() && j < brk.len() {
            let (bk, bp) = brk[j];
            match nb[i].0.cmp(&bk) {
                Ordering::Less => {
                    if nb[i].1.norm_sqr() > cutoff_sq {
                        out.push((nb[i].1, nb[i].0));
                    }
                    i += 1;
                }
                Ordering::Greater => {
                    let v = brv[bp as usize];
                    if v.norm_sqr() > cutoff_sq {
                        out.push((v, bk));
                    }
                    j += 1;
                }
                Ordering::Equal => {
                    let mut sv = nb[i].1;
                    sv += brv[bp as usize];
                    if sv.norm_sqr() > cutoff_sq {
                        out.push((sv, nb[i].0));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        while j < brk.len() {
            let (bk, bp) = brk[j];
            let v = brv[bp as usize];
            if v.norm_sqr() > cutoff_sq {
                out.push((v, bk));
            }
            j += 1;
        }
    }
    while i < nb.len() {
        if nb[i].1.norm_sqr() > cutoff_sq {
            out.push((nb[i].1, nb[i].0));
        }
        i += 1;
    }
    out
}

/// Build a real coefficient vector of size `m = 2^j` plus the decomposition
/// of one more T gate. `fresh_target` selects the collision regime:
/// * `true`  — the next T flips a *fresh* qubit bit ⇒ doubling (output `2·m`).
/// * `false` — the next T flips an *existing* bit  ⇒ all-merge (output `m`).
fn build(n_qubits: usize, j: usize, fresh_target: bool) -> (Vec<(Complex64, u128)>, Params) {
    let mut tab: Tab = GeneralizedTableau::new(n_qubits, 0.0);
    // `j` branching T gates on distinct fresh qubits ⇒ exactly 2^j branches.
    for i in 0..j {
        tab.h(i);
        tab.t(i);
    }
    let target = if fresh_target {
        tab.h(j); // a qubit not yet branched on
        j
    } else {
        0 // qubit 0 is already branched — the set is closed under its bit
    };
    let (phase_decomp, stab_bits, destab_bits) = tab.compute_decomposition(target, Pauli::Z);
    let odd_mask = tab.odd_phase_destabilizer_mask();
    let input = tab.coefficients.clone();
    (
        input,
        Params {
            stab_bits,
            destab_bits,
            odd_mask,
            phase_decomp,
            coefficient_factor: COS,
            branch_factor: SIN,
            cutoff_sq: 0.0,
        },
    )
}

/// Assert the two coalesce routines agree on a real input — same set of
/// indices, same coefficients to 1e-9. Guards against the ports drifting from
/// their sources.
fn assert_equivalent(label: &str, a: &[(Complex64, u128)], b: &[(Complex64, u128)]) {
    assert_eq!(a.len(), b.len(), "{label}: output sizes differ");
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    a.sort_unstable_by_key(|x| x.1);
    b.sort_unstable_by_key(|x| x.1);
    for (x, y) in a.iter().zip(b.iter()) {
        assert_eq!(x.1, y.1, "{label}: index mismatch");
        assert!(
            (x.0 - y.0).norm() < 1e-9,
            "{label}: coeff mismatch at index {}: {:?} vs {:?}",
            x.1,
            x.0,
            y.0
        );
    }
}

fn verify_equivalence() {
    for (regime, fresh) in [("doubling", true), ("merge", false)] {
        let (input, params) = build(N_QUBITS, 6, fresh);
        let hm = coalesce_hashmap(&input, &params);
        let sm = coalesce_sortmerge(&input, &params);
        assert_equivalent(regime, &hm, &sm);
    }
}

/// Exponents `j` (with `m = 2^j`) chosen to bracket the packed-path cutoff
/// (`m ≤ 65535`, i.e. j ≤ 15) and span small → large branch counts.
fn exponents() -> Vec<usize> {
    let max = std::env::var("PPVM_BRANCH_MAX_EXP")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(20);
    [2usize, 5, 8, 11, 14, 15, 16, 18, 20]
        .into_iter()
        .filter(|&j| j <= max)
        .collect()
}

fn bench_scenario(c: &mut Criterion, group_name: &str, fresh_target: bool) {
    let mut group = c.benchmark_group(group_name);
    for j in exponents() {
        let (input, params) = build(N_QUBITS, j, fresh_target);
        let m = input.len() as u64;
        // ns/element is the cleaner scaling readout than raw wall time.
        group.throughput(Throughput::Elements(m));

        group.bench_with_input(BenchmarkId::new("hashmap", m), &m, |b, _| {
            b.iter(|| black_box(coalesce_hashmap(black_box(&input), &params)));
        });
        group.bench_with_input(BenchmarkId::new("sortmerge", m), &m, |b, _| {
            b.iter(|| black_box(coalesce_sortmerge(black_box(&input), &params)));
        });
    }
    group.finish();
}

fn bench_branch_coalesce(c: &mut Criterion) {
    verify_equivalence();
    bench_scenario(c, "branch-coalesce-doubling", true);
    bench_scenario(c, "branch-coalesce-merge", false);
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(30);
    targets = bench_branch_coalesce
}
criterion_main!(benches);
