// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Statistical comparison between `GeneralizedTableauSum` (probability-
//! weighted branched representation, sampled at the end) and the pure
//! `GeneralizedTableau` (one stochastic trajectory per shot) on the noisy
//! magic-state-distillation circuit from `msd-noisy.rs`.
//!
//! Both backends model the same noisy quantum channel. The sum backend
//! drops branches with probability below `sum_cutoff`, which is the only
//! place the two distributions can diverge. With per-channel error rate
//! `p = 1e-4`, branches built from two independent errors have weight
//! ~`p^2 = 1e-8`; any `sum_cutoff > p^2` keeps the branch count tractable
//! while retaining the dominant (zero-error and first-order-error)
//! contributions.
//!
//! The script:
//!   1. Runs `n_shots` independent pure trajectories as the ground truth
//!      (rayon-parallelised, since shots are independent).
//!   2. For each cutoff in the sweep, builds the sum, samples `n_shots`
//!      shots (also rayon-parallelised), and compares the per-qubit
//!      marginals to the pure baseline.
//!
//! The marginal-L1 metric is `Σ_o |P_sum(o|q) - P_pure(o|q)|` per qubit,
//! summed over outcomes `o ∈ {0, 1, lost}`. Two independent finite-sample
//! estimates of the *same* distribution already differ by ~3·sqrt(0.5/N)
//! per qubit, which is the shot-noise floor below which a cutoff-induced
//! difference is not resolvable. The script picks `n_shots` large enough
//! that the floor is well below the smallest expected cutoff effect at
//! `p = 1e-4` (~few × 1e-3 from per-qubit noise contributions).
//!
//! Requires the `rayon` feature for tractable sum sampling at high
//! `n_shots`. Run with:
//!     cargo run --release --features rayon --example msd-noisy-compare

use std::time::Instant;

use rayon::prelude::*;

use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::measure_all::LossyMeasureAll;
use ppvm_tableau::prelude::*;
use ppvm_tableau_sum::data::GeneralizedTableauSum;
use ppvm_traits::traits::{Clifford, CliffordExtensions, Depolarizing, LossChannel, TGate};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

type Cfg = Byte8F64<2>;
type Idx = u128;
type Tab = GeneralizedTableau<Cfg, Idx>;
type TabSum = GeneralizedTableauSum<Cfg, Idx>;

fn encode_block<B>(tab: &mut B, qubits: &[usize], p_loss: f64, p_depolarize: f64)
where
    B: Clifford + CliffordExtensions + Depolarizing<Cfg> + LossChannel<Cfg>,
{
    if qubits.len() != 17 {
        return;
    }
    for i in [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16] {
        tab.sqrt_y(qubits[i]);
        tab.loss_channel(qubits[i], p_loss);
        tab.depolarize1(qubits[i], p_depolarize);
    }
    for [i, j] in [[1, 3], [7, 10], [12, 14], [13, 16]] {
        tab.cz(qubits[i], qubits[j]);
        tab.loss_channel(qubits[i], p_loss);
        tab.loss_channel(qubits[j], p_loss);
        tab.depolarize1(qubits[i], p_depolarize);
        tab.depolarize1(qubits[j], p_depolarize);
    }
    for i in [7, 16] {
        tab.sqrt_y_dag(qubits[i]);
        tab.loss_channel(qubits[i], p_loss);
        tab.depolarize1(qubits[i], p_depolarize);
    }
    for [i, j] in [[4, 7], [8, 10], [11, 14], [15, 16]] {
        tab.cz(qubits[i], qubits[j]);
        tab.loss_channel(qubits[i], p_loss);
        tab.loss_channel(qubits[j], p_loss);
        tab.depolarize1(qubits[i], p_depolarize);
        tab.depolarize1(qubits[j], p_depolarize);
    }
    for i in [4, 10, 14, 16] {
        tab.sqrt_y_dag(qubits[i]);
        tab.loss_channel(qubits[i], p_loss);
        tab.depolarize1(qubits[i], p_depolarize);
    }
    for [i, j] in [[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]] {
        tab.cz(qubits[i], qubits[j]);
        tab.loss_channel(qubits[i], p_loss);
        tab.loss_channel(qubits[j], p_loss);
        tab.depolarize1(qubits[i], p_depolarize);
        tab.depolarize1(qubits[j], p_depolarize);
    }
    for i in [3, 6, 9, 10, 12, 13] {
        tab.sqrt_y(qubits[i]);
        tab.loss_channel(qubits[i], p_loss);
        tab.depolarize1(qubits[i], p_depolarize);
    }
    for [i, j] in [[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]] {
        tab.cz(qubits[i], qubits[j]);
        tab.loss_channel(qubits[i], p_loss);
        tab.loss_channel(qubits[j], p_loss);
        tab.depolarize1(qubits[i], p_depolarize);
        tab.depolarize1(qubits[j], p_depolarize);
    }
    for i in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14] {
        tab.sqrt_y(qubits[i]);
        tab.loss_channel(qubits[i], p_loss);
        tab.depolarize1(qubits[i], p_depolarize);
    }
    for [i, j] in [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]] {
        tab.cz(qubits[i], qubits[j]);
        tab.loss_channel(qubits[i], p_loss);
        tab.loss_channel(qubits[j], p_loss);
        tab.depolarize1(qubits[i], p_depolarize);
        tab.depolarize1(qubits[j], p_depolarize);
    }
    for i in [0, 2, 5, 6, 8, 10, 12] {
        tab.sqrt_y_dag(qubits[i]);
        tab.loss_channel(qubits[i], p_loss);
        tab.depolarize1(qubits[i], p_depolarize);
    }
}

fn build_msd<B>(tab: &mut B, ql: &[&[usize]], p_loss: f64, p_depolarize: f64)
where
    B: Clifford + CliffordExtensions + TGate<Cfg> + Depolarizing<Cfg> + LossChannel<Cfg>,
{
    // Phase 1: Encoding (H + T + encode per block).
    for q in ql.iter() {
        let encoding_qubit = q[7];

        tab.h(encoding_qubit);
        tab.loss_channel(encoding_qubit, p_loss);
        tab.depolarize1(encoding_qubit, p_depolarize);

        tab.t(encoding_qubit);
        tab.loss_channel(encoding_qubit, p_loss);
        tab.depolarize1(encoding_qubit, p_depolarize);

        encode_block(tab, q, p_loss, p_depolarize);
    }

    // Phase 2: Middle gates (sqrt_x, cz, sqrt_y, sqrt_x_dag layers).
    for i in [0, 1, 4] {
        for q in ql[i] {
            tab.sqrt_x(*q);
            tab.loss_channel(*q, p_loss);
            tab.depolarize1(*q, p_depolarize);
        }
    }
    for (control, target) in ql[0].iter().zip(ql[1]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize1(*control, p_depolarize);
        tab.depolarize1(*target, p_depolarize);
    }
    for (control, target) in ql[2].iter().zip(ql[3]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize1(*control, p_depolarize);
        tab.depolarize1(*target, p_depolarize);
    }
    for q in ql[0] {
        tab.sqrt_y(*q);
        tab.loss_channel(*q, p_loss);
        tab.depolarize1(*q, p_depolarize);
    }
    for q in ql[3] {
        tab.sqrt_y(*q);
        tab.loss_channel(*q, p_loss);
        tab.depolarize1(*q, p_depolarize);
    }
    for (control, target) in ql[0].iter().zip(ql[2]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize1(*control, p_depolarize);
        tab.depolarize1(*target, p_depolarize);
    }
    for (control, target) in ql[3].iter().zip(ql[4]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize1(*control, p_depolarize);
        tab.depolarize1(*target, p_depolarize);
    }
    for q in ql[0] {
        tab.sqrt_x_dag(*q);
        tab.loss_channel(*q, p_loss);
        tab.depolarize1(*q, p_depolarize);
    }
    for (control, target) in ql[0].iter().zip(ql[4]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize1(*control, p_depolarize);
        tab.depolarize1(*target, p_depolarize);
    }
    for (control, target) in ql[1].iter().zip(ql[3]) {
        tab.cz(*control, *target);
        tab.loss_channel(*control, p_loss);
        tab.loss_channel(*target, p_loss);
        tab.depolarize1(*control, p_depolarize);
        tab.depolarize1(*target, p_depolarize);
    }
    for block in ql.iter().take(5) {
        for q in *block {
            tab.sqrt_x_dag(*q);
            tab.loss_channel(*q, p_loss);
            tab.depolarize1(*q, p_depolarize);
        }
    }
}

/// Per-qubit empirical distribution over `{0, 1, lost}`.
fn per_qubit_marginals(shots: &[Vec<Option<bool>>], n_qubits: usize) -> Vec<[f64; 3]> {
    let mut counts = vec![[0u64; 3]; n_qubits];
    for shot in shots {
        for (q, outcome) in shot.iter().enumerate() {
            let bin = match outcome {
                Some(false) => 0,
                Some(true) => 1,
                None => 2,
            };
            counts[q][bin] += 1;
        }
    }
    let n = shots.len() as f64;
    counts
        .iter()
        .map(|c| [c[0] as f64 / n, c[1] as f64 / n, c[2] as f64 / n])
        .collect()
}

/// (max, mean) of per-qubit L1 distance between two marginal sets.
fn l1_distance_stats(a: &[[f64; 3]], b: &[[f64; 3]]) -> (f64, f64) {
    let per_qubit: Vec<f64> = a
        .iter()
        .zip(b)
        .map(|(x, y)| x.iter().zip(y).map(|(xi, yi)| (xi - yi).abs()).sum::<f64>())
        .collect();
    let max = per_qubit.iter().copied().fold(0.0_f64, f64::max);
    let mean = per_qubit.iter().sum::<f64>() / per_qubit.len() as f64;
    (max, mean)
}

/// Run one sweep at fixed noise rate `p`. Builds the pure baseline and an
/// alt-seed pure run (for the shot-noise floor), then sweeps `cutoffs` on
/// the sum backend. Prints a comparison table.
fn run_sweep(
    label: &str,
    n_qubits: usize,
    n_shots: usize,
    p: f64,
    cutoffs: &[f64],
    pure_seed: u64,
    sum_seed: u64,
    ql: &[&[usize]],
) {
    println!("\n========================================================");
    println!("Sweep: {label}");
    println!("  p_loss = p_depolarize = {p:.0e}");
    println!(
        "  p^2 = {:.1e}   (cutoff must exceed this to keep branch count bounded)",
        p * p
    );

    let pure_start = Instant::now();
    let pure_shots: Vec<Vec<Option<bool>>> = (0..n_shots as u64)
        .into_par_iter()
        .map(|i| {
            let mut tab: Tab =
                GeneralizedTableau::new_with_seed(n_qubits, 1e-10, pure_seed.wrapping_add(i));
            build_msd(&mut tab, ql, p, p);
            tab.measure_all()
        })
        .collect();
    let pure_time = pure_start.elapsed();
    let pure_marginals = per_qubit_marginals(&pure_shots, n_qubits);
    println!("  pure baseline:  {} ms total", pure_time.as_millis());

    // Reference: pure-vs-pure with an independent seed. Any sum row that
    // sits at or below these L1 numbers is statistically indistinguishable
    // from pure; rows materially above them are showing a real cutoff
    // effect.
    let alt_shots: Vec<Vec<Option<bool>>> = (0..n_shots as u64)
        .into_par_iter()
        .map(|i| {
            let mut tab: Tab = GeneralizedTableau::new_with_seed(
                n_qubits,
                1e-10,
                pure_seed
                    .wrapping_add(i)
                    .wrapping_add(0xA5A5_A5A5_A5A5_A5A5),
            );
            build_msd(&mut tab, ql, p, p);
            tab.measure_all()
        })
        .collect();
    let alt_marginals = per_qubit_marginals(&alt_shots, n_qubits);
    let (alt_max, alt_mean) = l1_distance_stats(&alt_marginals, &pure_marginals);
    println!(
        "  pure-vs-pure:   max L1 = {:.4}   mean L1 = {:.4}   (shot-noise reference)",
        alt_max, alt_mean
    );

    println!(
        "\n  {:>10}  {:>9}  {:>11}  {:>12}  {:>9}  {:>9}",
        "cutoff", "branches", "build (ms)", "sample (ms)", "max L1", "mean L1"
    );
    for &cutoff in cutoffs {
        let build_start = Instant::now();
        let mut tab: TabSum =
            GeneralizedTableauSum::new_with_seed(n_qubits, 1e-10, cutoff, sum_seed);
        build_msd(&mut tab, ql, p, p);
        let build_time = build_start.elapsed();
        let branches = tab.len();

        let sample_start = Instant::now();
        let sum_shots = tab.sampler().sample_shots(n_shots);
        let sample_time = sample_start.elapsed();

        let sum_marginals = per_qubit_marginals(&sum_shots, n_qubits);
        let (max_l1, mean_l1) = l1_distance_stats(&sum_marginals, &pure_marginals);

        println!(
            "  {:>10.0e}  {:>9}  {:>11}  {:>12}  {:>9.4}  {:>9.4}",
            cutoff,
            branches,
            build_time.as_millis(),
            sample_time.as_millis(),
            max_l1,
            mean_l1
        );
    }
}

fn main() {
    let n_qubits = 85;
    // 1M shots puts the shot-noise floor around 7e-4 (per-bin std), well
    // below the expected per-qubit noise signal at p = 1e-4 (~few × 1e-3
    // from per-qubit error accumulation). With the rayon feature on,
    // sum sampling auto-dispatches to a parallel implementation so this
    // stays tractable.
    let n_shots = 1_000_000_usize;
    let pure_seed: u64 = 0xDEAD_BEEF;
    let sum_seed: u64 = 0x00C0_FFEE;

    let qubit_addrs: Vec<usize> = (0..n_qubits).collect();
    let ql: Vec<&[usize]> = qubit_addrs.chunks_exact(17).collect();

    println!("MSD-noisy comparison: GeneralizedTableauSum vs GeneralizedTableau");
    println!("  n_qubits         = {n_qubits}");
    println!("  n_shots          = {n_shots}");
    println!(
        "  shot-noise floor = ~{:.4} per-qubit L1 between two finite samples of the same dist",
        3.0 * (0.5_f64 / n_shots as f64).sqrt()
    );

    // Two regimes, same circuit shape:
    //
    // 1. p = 1e-4 matches msd-noisy.rs. At this noise level the per-qubit
    //    marginal differs from a noiseless run by ~few × 1e-3, which is
    //    above the 1M-shot per-qubit-L1 shot-noise floor (~2e-3, printed
    //    below). So we expect:
    //      - cutoffs >= p (1 branch, noiseless) to sit visibly above the
    //        floor (the missing-noise signal),
    //      - cutoffs in (p^2, p) (2025 branches) to sit at the floor.
    //
    // 2. p = 1e-2 amplifies the noise — the sum's truncation bias becomes
    //    the dominant signal (~0.25 max L1) regardless of cutoff in the
    //    safe regime. The convergence in this regime requires cutoff far
    //    below p^2 and is combinatorially out of reach; included for
    //    comparison only.
    run_sweep(
        "msd-noisy parameters (p = 1e-4)",
        n_qubits,
        n_shots,
        1e-4,
        &[1e-3, 5e-4, 1e-4, 1e-5, 1e-6, 1e-7, 5e-8],
        pure_seed,
        sum_seed,
        &ql,
    );
    run_sweep(
        "amplified noise (p = 1e-2)",
        n_qubits,
        n_shots,
        1e-2,
        // p^2 = 1e-4 is the boundary mentioned in the script's preamble:
        //   - cutoffs > p (1e-1, 1e-2) drop single-error structure → 1 branch
        //   - cutoffs in (p^2, p) keep first-order branches (~789..2025)
        //   - cutoffs at/below p^2 progressively pick up two-error branches
        //
        // 1e-5 sits just below (p/3)^2 = 1.1e-5, so two-error
        // depolarize-depolarize products start surviving — branches jump
        // to ~6700 and the build slows by ~20x. The marginal-L1 against
        // pure stays at ~0.26 either way: convergence here would require
        // cutoff ≪ (p/3)^3 ≈ 4e-8 to also capture three-error events,
        // which is combinatorially out of reach. With this list the full
        // example takes ~90s in release on an 8-core laptop.
        &[1e-1, 1e-2, 1e-3, 1e-4, 5e-5, 1e-5],
        pure_seed,
        sum_seed,
        &ql,
    );
}
