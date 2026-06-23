// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Truncation-error scaling study on a small random noisy circuit.
//!
//! Builds one seeded brickwork random circuit (single-qubit Cliffords + T,
//! interleaved with CZ in a brickwork pattern), with `loss_channel(p)` and
//! `depolarize1(p)` after every gate. The circuit is small enough that the
//! `GeneralizedTableauSum` can be built essentially exactly (cutoff far
//! below the smallest physically relevant branch), and large enough that
//! several orders of `p` show up as distinct branch populations.
//!
//! The same circuit is then re-run at a sweep of `sum_cutoff` values. For
//! each cutoff we report:
//!
//!   * Branch count after the build.
//!   * Build time.
//!   * **Dropped probability mass** — analytic, exact, computed from the
//!     reference (near-zero-cutoff) sum as `Σ_{p_i < cutoff} p_i`. This is
//!     a tight bound on the TVD between the reference and truncated
//!     distributions and is independent of shot noise.
//!   * **Sampled per-qubit L1** between `n_shots` shots taken from the
//!     truncated sum and the same number of shots from the reference. This
//!     is the empirical cross-check; it floors at the shot-noise level
//!     (~sqrt(1/n_shots)) and saturates at the dropped-mass bound.
//!
//! The "order k" of a cutoff is `k = -ln(cutoff)/ln(p)` rounded to one
//! decimal. At cutoff ≈ `p^k` we expect dropped mass to drop another
//! factor of ~`p` per integer increase in `k` (each new order has
//! branches with weight ~`p^k`).
//!
//! After the main p=0.01 sweep, a secondary analytic-only block runs at
//! p=1e-3 on the same circuit. At p=1e-3 the effects are far below any
//! affordable shot-noise floor, but the dropped-mass numbers still tell
//! the same scaling story — and show that for realistic noise rates a
//! loose cutoff is enormously safe.
//!
//! Requires the `rayon` feature for tractable shot sampling. Run with:
//!     cargo run --release --features rayon --example truncation-scaling

use std::time::Instant;

use rayon::prelude::*;

use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::measure_all::LossyMeasureAll;
use ppvm_tableau::prelude::*;
use ppvm_tableau_sum::data::GeneralizedTableauSum;
use ppvm_tableau_sum::storage::EntryStore;
use ppvm_traits::traits::{Clifford, CliffordExtensions, Depolarizing, LossChannel, TGate};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

type Cfg = Byte8F64<1>;
type Idx = u128;
type Tab = GeneralizedTableau<Cfg, Idx>;
type TabSum = GeneralizedTableauSum<Cfg, Idx>;

/// One layer of single-qubit gates followed by a CZ sublayer in brickwork
/// pattern. After every gate on a qubit we apply loss + depolarize on that
/// qubit. Single-qubit gate identities are drawn from a fixed 6-element
/// pool {H, S, S_dag, sqrt_x, sqrt_y, T}. CZ pair offsets alternate by
/// layer (even: (0,1),(2,3); odd: (1,2)) to give a standard brickwork.
fn apply_layer<B>(
    tab: &mut B,
    rng: &mut SmallRng,
    n_qubits: usize,
    layer_idx: usize,
    p_loss: f64,
    p_depolarize: f64,
) where
    B: Clifford + CliffordExtensions + Depolarizing<Cfg> + LossChannel<Cfg> + TGate<Cfg>,
{
    for q in 0..n_qubits {
        let g = rng.random::<u32>() % 6;
        match g {
            0 => tab.h(q),
            1 => tab.s(q),
            2 => tab.s_dag(q),
            3 => tab.sqrt_x(q),
            4 => tab.sqrt_y(q),
            _ => tab.t(q),
        }
        tab.loss_channel(q, p_loss);
        tab.depolarize1(q, p_depolarize);
    }

    let pairs: &[(usize, usize)] = if layer_idx % 2 == 0 {
        &[(0, 1), (2, 3)]
    } else {
        &[(1, 2)]
    };
    for &(a, b) in pairs {
        if a < n_qubits && b < n_qubits {
            tab.cz(a, b);
            tab.loss_channel(a, p_loss);
            tab.loss_channel(b, p_loss);
            tab.depolarize1(a, p_depolarize);
            tab.depolarize1(b, p_depolarize);
        }
    }
}

fn build_circuit<B>(
    tab: &mut B,
    circuit_seed: u64,
    n_qubits: usize,
    depth: usize,
    p_loss: f64,
    p_depolarize: f64,
) where
    B: Clifford + CliffordExtensions + Depolarizing<Cfg> + LossChannel<Cfg> + TGate<Cfg>,
{
    // Per-circuit RNG so gate choices reproduce regardless of the sum's
    // own RNG state.
    let mut rng = SmallRng::seed_from_u64(circuit_seed);
    for layer in 0..depth {
        apply_layer(tab, &mut rng, n_qubits, layer, p_loss, p_depolarize);
    }
}

/// Per-qubit empirical marginal over {0, 1, lost}.
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

/// Sum of probabilities of reference-sum branches whose individual weight
/// is below `cutoff`. This is the mass that a cutoff-`cutoff` build would
/// have dropped at the *final* step; in the perturbative regime this is
/// also an excellent estimate of the cumulative mass dropped throughout
/// the build (the difference is second-order in the dropped mass itself).
fn dropped_mass_below(reference: &TabSum, cutoff: f64) -> f64 {
    reference
        .entries
        .iter()
        .filter_map(|(_, p)| if *p < cutoff { Some(*p) } else { None })
        .sum()
}

/// Independent stochastic-trajectory shots of the same circuit, parallelised
/// across rayon. Each shot does its own full build + measure on a plain
/// `GeneralizedTableau`. Used only to print the empirical shot-noise floor.
fn pure_trajectory_shots(
    seeds: &[u64],
    n_qubits: usize,
    p: f64,
    circuit_seed: u64,
    depth: usize,
) -> Vec<Vec<Option<bool>>> {
    seeds
        .par_iter()
        .map(|&s| {
            let mut tab: Tab = GeneralizedTableau::new_with_seed(n_qubits, 1e-12, s);
            build_circuit(&mut tab, circuit_seed, n_qubits, depth, p, p);
            tab.measure_all()
        })
        .collect()
}

fn run_main_sweep(
    n_qubits: usize,
    depth: usize,
    p: f64,
    n_shots: usize,
    circuit_seed: u64,
    sum_seed: u64,
    reference_cutoff: f64,
    cutoffs: &[f64],
) {
    println!("\n========================================================");
    println!("Main sweep: random brickwork (Clifford + T), n={n_qubits}, depth={depth}, p={p:.0e}");
    println!(
        "  reference cutoff = {:.0e}   (effective 'exact' branched sum)",
        reference_cutoff
    );
    println!(
        "  shot-noise floor (per-qubit L1, 3 bins, {} shots) ≈ {:.4}",
        n_shots,
        3.0 * (0.5_f64 / n_shots as f64).sqrt()
    );

    // Build the reference (exact) branched sum once. All analytic
    // dropped-mass values are read directly off this single object.
    let ref_start = Instant::now();
    let mut reference: TabSum =
        GeneralizedTableauSum::new_with_seed(n_qubits, 1e-12, reference_cutoff, sum_seed);
    build_circuit(&mut reference, circuit_seed, n_qubits, depth, p, p);
    let ref_time = ref_start.elapsed();
    let ref_branches = reference.len();
    println!(
        "  reference build: {} branches in {} ms",
        ref_branches,
        ref_time.as_millis()
    );

    // Sample the reference once; reuse for every L1 comparison.
    let ref_sample_start = Instant::now();
    let ref_shots = reference.sampler().sample_shots(n_shots);
    let ref_sample_time = ref_sample_start.elapsed();
    let ref_marginals = per_qubit_marginals(&ref_shots, n_qubits);
    println!(
        "  reference sampling: {} ms ({:.0} ns/shot)",
        ref_sample_time.as_millis(),
        ref_sample_time.as_nanos() as f64 / n_shots as f64
    );

    // Independent pure-trajectory shots of the same circuit, to print the
    // actual two-sample shot-noise floor for this n_shots. Anything in the
    // sweep table at or below this is statistically indistinguishable from
    // the reference.
    let alt_seeds: Vec<u64> = (0..n_shots)
        .map(|i| {
            sum_seed
                .wrapping_add(i as u64)
                .wrapping_add(0xA5A5_A5A5_A5A5_A5A5)
        })
        .collect();
    let alt_shots = pure_trajectory_shots(&alt_seeds, n_qubits, p, circuit_seed, depth);
    let alt_marginals = per_qubit_marginals(&alt_shots, n_qubits);
    let (alt_max, alt_mean) = l1_distance_stats(&alt_marginals, &ref_marginals);
    println!(
        "  pure-vs-sum-reference:  max L1 = {:.5}   mean L1 = {:.5}   (empirical shot-noise floor)",
        alt_max, alt_mean
    );

    println!(
        "\n  {:>10}  {:>5}  {:>9}  {:>11}  {:>13}  {:>11}  {:>9}  {:>9}",
        "cutoff", "k", "branches", "build (ms)", "dropped mass", "sample(ms)", "max L1", "mean L1"
    );
    println!("  {}", "-".repeat(95));

    let mut rows = Vec::new();
    for &cutoff in cutoffs {
        let k = cutoff.ln() / p.ln();
        let dropped = dropped_mass_below(&reference, cutoff);

        let build_start = Instant::now();
        let mut tab: TabSum =
            GeneralizedTableauSum::new_with_seed(n_qubits, 1e-12, cutoff, sum_seed);
        build_circuit(&mut tab, circuit_seed, n_qubits, depth, p, p);
        let build_time = build_start.elapsed();
        let branches = tab.len();

        let sample_start = Instant::now();
        let shots = tab.sampler().sample_shots(n_shots);
        let sample_time = sample_start.elapsed();

        let marginals = per_qubit_marginals(&shots, n_qubits);
        let (max_l1, mean_l1) = l1_distance_stats(&marginals, &ref_marginals);

        println!(
            "  {:>10.0e}  {:>5.1}  {:>9}  {:>11}  {:>13.4e}  {:>11}  {:>9.5}  {:>9.5}",
            cutoff,
            k,
            branches,
            build_time.as_millis(),
            dropped,
            sample_time.as_millis(),
            max_l1,
            mean_l1,
        );

        rows.push((cutoff, k, dropped, max_l1));
    }

    // Headline answer: smallest order k whose dropped mass falls below a
    // target accuracy ε. We round k UP because we want a *sufficient*
    // cutoff, i.e. the next coarser order would already exceed ε.
    // Loosest cutoff in the sweep whose dropped mass is already below ε.
    // Iterate from largest cutoff to smallest; the first match is the
    // loosest sufficient one (and every smaller cutoff is also sufficient).
    println!("\n  Sufficient cutoff vs target accuracy ε (analytic, from dropped mass):");
    for &eps in &[1e-2_f64, 1e-3, 1e-4, 1e-5, 1e-6] {
        let candidate = rows.iter().find(|(_, _, dropped, _)| *dropped <= eps);
        match candidate {
            Some((c, k, dropped, _)) => println!(
                "    ε = {:.0e}  →  cutoff = {:.0e} suffices  (order k = {:.1}, dropped mass = {:.2e})",
                eps, c, k, dropped
            ),
            None => println!(
                "    ε = {:.0e}  →  not reached within the sweep (need smaller cutoff)",
                eps
            ),
        }
    }
}

fn run_secondary_sweep(
    n_qubits: usize,
    depth: usize,
    p: f64,
    circuit_seed: u64,
    sum_seed: u64,
    reference_cutoff: f64,
    cutoffs: &[f64],
) {
    println!("\n========================================================");
    println!("Secondary sweep (analytic-only): same circuit, p={p:.0e}");
    println!(
        "  reference cutoff = {:.0e};  no sampling — at p={p:.0e} all effects sit far below any affordable shot-noise floor.",
        reference_cutoff
    );

    let ref_start = Instant::now();
    let mut reference: TabSum =
        GeneralizedTableauSum::new_with_seed(n_qubits, 1e-14, reference_cutoff, sum_seed);
    build_circuit(&mut reference, circuit_seed, n_qubits, depth, p, p);
    let ref_time = ref_start.elapsed();
    println!(
        "  reference build: {} branches in {} ms",
        reference.len(),
        ref_time.as_millis()
    );

    println!("\n  {:>10}  {:>5}  {:>13}", "cutoff", "k", "dropped mass");
    println!("  {}", "-".repeat(40));
    for &cutoff in cutoffs {
        let k = cutoff.ln() / p.ln();
        let dropped = dropped_mass_below(&reference, cutoff);
        println!("  {:>10.0e}  {:>5.1}  {:>13.4e}", cutoff, k, dropped);
    }
}

fn main() {
    let n_qubits = 4;
    let depth = 6;
    let n_shots = 1_000_000;
    let circuit_seed: u64 = 0xC1FF_0177;
    let sum_seed: u64 = 0x5EED_BEEF;

    println!("Truncation-scaling study");
    println!("  n_qubits = {n_qubits}, depth = {depth}");
    println!("  random brickwork: single-qubit ∈ {{H,S,S†,√X,√Y,T}}, CZ on brickwork pairs");
    println!("  noise after every gate: loss_channel(p) then depolarize1(p)");

    // p = 0.01: main sweep — orders 1..~5 visible analytically, orders 1..2
    // resolvable by sampling at 1e6 shots.
    run_main_sweep(
        n_qubits,
        depth,
        1e-2,
        n_shots,
        circuit_seed,
        sum_seed,
        1e-14,
        &[1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10],
    );

    // p = 1e-3: realistic regime, analytic only. At p=1e-3 the dropped
    // mass at cutoff = 1e-4 is already ~p^? ; see what the sweep shows.
    run_secondary_sweep(
        n_qubits,
        depth,
        1e-3,
        circuit_seed,
        sum_seed,
        1e-18,
        &[1e-3, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16],
    );
}
