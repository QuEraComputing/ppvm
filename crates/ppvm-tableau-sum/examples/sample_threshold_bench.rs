// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Quick heuristic benchmark to find the n_shots crossover between
//! the serial and parallel implementations of `Sampler::sample_shots`.

use std::time::Instant;

use num::complex::Complex64;

use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau_sum::data::GeneralizedTableauSum;
use ppvm_traits::traits::{Clifford, CliffordExtensions, Depolarizing, LossChannel, TGate};

type Cfg = Byte8F64<2>;
type Idx = u128;
type Coef = Vec<(Complex64, Idx)>;
type GTabSum = GeneralizedTableauSum<Cfg, Idx, Coef>;

fn time<R>(reps: u32, mut f: impl FnMut() -> R) -> f64 {
    let _ = f();
    let start = Instant::now();
    for _ in 0..reps {
        std::hint::black_box(f());
    }
    start.elapsed().as_secs_f64() / reps as f64
}

fn light_workload() -> GTabSum {
    let mut tab: GTabSum = GeneralizedTableauSum::new_with_seed(2, 1e-12, 1e-10, 42);
    tab.h(0);
    tab.cnot(0, 1);
    tab
}

fn medium_workload() -> GTabSum {
    let mut tab: GTabSum = GeneralizedTableauSum::new_with_seed(17, 1e-12, 1e-6, 42);
    for q in 0..17 {
        tab.sqrt_y(q);
    }
    for q in [0, 7, 12] {
        tab.t(q);
        tab.loss_channel(q, 1e-3);
        tab.depolarize1(q, 1e-3);
    }
    for [i, j] in [[1, 3], [7, 10], [12, 14], [13, 16], [4, 7]] {
        tab.cz(i, j);
    }
    tab
}

fn heavy_workload() -> GTabSum {
    let mut tab: GTabSum = GeneralizedTableauSum::new_with_seed(30, 1e-12, 1e-5, 42);
    for q in 0..30 {
        tab.sqrt_y(q);
    }
    for q in [0, 5, 10, 15, 20, 25] {
        tab.t(q);
        tab.depolarize1(q, 1e-3);
    }
    for [i, j] in [
        [1, 3],
        [7, 10],
        [12, 14],
        [13, 16],
        [4, 7],
        [20, 22],
        [24, 27],
    ] {
        tab.cz(i, j);
    }
    for q in 0..30 {
        tab.loss_channel(q, 1e-4);
    }
    tab
}

fn measure_crossover(label: &str, mut tab: GTabSum, shots: &[usize]) {
    println!("\n== Workload: {label} ==");
    println!("  n_entries = {}", tab.len());
    let mut sampler = tab.sampler();

    println!(
        "  {:>8}  {:>13}  {:>13}  {:>8}",
        "n_shots", "serial (us)", "parallel (us)", "speedup"
    );
    for &n in shots {
        let reps = (200_000 / n.max(1)).clamp(5, 500) as u32;
        let t_ser = time(reps, || sampler.sample_shots_serial(n));
        let t_par = time(reps, || sampler.sample_shots_parallel(n));
        let speedup = t_ser / t_par;
        println!(
            "  {:>8}  {:>13.2}  {:>13.2}  {:>7.2}x",
            n,
            t_ser * 1e6,
            t_par * 1e6,
            speedup
        );
    }
}

fn main() {
    let n_threads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(0);
    println!("available_parallelism = {n_threads}");

    let shots = [1usize, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096];

    measure_crossover("light (Bell pair, 2q, 1 entry)", light_workload(), &shots);
    measure_crossover("medium (17q + noise)", medium_workload(), &shots);
    measure_crossover("heavy (30q + many T)", heavy_workload(), &shots);
}
