//! Profile the time breakdown inside branch_with_coefficients:
//! computation (phase + multiply) vs HashMap accumulation.

use std::time::Instant;

use fxhash::FxHashMap as HashMap;
use num::Zero;
use num::complex::{Complex, Complex64};

fn main() {
    // Simulate the coefficient branching loop at various sizes
    println!(
        "{:>10}  {:>10}  {:>10}  {:>10}  {:>6}",
        "N", "compute", "accum", "total", "%comp"
    );

    for n in [2048, 8192, 32768, 131072] {
        // Create fake coefficient data
        let items: Vec<(Complex<f64>, u128)> = (0..n)
            .map(|i| (Complex::new(1.0 / n as f64, 0.0), i as u128))
            .collect();

        let stab_anticomm_bits: u128 = 0b1010_1010;
        let destab_anticomm_bits: u128 = 0b0101_0101;
        let odd_phase_mask: u128 = 0b1100_1100;
        let phase_decomp: u8 = 1;
        let coefficient_factor = Complex::new(0.85, 0.35);
        let branch_factor = Complex::new(0.15, -0.35);

        let phase_table: [Complex64; 4] = [
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 1.0),
            Complex64::new(-1.0, 0.0),
            Complex64::new(0.0, -1.0),
        ];

        let n_runs = 10;

        // Phase 1: compute only (write to Vec, no HashMap)
        let mut compute_ns = 0u128;
        for _ in 0..n_runs {
            let mut results = Vec::with_capacity(items.len());
            let t0 = Instant::now();
            for &(coeff, idx) in &items {
                let branch_index = idx ^ stab_anticomm_bits;
                let symplectic = (destab_anticomm_bits & idx).count_ones();
                let mut phase = (2 * symplectic as u8) % 4;
                let active = idx & stab_anticomm_bits;
                let parity = (active & odd_phase_mask).count_ones() % 2;
                phase = (phase + 2 * parity as u8) % 4;
                let branch_phase = (phase + phase_decomp) % 4;
                let phase_factor: Complex<f64> = phase_table[branch_phase as usize].into();
                let branch_coeff = phase_factor * coeff * branch_factor;
                let nonbranch_coeff = coeff * coefficient_factor;
                results.push((branch_index, branch_coeff, idx, nonbranch_coeff));
            }
            compute_ns += t0.elapsed().as_nanos();
            std::hint::black_box(&results);
        }

        // Phase 2: accumulate only (from pre-computed data)
        let mut precomputed: Vec<(u128, Complex<f64>, u128, Complex<f64>)> =
            Vec::with_capacity(items.len());
        for &(coeff, idx) in &items {
            let branch_index = idx ^ stab_anticomm_bits;
            let symplectic = (destab_anticomm_bits & idx).count_ones();
            let mut phase = (2 * symplectic as u8) % 4;
            let active = idx & stab_anticomm_bits;
            let parity = (active & odd_phase_mask).count_ones() % 2;
            phase = (phase + 2 * parity as u8) % 4;
            let branch_phase = (phase + phase_decomp) % 4;
            let phase_factor: Complex<f64> = phase_table[branch_phase as usize].into();
            precomputed.push((
                branch_index,
                phase_factor * coeff * branch_factor,
                idx,
                coeff * coefficient_factor,
            ));
        }

        let mut accum_ns = 0u128;
        for _ in 0..n_runs {
            let mut map: HashMap<u128, Complex<f64>> = HashMap::default();
            map.reserve(2 * items.len());
            let t0 = Instant::now();
            for &(branch_idx, branch_coeff, idx, nonbranch_coeff) in &precomputed {
                *map.entry(branch_idx).or_insert(Complex::zero()) += branch_coeff;
                *map.entry(idx).or_insert(Complex::zero()) += nonbranch_coeff;
            }
            accum_ns += t0.elapsed().as_nanos();
            std::hint::black_box(&map);
        }

        // Phase 3: combined (original pattern)
        let mut total_ns = 0u128;
        for _ in 0..n_runs {
            let mut map: HashMap<u128, Complex<f64>> = HashMap::default();
            let t0 = Instant::now();
            for &(coeff, idx) in &items {
                let branch_index = idx ^ stab_anticomm_bits;
                let symplectic = (destab_anticomm_bits & idx).count_ones();
                let mut phase = (2 * symplectic as u8) % 4;
                let active = idx & stab_anticomm_bits;
                let parity = (active & odd_phase_mask).count_ones() % 2;
                phase = (phase + 2 * parity as u8) % 4;
                let branch_phase = (phase + phase_decomp) % 4;
                let phase_factor: Complex<f64> = phase_table[branch_phase as usize].into();
                let branch_coeff = phase_factor * coeff * branch_factor;
                let nonbranch_coeff = coeff * coefficient_factor;
                *map.entry(branch_index).or_insert(Complex::zero()) += branch_coeff;
                *map.entry(idx).or_insert(Complex::zero()) += nonbranch_coeff;
            }
            total_ns += t0.elapsed().as_nanos();
            std::hint::black_box(&map);
        }

        let compute_us = compute_ns as f64 / n_runs as f64 / 1000.0;
        let accum_us = accum_ns as f64 / n_runs as f64 / 1000.0;
        let total_us = total_ns as f64 / n_runs as f64 / 1000.0;
        println!(
            "{:>10}  {:>9.0}µ  {:>9.0}µ  {:>9.0}µ  {:>5.1}%",
            n,
            compute_us,
            accum_us,
            total_us,
            compute_us / total_us * 100.0,
        );
    }
}
