// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Per-step cost of the Kossakowski-form dissipator vs the equivalent
//! eigenmode-jump representation, on the subwavelength superradiance chain
//! (free-space photon-mediated collective σ⁻ decay, d = 0.1 λ₀).
//!
//! Both specs generate the identical adjoint action; the difference is
//! representation cost: eigenmode jumps pay `N · (2N)²` Pauli products per
//! dissipator evaluation, Kossakowski pairs pay `4·nnz(Γ) = 4N²`.
//!
//! The benchmark grows a realistic working basis with a few capped
//! `pc_step` calls, then measures one full `pc_step` (two leakage passes +
//! predictor/corrector expm) from a cloned copy of that basis.

use criterion::{Criterion, criterion_group, criterion_main};
use num::Complex;
use ppvm_lindblad::{JumpInput, LindbladSpec, PcStepConfig, Word, parse_pauli_string};
use std::f64::consts::PI;
use std::hint::black_box;

const G0: f64 = 1.0;
const D_OVER_LAM: f64 = 0.1;
const B: usize = 4096;
const GROW_STEPS: usize = 3;
const DT: f64 = 0.01;

/// Free-space couplings `(J, Γ)` of a chain along x with spacing `d·λ₀`,
/// circular polarization `(1, i, 0)/√2`.
fn chain_couplings(n: usize) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let k0 = 2.0 * PI;
    let mut j = vec![vec![0.0; n]; n];
    let mut gam = vec![vec![0.0; n]; n];
    #[allow(clippy::needless_range_loop)]
    for a in 0..n {
        for b in 0..n {
            if a == b {
                gam[a][b] = G0;
                continue;
            }
            let r = (a as f64 - b as f64).abs() * D_OVER_LAM;
            let kr = k0 * r;
            let e = Complex::from_polar(1.0, kr);
            let pref = e / (4.0 * PI * k0 * k0 * r * r * r);
            // p†·G·p with p = (1, i, 0)/√2 and r̂ = x̂:
            // p†·(kr²+ikr−1)·1·p = (kr²+ikr−1); p†·r̂r̂·p = 1/2.
            let g = pref
                * (Complex::new(kr * kr - 1.0, kr) - Complex::new(kr * kr - 3.0, 3.0 * kr) * 0.5);
            j[a][b] = -3.0 * PI * G0 / k0 * g.re;
            gam[a][b] = 6.0 * PI * G0 / k0 * g.im;
        }
    }
    (j, gam)
}

fn pstr(n: usize, sites: &[(usize, char)]) -> String {
    let mut s = vec!['I'; n];
    for &(q, c) in sites {
        s[q] = c;
    }
    s.into_iter().collect()
}

#[allow(clippy::needless_range_loop)]
fn hamiltonian_terms(n: usize, j: &[Vec<f64>]) -> Vec<(String, f64)> {
    let mut h = Vec::new();
    for a in 0..n {
        for b in (a + 1)..n {
            if j[a][b].abs() > 1e-14 {
                h.push((pstr(n, &[(a, 'X'), (b, 'X')]), j[a][b] / 2.0));
                h.push((pstr(n, &[(a, 'Y'), (b, 'Y')]), j[a][b] / 2.0));
            }
        }
    }
    h
}

fn sigma_minus(site: usize, n: usize) -> Vec<(String, Complex<f64>)> {
    vec![
        (pstr(n, &[(site, 'X')]), Complex::new(0.5, 0.0)),
        (pstr(n, &[(site, 'Y')]), Complex::new(0.0, -0.5)),
    ]
}

/// Eigenmode jumps `L_ν = √γ_ν Σ_j V_jν σ⁻_j` from `Γ = V diag(γ) Vᵀ`.
fn eigenmode_jumps(n: usize, gam: &[Vec<f64>]) -> Vec<JumpInput> {
    let mat = nalgebra::DMatrix::from_fn(n, n, |a, b| gam[a][b]);
    let eig = nalgebra::SymmetricEigen::new(mat);
    let mut jumps = Vec::new();
    for nu in 0..n {
        let g = eig.eigenvalues[nu];
        if g < 1e-12 {
            continue;
        }
        let mut lin = Vec::new();
        for j in 0..n {
            let v = eig.eigenvectors[(j, nu)];
            if v.abs() > 1e-14 {
                for (p, c) in sigma_minus(j, n) {
                    lin.push((p, c * v));
                }
            }
        }
        jumps.push(JumpInput {
            lincomb: lin,
            rate: g,
        });
    }
    jumps
}

/// `O = Σ_nm Γ_nm σ⁺_n σ⁻_m` as a real Pauli sum.
#[allow(clippy::needless_range_loop)]
fn observable(n: usize, gam: &[Vec<f64>]) -> (Vec<Word>, Vec<f64>) {
    let mut basis = Vec::new();
    let mut coeffs = Vec::new();
    let mut push = |s: String, c: f64| {
        basis.push(parse_pauli_string(&s, n).unwrap().0);
        coeffs.push(c);
    };
    push(pstr(n, &[]), n as f64 * G0 / 2.0);
    for a in 0..n {
        push(pstr(n, &[(a, 'Z')]), G0 / 2.0);
        for b in (a + 1)..n {
            push(pstr(n, &[(a, 'X'), (b, 'X')]), gam[a][b] / 2.0);
            push(pstr(n, &[(a, 'Y'), (b, 'Y')]), gam[a][b] / 2.0);
        }
    }
    (basis, coeffs)
}

fn bench_kossakowski(c: &mut Criterion) {
    let mut group = c.benchmark_group("pc_step_superradiance");
    group.sample_size(10);
    for n in [10usize, 20, 30] {
        let (j, gam) = chain_couplings(n);
        let h = hamiltonian_terms(n, &j);

        let spec_eig = LindbladSpec::new(n, &h, &eigenmode_jumps(n, &gam)).unwrap();
        let mut spec_koss = LindbladSpec::new(n, &h, &[]).unwrap();
        let ops: Vec<_> = (0..n).map(|q| sigma_minus(q, n)).collect();
        let k: Vec<Vec<Complex<f64>>> = gam
            .iter()
            .map(|row| row.iter().map(|&v| Complex::new(v, 0.0)).collect())
            .collect();
        spec_koss.add_kossakowski(&ops, &k).unwrap();

        // Grow a realistic capped working basis once (shared by both).
        let cfg = PcStepConfig {
            max_basis: B,
            admit_basis: Some(3 * B),
            ..Default::default()
        };
        let (mut basis, mut coeffs) = observable(n, &gam);
        for _ in 0..GROW_STEPS {
            spec_koss
                .pc_step(&mut basis, &mut coeffs, DT, &[], &cfg)
                .unwrap();
        }

        for (label, spec) in [("eigenmode", &spec_eig), ("kossakowski", &spec_koss)] {
            group.bench_function(format!("{label}_n{n}"), |bch| {
                bch.iter(|| {
                    let mut b = basis.clone();
                    let mut cf = coeffs.clone();
                    spec.pc_step(&mut b, &mut cf, DT, &[], &cfg).unwrap();
                    black_box(cf.len())
                })
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_kossakowski);
criterion_main!(benches);
