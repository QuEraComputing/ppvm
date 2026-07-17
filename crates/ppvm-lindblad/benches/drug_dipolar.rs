// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Per-step cost of the Kossakowski-form dissipator on a *molecular dipolar
//! relaxation* workload (the ZULF-NMR drug-FID application), which stresses
//! the path differently from the superradiance chain in `kossakowski.rs`:
//!
//! - operators are **2-local rank-2 tensors** with ~4 Pauli terms each (one
//!   channel per site pair, per spatial harmonic `m`), not single-site σ⁻;
//! - `K` is **block-diagonal** over the 5 spatial components `m ∈ −2..2`,
//!   each block a dense `P×P` Gram matrix over the `P = C(N,2)` pairs;
//! - the basis strings are dense/high-weight, so most candidate pairs hit
//!   the **both-sided sandwich** (12-product) path — the arm the 2026-07-16
//!   ledger flagged as remaining headroom.
//!
//! Both specs (eigenmode jumps vs Kossakowski) generate the identical action;
//! the benchmark measures representation cost on one full `pc_step`.

use criterion::{Criterion, criterion_group, criterion_main};
use num::Complex;
use ppvm_lindblad::{JumpInput, LindbladSpec, PcStepConfig, Word, parse_pauli_string};
use std::hint::black_box;

const B: usize = 4096;
const GROW_STEPS: usize = 3;
const DT: f64 = 1e-3;
const N_M: usize = 5; // spatial harmonics m = -2..2

fn pstr(n: usize, sites: &[(usize, char)]) -> String {
    let mut s = vec!['I'; n];
    for &(q, c) in sites {
        s[q] = c;
    }
    s.into_iter().collect()
}

/// Deterministic pseudo-random 3D unit direction + distance for pair (a,b),
/// so the model is reproducible without an RNG dependency in the bench.
fn hashf(mut x: u64) -> f64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    (x >> 11) as f64 / (1u64 << 53) as f64
}

/// Dipolar coupling `b` and unit vector for every pair, from placing spins on
/// a jittered chain (real molecules: `b ∝ 1/r³`, generic directions).
#[allow(clippy::type_complexity)]
fn geometry(n: usize) -> (Vec<(usize, usize)>, Vec<f64>, Vec<[f64; 3]>) {
    let pos: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            [
                i as f64 + 0.3 * hashf(i as u64 * 3 + 1),
                0.4 * hashf(i as u64 * 3 + 2),
                0.4 * hashf(i as u64 * 3 + 3),
            ]
        })
        .collect();
    let mut pairs = Vec::new();
    let mut bmag = Vec::new();
    let mut dir = Vec::new();
    for a in 0..n {
        for b in (a + 1)..n {
            let d = [
                pos[a][0] - pos[b][0],
                pos[a][1] - pos[b][1],
                pos[a][2] - pos[b][2],
            ];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            pairs.push((a, b));
            bmag.push(1.0 / (r * r * r));
            dir.push([d[0] / r, d[1] / r, d[2] / r]);
        }
    }
    (pairs, bmag, dir)
}

/// Real rank-2 spherical harmonics (up to normalization) of a unit vector,
/// ordered m = -2,-1,0,1,2 — the spatial factors that make `Γ` rank 5.
fn y2(u: &[f64; 3]) -> [f64; N_M] {
    let (x, y, z) = (u[0], u[1], u[2]);
    [
        x * y,
        y * z,
        (3.0 * z * z - 1.0) / 2.0,
        x * z,
        (x * x - y * y) / 2.0,
    ]
}

/// The 2-local rank-2 tensor operator on pair (a,b) for tensor component
/// `mt` — a representative 4-term Pauli lincomb matching the high-field
/// dressed-tensor forms (T^(2,±2): XX∓YY ± i(XY±YX), etc.). The exact
/// coefficients are immaterial to the cost profile; the term *count* and
/// 2-locality are what matter.
fn tensor_op(n: usize, a: usize, b: usize, mt: usize) -> Vec<(String, Complex<f64>)> {
    let (i, j) = (Complex::new(0.0, 1.0), Complex::new(1.0, 0.0));
    match mt {
        0 | 4 => {
            let s = if mt == 0 { i } else { -i }; // ±2 components
            vec![
                (pstr(n, &[(a, 'X'), (b, 'X')]), j),
                (pstr(n, &[(a, 'Y'), (b, 'Y')]), -j),
                (pstr(n, &[(a, 'X'), (b, 'Y')]), s),
                (pstr(n, &[(a, 'Y'), (b, 'X')]), s),
            ]
        }
        1 | 3 => {
            let s = if mt == 1 { i } else { -i }; // ±1 components
            vec![
                (pstr(n, &[(a, 'X'), (b, 'Z')]), j),
                (pstr(n, &[(a, 'Y'), (b, 'Z')]), s),
                (pstr(n, &[(a, 'Z'), (b, 'X')]), j),
                (pstr(n, &[(a, 'Z'), (b, 'Y')]), s),
            ]
        }
        _ => vec![
            // m = 0
            (pstr(n, &[(a, 'X'), (b, 'X')]), j),
            (pstr(n, &[(a, 'Y'), (b, 'Y')]), j),
            (pstr(n, &[(a, 'Z'), (b, 'Z')]), Complex::new(2.0, 0.0)),
        ],
    }
}

fn hamiltonian_terms(n: usize, pairs: &[(usize, usize)], bmag: &[f64]) -> Vec<(String, f64)> {
    let mut h = Vec::new();
    for (k, &(a, b)) in pairs.iter().enumerate() {
        let jc = 0.1 * bmag[k]; // scalar J-coupling, XX+YY
        h.push((pstr(n, &[(a, 'X'), (b, 'X')]), jc));
        h.push((pstr(n, &[(a, 'Y'), (b, 'Y')]), jc));
    }
    h
}

/// Kossakowski ops (`N_M` blocks of `P` pair tensors) and the block-diagonal
/// `K = blockdiag(Γ_m)`, `Γ_m[μν] = Σ_{m'} c_μ^{m'} c_ν^{m'}` with
/// `c_μ^{m'} = b_μ Y_2^{m'}(r̂_μ)` — a rank-5 Gram block, exactly as the drug
/// pickles decompose.
#[allow(clippy::type_complexity)]
fn kossakowski_model(
    n: usize,
) -> (Vec<Vec<(String, Complex<f64>)>>, Vec<Vec<Complex<f64>>>) {
    let (pairs, bmag, dir) = geometry(n);
    let p = pairs.len();
    let c: Vec<[f64; N_M]> = (0..p)
        .map(|k| {
            let y = y2(&dir[k]);
            std::array::from_fn(|mp| bmag[k] * y[mp])
        })
        .collect();
    let mut ops = Vec::with_capacity(N_M * p);
    for mt in 0..N_M {
        for &(a, b) in &pairs {
            ops.push(tensor_op(n, a, b, mt));
        }
    }
    let m_ops = N_M * p;
    let mut k = vec![vec![Complex::new(0.0, 0.0); m_ops]; m_ops];
    for mt in 0..N_M {
        let off = mt * p;
        for mu in 0..p {
            for nu in 0..p {
                let g: f64 = (0..N_M).map(|mp| c[mu][mp] * c[nu][mp]).sum();
                k[off + mu][off + nu] = Complex::new(g, 0.0);
            }
        }
    }
    (ops, k)
}

/// Eigenmode jumps of the block-diagonal `K` (the dense representation the
/// Kossakowski path replaces): per block, `L_ν = √γ_ν Σ_μ V_μν T_μ`.
fn eigenmode_jumps(
    ops: &[Vec<(String, Complex<f64>)>],
    k: &[Vec<Complex<f64>>],
) -> Vec<JumpInput> {
    let m_ops = ops.len();
    let p = m_ops / N_M;
    let mut jumps = Vec::new();
    for mt in 0..N_M {
        let off = mt * p;
        let block = nalgebra::DMatrix::from_fn(p, p, |a, b| k[off + a][off + b].re);
        let eig = nalgebra::SymmetricEigen::new(block);
        for nu in 0..p {
            let g = eig.eigenvalues[nu];
            if g < 1e-12 {
                continue;
            }
            let mut lin = Vec::new();
            for mu in 0..p {
                let v = eig.eigenvectors[(mu, nu)];
                if v.abs() > 1e-14 {
                    for (s, cc) in &ops[off + mu] {
                        lin.push((s.clone(), cc * Complex::new(v, 0.0)));
                    }
                }
            }
            jumps.push(JumpInput { lincomb: lin, rate: g });
        }
    }
    jumps
}

/// Initial observable: γ-weighted transverse magnetization Σ_i X_i (the coil
/// quadrature), a sparse single-site sum like the drug FID initial operator.
fn observable(n: usize) -> (Vec<Word>, Vec<f64>) {
    let mut basis = Vec::new();
    let mut coeffs = Vec::new();
    for a in 0..n {
        basis.push(parse_pauli_string(&pstr(n, &[(a, 'X')]), n).unwrap().0);
        coeffs.push(1.0);
    }
    (basis, coeffs)
}

fn bench_drug(c: &mut Criterion) {
    let mut group = c.benchmark_group("pc_step_drug_dipolar");
    group.sample_size(10);
    for n in [10usize, 20, 32] {
        let (pairs, bmag, _) = geometry(n);
        let h = hamiltonian_terms(n, &pairs, &bmag);
        let (ops, k) = kossakowski_model(n);

        let spec_eig = LindbladSpec::new(n, &h, &eigenmode_jumps(&ops, &k)).unwrap();
        let mut spec_koss = LindbladSpec::new(n, &h, &[]).unwrap();
        spec_koss.add_kossakowski(&ops, &k).unwrap();

        let cfg = PcStepConfig {
            max_basis: B,
            admit_basis: Some(3 * B),
            ..Default::default()
        };
        let (mut basis, mut coeffs) = observable(n);
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

criterion_group!(benches, bench_drug);
criterion_main!(benches);
