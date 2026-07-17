// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Fast single-config profiler for the Kossakowski-form dissipator on the
//! molecular dipolar-relaxation (ZULF drug-FID) workload — the A/B harness
//! for the 2026-07-17-drug-kossakowski autotune campaign.
//!
//! Prints per-phase `pc_step_timed` breakdown (median of N steps) for the
//! Kossakowski path only. The eigenmode representation is *not* measured here
//! (its O(N³) dense-jump blowup is the thing this path removes — see the
//! `drug_dipolar` criterion bench for the documented representation ratio).
//!
//! Usage: `cargo run --release --example drug_profile -- [N] [B] [STEPS]`

use num::Complex;
use ppvm_lindblad::{LindbladSpec, PcStepConfig, Word, parse_pauli_string};
use std::time::Instant;

const N_M: usize = 5;

fn pstr(n: usize, sites: &[(usize, char)]) -> String {
    let mut s = vec!['I'; n];
    for &(q, c) in sites {
        s[q] = c;
    }
    s.into_iter().collect()
}

fn hashf(mut x: u64) -> f64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    (x >> 11) as f64 / (1u64 << 53) as f64
}

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
    let (mut pairs, mut bmag, mut dir) = (Vec::new(), Vec::new(), Vec::new());
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

fn y2(u: &[f64; 3]) -> [f64; N_M] {
    let (x, y, z) = (u[0], u[1], u[2]);
    [x * y, y * z, (3.0 * z * z - 1.0) / 2.0, x * z, (x * x - y * y) / 2.0]
}

fn tensor_op(n: usize, a: usize, b: usize, mt: usize) -> Vec<(String, Complex<f64>)> {
    let (i, j) = (Complex::new(0.0, 1.0), Complex::new(1.0, 0.0));
    match mt {
        0 | 4 => {
            let s = if mt == 0 { i } else { -i };
            vec![
                (pstr(n, &[(a, 'X'), (b, 'X')]), j),
                (pstr(n, &[(a, 'Y'), (b, 'Y')]), -j),
                (pstr(n, &[(a, 'X'), (b, 'Y')]), s),
                (pstr(n, &[(a, 'Y'), (b, 'X')]), s),
            ]
        }
        1 | 3 => {
            let s = if mt == 1 { i } else { -i };
            vec![
                (pstr(n, &[(a, 'X'), (b, 'Z')]), j),
                (pstr(n, &[(a, 'Y'), (b, 'Z')]), s),
                (pstr(n, &[(a, 'Z'), (b, 'X')]), j),
                (pstr(n, &[(a, 'Z'), (b, 'Y')]), s),
            ]
        }
        _ => vec![
            (pstr(n, &[(a, 'X'), (b, 'X')]), j),
            (pstr(n, &[(a, 'Y'), (b, 'Y')]), j),
            (pstr(n, &[(a, 'Z'), (b, 'Z')]), Complex::new(2.0, 0.0)),
        ],
    }
}

#[allow(clippy::type_complexity)]
fn model(n: usize) -> (Vec<(String, f64)>, Vec<Vec<(String, Complex<f64>)>>, Vec<Vec<Complex<f64>>>) {
    let (pairs, bmag, dir) = geometry(n);
    let p = pairs.len();
    let mut h = Vec::new();
    for (k, &(a, b)) in pairs.iter().enumerate() {
        let jc = 0.1 * bmag[k];
        h.push((pstr(n, &[(a, 'X'), (b, 'X')]), jc));
        h.push((pstr(n, &[(a, 'Y'), (b, 'Y')]), jc));
    }
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
    (h, ops, k)
}

fn observable(n: usize) -> (Vec<Word>, Vec<f64>) {
    let mut basis = Vec::new();
    let mut coeffs = Vec::new();
    for a in 0..n {
        basis.push(parse_pauli_string(&pstr(n, &[(a, 'X')]), n).unwrap().0);
        coeffs.push(1.0);
    }
    (basis, coeffs)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(20);
    let b: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4096);
    let steps: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(8);
    let dt = 1e-3;

    let (h, ops, k) = model(n);
    let mut spec = LindbladSpec::new(n, &h, &[]).unwrap();
    let t0 = Instant::now();
    spec.add_kossakowski(&ops, &k).unwrap();
    let build_ms = t0.elapsed().as_secs_f64() * 1e3;

    let cfg = PcStepConfig { max_basis: b, admit_basis: Some(3 * b), ..Default::default() };
    let (mut basis, mut coeffs) = observable(n);
    // Grow into a realistic capped basis (not timed).
    for _ in 0..3 {
        spec.pc_step(&mut basis, &mut coeffs, dt, &[], &cfg).unwrap();
    }

    let mut totals = Vec::new();
    let (mut l1, mut e1, mut x1, mut l2, mut e2, mut x2) = (0u64, 0u64, 0u64, 0u64, 0u64, 0u64);
    for _ in 0..steps {
        let mut bb = basis.clone();
        let mut cf = coeffs.clone();
        let t = spec.pc_step_timed(&mut bb, &mut cf, dt, &[], &cfg).unwrap();
        totals.push(t.total_us());
        l1 += t.leakage1_us;
        e1 += t.expand1_us;
        x1 += t.expm1_us;
        l2 += t.leakage2_us;
        e2 += t.expand2_us;
        x2 += t.expm2_us;
    }
    totals.sort_unstable();
    let med = totals[totals.len() / 2] as f64 / 1e3;
    let s = steps as f64;
    println!(
        "N={n} B={b} pairs={} ops={} nnz(K)={}",
        n * (n - 1) / 2,
        ops.len(),
        N_M * (n * (n - 1) / 2) * (n * (n - 1) / 2)
    );
    println!("  add_kossakowski build: {build_ms:.0} ms");
    println!("  median total/step: {med:.1} ms  (over {steps} steps)");
    println!(
        "  phase avg (ms): leak1 {:.1}  expm1 {:.1}  leak2 {:.1}  expm2 {:.1}  expand {:.1}",
        l1 as f64 / s / 1e3,
        x1 as f64 / s / 1e3,
        l2 as f64 / s / 1e3,
        x2 as f64 / s / 1e3,
        (e1 + e2) as f64 / s / 1e3,
    );
    let diss = (l1 + l2) as f64;
    let tot = (l1 + e1 + x1 + l2 + e2 + x2) as f64;
    println!("  leakage(action) share: {:.0}%", 100.0 * diss / tot);
}
