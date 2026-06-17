// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Per-step orbit-representative evolution under translation symmetry.
//!
//! The state lives entirely in **orbit-rep form** throughout: `basis`
//! contains only canonical translation-orbit representatives, and
//! `coeffs` are complex (one per rep). The dynamics `L*` is computed
//! with **phase-aware action** — for each output Pauli `q`, we
//! canonicalize `q` to its orbit rep `r_q` with shift counter `cnt_q`,
//! and accumulate `χ_k(g_{cnt_q}) · v · c_r` (where `v` is the matrix
//! element of `L*` between input rep `r` and output `q`).
//!
//! This gives the user the per-step memory reduction promised by the
//! symmetry-merging paper: basis is ~|G|× smaller than the full-basis
//! representation, throughout the entire evolution.
//!
//! The CSR matrix that the Krylov-expm machinery operates on has
//! **complex entries** (because of the `χ_k(g)` phase factors), so the
//! complex CSR-backed `mf_expm::expm_apply_csr_cx` (driving the external
//! `quspin-expm` engine) is used instead of the real path.
//!
//! ## Limitations
//!
//! - Caller is responsible for ensuring the input basis is in orbit-rep
//!   form (i.e. each entry is the canonical representative of its
//!   translation orbit). Use [`canonicalize_basis_to_rep`] if needed.
//! - The momentum sector `k_modes` is fixed for the duration of one
//!   pc_step call. To compute a full site-resolved profile, call
//!   `pc_step_orbit_rep` once per momentum mode and inverse-Fourier
//!   the results.

use crate::expm::CsrCx;
use crate::{Error, LindbladSpec};
use fxhash::{FxBuildHasher, FxHashMap};
use num::Complex;
use ppvm_runtime::symmetry::TranslationGroup;
use rayon::prelude::*;
use sprs::CsMatI;

// Word type re-exported from lib.rs.
use crate::Word;

/// Replace each entry of `basis` with its canonical orbit
/// representative under `group`. Pure rewrite; coefficients are
/// untouched. Useful to enforce the orbit-rep invariant before calling
/// [`pc_step_orbit_rep`].
///
/// Does NOT deduplicate — if multiple input entries collapse to the
/// same rep, both are kept (caller should run a merge afterwards).
pub fn canonicalize_basis_to_rep(basis: &mut [Word], group: &TranslationGroup) {
    for w in basis.iter_mut() {
        let canon = group.canonicalize(w);
        *w = canon;
    }
}

/// Build the in-basis-restricted complex CSR of the L* generator,
/// **in orbit-rep form** at momentum sector `k_modes`.
///
/// For each input rep `r` (column index) and each output Pauli `q` of
/// `L*(r) = Σ_q v_q · q`:
/// 1. Canonicalize `q` → `(r_q, cnt_q)`.
/// 2. If `r_q` is in `basis` at row index `i`, add `χ_k(g_{cnt_q}) · v_q`
///    to `M[i, col_r]`.
/// 3. Outputs not in `basis` are silently dropped (they would be
///    handled by [`leakage_orbit_rep`] in the surrounding pc_step).
///
/// Returns a `(basis.len() × basis.len())` complex CSR.
pub fn generator_csr_orbit_rep(
    spec: &LindbladSpec,
    basis: &[Word],
    group: &TranslationGroup,
    k_modes: &[i32],
) -> CsrCx {
    let n = basis.len();
    if n == 0 {
        return CsMatI::new_from_unsorted((0, 0), vec![0], Vec::new(), Vec::new())
            .map_err(|(_, _, _, e)| e)
            .expect("empty CSR");
    }

    // Build a Word → row index map.
    let index = build_basis_index(basis);

    // Per-column accumulation: for each rep r, list of (row, phase × v).
    let cols: Vec<Vec<(u32, Complex<f64>)>> = basis
        .par_iter()
        .map_init(
            || {
                (
                    Vec::<u32>::with_capacity(spec.n_qubits()),
                    Vec::<u32>::with_capacity(128),
                    FxHashMap::<Word, Complex<f64>>::with_capacity_and_hasher(
                        128,
                        FxBuildHasher::default(),
                    ),
                )
            },
            |(s1, s2, lm), r| {
                let terms = spec.compute_action_terms(r, s1, s2, lm);
                let mut out = Vec::with_capacity(terms.len());
                for (q, v) in terms.iter() {
                    let (r_q, cnt_q) = group.canonicalize_with_shift(q);
                    if let Some(&row) = index.get(&r_q) {
                        let phase = group.character(k_modes, &cnt_q);
                        out.push((row, phase * (*v as f64)));
                    }
                }
                out
            },
        )
        .collect();

    // Build per-row (col, value) lists. Multiple (col, row) entries
    // for the same (row, col) cell get summed (sprs's `new_from_unsorted`
    // doesn't dedup, and would also choke on unsorted within-row
    // indices); we do both up front.
    let mut row_data: Vec<Vec<(u32, Complex<f64>)>> = vec![Vec::new(); n];
    for (col, col_data) in cols.iter().enumerate() {
        for &(row, v) in col_data.iter() {
            row_data[row as usize].push((col as u32, v));
        }
    }
    let mut row_ptr = vec![0usize; n + 1];
    let mut indices: Vec<u32> = Vec::new();
    let mut data: Vec<Complex<f64>> = Vec::new();
    for (i, mut rd) in row_data.into_iter().enumerate() {
        // sort by col, then sum duplicates.
        rd.sort_by_key(|&(c, _)| c);
        let mut last_col: Option<u32> = None;
        for (c, v) in rd {
            if last_col == Some(c) {
                let n_so_far = data.len();
                data[n_so_far - 1] += v;
            } else {
                indices.push(c);
                data.push(v);
                last_col = Some(c);
            }
        }
        row_ptr[i + 1] = indices.len();
    }

    CsMatI::new((n, n), row_ptr, indices, data)
}

/// Phase-aware leakage: out-of-basis component of `L*(O_k)` where `O_k`
/// is the operator represented by `basis` (orbit reps) and `coeffs`
/// (complex coefficients in momentum sector `k_modes`).
///
/// For each input rep `r` with coefficient `c_r`, and each output `q`
/// of `L*(r) = Σ_q v_q · q`:
/// 1. Canonicalize `q` → `(r_q, cnt_q)`.
/// 2. If `r_q` NOT in `basis` and NOT in `protected`:
///    `merged[r_q] += χ_k(g_{cnt_q}) · v_q · c_r`.
///
/// Returns `(r_q, sum)` pairs for all candidates with nonzero sum.
pub fn leakage_orbit_rep(
    spec: &LindbladSpec,
    basis: &[Word],
    coeffs: &[Complex<f64>],
    protected: &[Word],
    group: &TranslationGroup,
    k_modes: &[i32],
) -> Result<Vec<(Word, Complex<f64>)>, Error> {
    if basis.len() != coeffs.len() {
        return Err(Error::LengthMismatch {
            what: "basis and coeffs",
            a: basis.len(),
            b: coeffs.len(),
        });
    }
    let in_basis: FxHashMap<&Word, ()> = basis.iter().map(|w| (w, ())).collect();
    let protected_set: FxHashMap<&Word, ()> = protected.iter().map(|w| (w, ())).collect();

    const CHUNK_SIZE: usize = 4096;
    let mut merged: FxHashMap<Word, Complex<f64>> = FxHashMap::default();
    for chunk_start in (0..basis.len()).step_by(CHUNK_SIZE) {
        let chunk_end = (chunk_start + CHUNK_SIZE).min(basis.len());
        let chunk_basis = &basis[chunk_start..chunk_end];
        let chunk_coeffs = &coeffs[chunk_start..chunk_end];
        let local: Vec<Vec<(Word, Complex<f64>)>> = chunk_basis
            .par_iter()
            .zip(chunk_coeffs.par_iter())
            .map_init(
                || {
                    (
                        Vec::<u32>::with_capacity(spec.n_qubits()),
                        Vec::<u32>::with_capacity(128),
                        FxHashMap::<Word, Complex<f64>>::with_capacity_and_hasher(
                            128,
                            FxBuildHasher::default(),
                        ),
                    )
                },
                |(s1, s2, lm), (r, &c_r)| {
                    let terms = spec.compute_action_terms(r, s1, s2, lm);
                    let mut out = Vec::with_capacity(terms.len());
                    for (q, v) in terms.iter() {
                        let (r_q, cnt_q) = group.canonicalize_with_shift(q);
                        if !in_basis.contains_key(&r_q) && !protected_set.contains_key(&r_q) {
                            let phase = group.character(k_modes, &cnt_q);
                            out.push((r_q, phase * (*v as f64) * c_r));
                        }
                    }
                    out
                },
            )
            .collect();
        for v in local {
            for (k, val) in v {
                *merged.entry(k).or_insert(Complex::new(0.0, 0.0)) += val;
            }
        }
    }
    Ok(merged.into_iter().filter(|(_, c)| c.norm() > 0.0).collect())
}

/// Per-step orbit-rep predictor-corrector evolution.
///
/// All state lives in orbit-rep form throughout. Each pc step does:
/// 1. Phase-aware leakage from `(basis, coeffs)`; append leakage reps
///    with `|c| > tau_add`.
/// 2. Predictor: build complex orbit-rep CSR, run `expm_multiply_cx`.
/// 3. Phase-aware leakage from predicted state; append further reps.
/// 4. Corrector: rebuild CSR (basis grew), `expm_multiply_cx` from
///    pre-step coefficients (with the same trick as `pc_step_inner`).
/// 5. Prune `|c| < drop_tol`, protected reps never dropped.
///
/// `basis` is assumed to contain only canonical orbit representatives.
/// If not, [`canonicalize_basis_to_rep`] should be called first.
#[allow(clippy::too_many_arguments)]
pub fn pc_step_orbit_rep(
    spec: &LindbladSpec,
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<Complex<f64>>,
    dt: f64,
    tau_add: f64,
    drop_tol: f64,
    protected: &[Word],
    group: &TranslationGroup,
    k_modes: &[i32],
) -> Result<(), Error> {
    // 1. First-hop phase-aware leakage.
    let leak = leakage_orbit_rep(spec, basis, coeffs, protected, group, k_modes)?;
    for (w, v) in leak {
        if v.norm() > tau_add {
            basis.push(w);
            coeffs.push(Complex::new(0.0, 0.0));
        }
    }
    // 2. Predictor: build CSR + expm.
    let csr = generator_csr_orbit_rep(spec, basis, group, k_modes);
    let coeffs_predict = crate::mf_expm::expm_apply_csr_cx(&csr, dt, coeffs);
    drop(csr);
    // 3. Second-hop leakage from predicted state.
    let leak2 = leakage_orbit_rep(spec, basis, &coeffs_predict, protected, group, k_modes)?;
    drop(coeffs_predict);
    for (w, v) in leak2 {
        if v.norm() > tau_add {
            basis.push(w);
            coeffs.push(Complex::new(0.0, 0.0));
        }
    }
    // 4. Corrector: rebuild + expm from pre-step state.
    let csr = generator_csr_orbit_rep(spec, basis, group, k_modes);
    *coeffs = crate::mf_expm::expm_apply_csr_cx(&csr, dt, coeffs);
    drop(csr);
    // 5. Prune.
    if drop_tol > 0.0 {
        prune_basis_complex_local(basis, coeffs, drop_tol, protected);
    }
    Ok(())
}

fn prune_basis_complex_local(
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<Complex<f64>>,
    drop_tol: f64,
    protected: &[Word],
) {
    if drop_tol <= 0.0 {
        return;
    }
    let protected_set: fxhash::FxHashSet<&Word> = protected.iter().collect();
    let mut write = 0;
    for read in 0..basis.len() {
        if coeffs[read].norm() >= drop_tol || protected_set.contains(&basis[read]) {
            if write != read {
                basis.swap(write, read);
                coeffs.swap(write, read);
            }
            write += 1;
        }
    }
    basis.truncate(write);
    coeffs.truncate(write);
}

fn build_basis_index(basis: &[Word]) -> FxHashMap<Word, u32> {
    let mut index: FxHashMap<Word, u32> = FxHashMap::default();
    for (i, w) in basis.iter().enumerate() {
        let prev = index.insert(w.clone(), i as u32);
        debug_assert!(
            prev.is_none(),
            "orbit-rep basis contains duplicate at positions {} and {}",
            prev.unwrap(),
            i,
        );
    }
    index
}
