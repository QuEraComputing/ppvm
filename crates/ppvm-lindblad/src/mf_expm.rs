// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Matrix-free `exp(dt · L*) · b` for the real (`f64`) path, driven by the
//! external `quspin-expm` crate.
//!
//! Instead of materialising the in-basis-restricted generator as a CSR, the
//! per-column generator action is computed ONCE per expm call (via
//! [`build_mf_cols`]) and reused, CSC-style, across every Krylov/Taylor matvec
//! by [`CachedCscOp`] (a [`quspin_types::LinearOperator`]) fed to
//! [`quspin_expm::ExpmOp::from_parts`]. Caching the action is the dominant
//! win: the old fully-matrix-free op recomputed the Pauli-commutator action on
//! every column on every matvec (~15 matvecs/call), which was ~90% of the
//! per-step cost; building it once collapses each matvec to a cheap CSC
//! scatter. Using `from_parts` (rather than `ExpmOp::new`) supplies the
//! diagonal shift `μ`, the partition count `s`, and the truncation order `m*`
//! directly, which BYPASSES quspin's adaptive parameter selection — so the
//! 1-norm *estimator* and `dot_transpose` are never invoked on the
//! single-vector `apply` path; only [`LinearOperator::dot`] runs.
//!
//! `μ`, the trace, and the exact column 1-norm of `A − μ·I` are computed in
//! the same single action pass as the cache. The `(m, s)` Taylor partition is
//! picked with the tolerance-matched tables in [`crate::expm`]: a relaxed
//! `tol=1e-6` table when the PC prunes coarsely (`drop_tol ≥ 1e-4`), else the
//! double-precision table (keeping the exact-reference test paths bit-exact).

use crate::{LindbladSpec, Word, build_basis_index, expm};
use fxhash::{FxBuildHasher, FxHashMap};
use num::Complex;
use quspin_types::{ExpmComputation, LinearOperator, QuSpinError};
use rayon::prelude::*;

/// Per-column in-basis action of the real generator `M`, plus the data the
/// `(m, s)`/`μ` selection needs — all from ONE action pass over the basis.
///
/// Returns `(cols, per_col)` where `cols[c]` holds `(row, coeff)` for every
/// action output of `L*(basis[c])` that lands back in `basis` (CSC column
/// `c`), and `per_col[c] = (raw, diag)` with `raw = Σ|coeff|` over ALL action
/// outputs (in- and out-of-basis, matching the historical 1-norm) and `diag`
/// the coefficient of the term whose output Word equals the input Word. The
/// action is computed once here and reused by [`CachedCscOp`] across every
/// Krylov/Taylor matvec, so it is never recomputed per matvec (the dominant
/// cost of the old matrix-free op). Same reusable scratch pattern as
/// [`crate::orbit_rep::build_orbit_rep_cols`]; numerics identical.
fn build_mf_cols(
    spec: &LindbladSpec,
    basis: &[Word],
    index: &FxHashMap<Word, u32>,
) -> (Vec<Vec<(u32, f64)>>, Vec<(f64, f64)>) {
    basis
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
            |(s1, s2, lm), p| {
                let terms = spec.compute_action_terms(p, s1, s2, lm);
                let mut out = Vec::with_capacity(terms.len());
                let mut raw = 0.0;
                let mut diag = 0.0;
                for (w, c) in terms.iter() {
                    raw += c.abs();
                    if w == p {
                        diag = *c;
                    }
                    if let Some(&row) = index.get(w) {
                        out.push((row, *c));
                    }
                }
                (out, (raw, diag))
            },
        )
        .unzip()
}

/// Borrowed CSC-style view of the in-basis-restricted real generator `M`,
/// backed by the cached per-column action [`build_mf_cols`] (computed once
/// per expm call). `dot` does a CSC matvec `y = M·x` against the cache with
/// no per-matvec action recompute and no CSR materialisation. Mirrors
/// [`crate::orbit_rep::OrbitRepCscOp`] for the real `f64` path.
///
/// Borrowed, not owned: `quspin-types` provides a blanket `LinearOperator`
/// impl for `&T`, so `ExpmOp::from_parts(op, …)` accepts a `CachedCscOp` by
/// value while it keeps borrowing `cols`.
struct CachedCscOp<'a> {
    cols: &'a [Vec<(u32, f64)>],
    dim: usize,
}

impl LinearOperator<f64> for CachedCscOp<'_> {
    fn dim(&self) -> usize {
        self.dim
    }

    fn parallel_hint(&self) -> bool {
        // `dot` parallelises internally over column chunks, and we drive the
        // sequential single-vector `apply` path; never let quspin run its
        // persistent-thread pool on top of our rayon parallelism.
        false
    }

    fn dot(&self, overwrite: bool, input: &[f64], output: &mut [f64]) -> Result<(), QuSpinError> {
        let n = self.dim;
        if n == 0 {
            return Ok(());
        }
        let num_threads = rayon::current_num_threads().max(1);
        let chunk_size = n.div_ceil(num_threads);

        // Parallelise over column chunks; each thread accumulates into a dense
        // local `y` of length `dim`, reading the cached action.
        let partial_ys: Vec<Vec<f64>> = self
            .cols
            .par_chunks(chunk_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| {
                let c_offset = chunk_idx * chunk_size;
                let mut y_local = vec![0.0; n];
                for (c_local, col) in chunk.iter().enumerate() {
                    let xc = input[c_offset + c_local];
                    if xc == 0.0 {
                        continue;
                    }
                    for &(row, val) in col.iter() {
                        y_local[row as usize] += val * xc;
                    }
                }
                y_local
            })
            .collect();

        if overwrite {
            output.fill(0.0);
        }
        for partial in &partial_ys {
            for (oi, &pi) in output.iter_mut().zip(partial.iter()) {
                *oi += pi;
            }
        }
        Ok(())
    }

    fn trace(&self) -> f64 {
        // Computed eagerly in `expm_apply_mf`; never reached on the
        // `from_parts` + single-vector `apply` path.
        unreachable!("CachedCscOp::trace not used on the from_parts apply path")
    }

    fn onenorm(&self, _shift: f64) -> f64 {
        unreachable!("CachedCscOp::onenorm not used on the from_parts apply path")
    }

    fn dot_transpose(
        &self,
        _overwrite: bool,
        _input: &[f64],
        _output: &mut [f64],
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CachedCscOp: dot_transpose not used on the from_parts apply path".into(),
        ))
    }

    fn dot_many(
        &self,
        _overwrite: bool,
        _input: ndarray::ArrayView2<'_, f64>,
        _output: ndarray::ArrayViewMut2<'_, f64>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CachedCscOp: dot_many not used on the from_parts apply path".into(),
        ))
    }

    fn dot_chunk(
        &self,
        _overwrite: bool,
        _input: &[f64],
        _output_chunk: &mut [f64],
        _row_start: usize,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CachedCscOp: dot_chunk not used on the from_parts apply path".into(),
        ))
    }

    fn dot_transpose_chunk(
        &self,
        _input: &[f64],
        _output: &[<f64 as ExpmComputation>::Atomic],
        _rows: std::ops::Range<usize>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CachedCscOp: dot_transpose_chunk not used on the from_parts apply path".into(),
        ))
    }
}

/// Compute `exp(dt · M) · coeffs` for the in-basis-restricted generator
/// `M`, matrix-free, via `quspin-expm`. Returns a fresh `Vec<f64>` of length
/// `basis.len()`.
///
/// One matrix-free pass extracts the diagonal shift `μ = tr(M)/n` and the
/// exact column 1-norm of `M − μ·I`; from `‖dt·(M−μI)‖₁` we pick the Taylor
/// partition `(m*, s)` and hand everything to
/// [`quspin_expm::ExpmOp::from_parts`].
pub(crate) fn expm_apply_mf(
    spec: &LindbladSpec,
    basis: &[Word],
    dt: f64,
    coeffs: &[f64],
    drop_tol: f64,
) -> Vec<f64> {
    let n = basis.len();
    if n == 0 {
        return Vec::new();
    }

    // ONE action pass: build the CSC cache `cols` (reused across every matvec)
    // and, in the same pass, `per_col = (raw, diag)` for the `μ`/1-norm
    // selection. `raw = Σ|coeff|` (all outputs), `diag` = coeff of the
    // output == input term. From these: `trace = Σ diag`, `μ = trace/n`, and
    // the column 1-norm of `M − μ·I` is `raw − |diag| + |diag − μ|`.
    let index = build_basis_index(basis);
    let (cols, per_col) = build_mf_cols(spec, basis, &index);

    let trace: f64 = per_col.iter().map(|(_, d)| *d).sum();
    let mu = trace / n as f64;
    let onenorm = per_col
        .iter()
        .map(|(raw, diag)| raw - diag.abs() + (diag - mu).abs())
        .fold(0.0_f64, f64::max);

    // Pick the Taylor backward-error tolerance to match the basis truncation:
    // when the PC prunes coarsely (drop_tol >= 1e-4) a double-precision exp is
    // ~10 orders more accurate than the state it acts on, so the relaxed
    // (tol=1e-6, still >=100x tighter than the cut) table is used — it admits a
    // lower-degree Taylor polynomial and cuts the SpMV count with no effect on
    // the truncated result. At tight/zero drop_tol we keep double precision so
    // the exact-reference paths (orbit-rep / merged) still agree bit-for-bit.
    let t_norm = dt.abs() * onenorm;
    let (m_star, s, expm_tol) = if drop_tol >= 1e-4 {
        let (m, s) = expm::select_ms_loose(t_norm, None);
        (m, s, 1e-6_f64)
    } else {
        let (m, s) = expm::select_ms(t_norm, None);
        (m, s, 1e-12_f64)
    };

    let op = CachedCscOp { cols: &cols, dim: n };
    let expm = quspin_expm::ExpmOp::from_parts(op, dt, mu, s as usize, m_star as usize, expm_tol);

    let mut v = coeffs.to_vec();
    expm.apply(ndarray::ArrayViewMut1::from(v.as_mut_slice()))
        .expect("expm apply");
    v
}

/// `exp(dt · M) · b` where `M` is the REAL in-basis-restricted generator but
/// the input vector `b` is complex. Because `M` is real,
/// `exp(dt·M)·(re + i·im) = exp(dt·M)·re + i·exp(dt·M)·im`, so we split the
/// complex vector into its real and imaginary parts, run two real
/// matrix-free applies, and recombine. Fully matrix-free; no CSR.
pub(crate) fn expm_apply_mf_cxvec(
    spec: &LindbladSpec,
    basis: &[Word],
    dt: f64,
    b: &[Complex<f64>],
    drop_tol: f64,
) -> Vec<Complex<f64>> {
    let n = basis.len();
    if n == 0 {
        return Vec::new();
    }
    let re: Vec<f64> = b.iter().map(|z| z.re).collect();
    let im: Vec<f64> = b.iter().map(|z| z.im).collect();
    let re_out = expm_apply_mf(spec, basis, dt, &re, drop_tol);
    let im_out = expm_apply_mf(spec, basis, dt, &im, drop_tol);
    re_out
        .into_iter()
        .zip(im_out)
        .map(|(r, i)| Complex::new(r, i))
        .collect()
}


