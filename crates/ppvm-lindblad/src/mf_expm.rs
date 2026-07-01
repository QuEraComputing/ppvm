// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Matrix-free `exp(dt · L*) · b` for the real (`f64`) path, driven by the
//! external `quspin-expm` crate.
//!
//! Instead of materialising the in-basis-restricted generator as a CSR, the
//! per-column generator action is computed ONCE per expm call (via
//! [`build_mf_cols`]) and reused, CSC-style, across every Krylov/Taylor
//! matvec by [`CachedCscOp`] (a [`quspin_types::LinearOperator`]) which is
//! fed to [`quspin_expm::ExpmOp::from_parts`]. Using `from_parts` (rather
//! than `ExpmOp::new`) supplies the diagonal shift `μ`, the partition count
//! `s`, and the truncation order `m*` directly, which BYPASSES quspin's
//! adaptive parameter selection. As a consequence the 1-norm *estimator* and
//! `dot_transpose` are never invoked on the single-vector `apply` path — only
//! [`LinearOperator::dot`] runs, inside the Taylor/Horner loop.
//!
//! The shift `μ` and the exact column 1-norm of `A − μ·I` are computed in one
//! pass over the cached columns (the same `(m, s)` selection table as the
//! retired hand-rolled engine in [`crate::expm`]).

use crate::{LindbladSpec, Word, build_basis_index, expm};
use fxhash::{FxBuildHasher, FxHashMap};
use num::Complex;
use quspin_types::{ExpmComputation, LinearOperator, QuSpinError};
use rayon::prelude::*;

/// Build the per-column action of the in-basis-restricted real generator
/// `M` (the same matrix [`LindbladSpec::generator_csr`] would build, never
/// materialised).
///
/// Returns, for each input Pauli `basis[c]` (column `c`), the list of
/// `(row, coeff)` pairs for every action output term `(w, coeff)` of
/// `L*(basis[c])` whose output Pauli `w` is in `basis` at index `row`.
/// Outputs not in `basis` are dropped. Computed once per expm call and
/// reused, CSC-style, by [`CachedCscOp`] across every Krylov/Taylor matvec.
///
/// Uses `LindbladSpec::compute_action_terms` with the same reusable
/// scratch-buffer pattern as [`crate::orbit_rep::build_orbit_rep_cols`];
/// this is the exact same per-column action the retired matrix-free op fed
/// to `LindbladSpec::spmv_matrix_free`, so the numerics are unchanged.
fn build_mf_cols(
    spec: &LindbladSpec,
    basis: &[Word],
    index: &FxHashMap<Word, u32>,
) -> Vec<Vec<(u32, f64)>> {
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
                for (w, c) in terms.iter() {
                    if let Some(&row) = index.get(w) {
                        out.push((row, *c));
                    }
                }
                out
            },
        )
        .collect()
}

/// Borrowed CSC-style view of the in-basis-restricted real generator `M`,
/// backed by the cached per-column action [`build_mf_cols`] (computed once
/// per expm call).
///
/// `cols[c]` holds `(row, coeff)` pairs for column `c`; `dot` does a CSC
/// matvec `y = M·x` against the cached action, with no per-matvec action
/// recompute and no CSR materialisation. Mirrors
/// [`crate::orbit_rep::OrbitRepCscOp`] but for the real `f64` path.
///
/// Borrowed, not owned: `quspin-types` provides a blanket
/// `LinearOperator` impl for `&T`, so `ExpmOp::from_parts(op, …)` accepts a
/// `CachedCscOp` by value while it keeps borrowing `cols`.
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

        // Parallelise over column chunks: each thread accumulates into a
        // dense local `y` of length `dim`, reading the cached action.
        let partial_ys: Vec<Vec<f64>> = self
            .cols
            .par_chunks(chunk_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| {
                let c_offset = chunk_idx * chunk_size;
                let mut y_local = vec![0.0; n];
                for (c_local, col) in chunk.iter().enumerate() {
                    let c = c_offset + c_local;
                    let xc = input[c];
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
/// The generator action is computed ONCE here (via [`build_mf_cols`]) and
/// reused, CSC-style, across every Krylov/Taylor matvec (see
/// [`CachedCscOp`]) — instead of recomputing `spec.action` on every matvec.
/// One pass over the cached columns extracts the diagonal shift
/// `μ = tr(M)/n` and the exact column 1-norm of `M − μ·I`; from
/// `‖dt·(M−μI)‖₁` we pick the Taylor partition `(m*, s)` and hand everything
/// to [`quspin_expm::ExpmOp::from_parts`].
pub(crate) fn expm_apply_mf(spec: &LindbladSpec, basis: &[Word], dt: f64, coeffs: &[f64]) -> Vec<f64> {
    let n = basis.len();
    if n == 0 {
        return Vec::new();
    }

    let index = build_basis_index(basis);
    // Generator action: computed ONCE, reused across every matvec below.
    let cols = build_mf_cols(spec, basis, &index);

    // One pass over the cached columns: per column `c` gather
    // `raw = Σ|coeff|` and the diagonal entry `diag` (coeff of the term
    // whose row == the input column). From these we get `trace = Σ diag`,
    // `μ = trace/n`, and the column 1-norm of `M − μ·I`: per column it is
    // `raw − |diag| + |diag − μ|`.
    let per_col: Vec<(f64, f64)> = cols
        .par_iter()
        .enumerate()
        .map(|(c, col)| {
            let mut raw = 0.0;
            let mut diag = 0.0;
            for &(row, coeff) in col.iter() {
                raw += coeff.abs();
                if row as usize == c {
                    diag = coeff;
                }
            }
            (raw, diag)
        })
        .collect();

    let trace: f64 = per_col.iter().map(|(_, d)| *d).sum();
    let mu = trace / n as f64;
    let onenorm = per_col
        .iter()
        .map(|(raw, diag)| raw - diag.abs() + (diag - mu).abs())
        .fold(0.0_f64, f64::max);

    let (m_star, s) = expm::select_ms(dt.abs() * onenorm, None);

    let op = CachedCscOp {
        cols: &cols,
        dim: n,
    };
    let expm = quspin_expm::ExpmOp::from_parts(op, dt, mu, s as usize, m_star as usize, 1e-12_f64);

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
) -> Vec<Complex<f64>> {
    let n = basis.len();
    if n == 0 {
        return Vec::new();
    }
    let re: Vec<f64> = b.iter().map(|z| z.re).collect();
    let im: Vec<f64> = b.iter().map(|z| z.im).collect();
    let re_out = expm_apply_mf(spec, basis, dt, &re);
    let im_out = expm_apply_mf(spec, basis, dt, &im);
    re_out
        .into_iter()
        .zip(im_out)
        .map(|(r, i)| Complex::new(r, i))
        .collect()
}


