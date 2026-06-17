// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Matrix-free `exp(dt · L*) · b` for the real (`f64`) path, driven by the
//! external `quspin-expm` crate.
//!
//! Instead of materialising the in-basis-restricted generator as a CSR, we
//! wrap the Lindbladian action in a [`MfOp`] implementing
//! [`quspin_types::LinearOperator`] and feed it to
//! [`quspin_expm::ExpmOp::from_parts`]. Using `from_parts` (rather than
//! `ExpmOp::new`) supplies the diagonal shift `μ`, the partition count `s`,
//! and the truncation order `m*` directly, which BYPASSES quspin's adaptive
//! parameter selection. As a consequence the 1-norm *estimator* and
//! `dot_transpose` are never invoked on the single-vector `apply` path — only
//! [`LinearOperator::dot`] runs, inside the Taylor/Horner loop.
//!
//! The shift `μ` and the exact column 1-norm of `A − μ·I` are computed in one
//! matrix-free pass over the basis (the same `(m, s)` selection table as the
//! retired hand-rolled engine in [`crate::expm`]).

use crate::{LindbladSpec, Word, build_basis_index, expm};
use fxhash::FxHashMap;
use quspin_types::{ExpmComputation, LinearOperator, QuSpinError};
use rayon::prelude::*;

/// Borrowed matrix-free view of the in-basis-restricted Lindbladian
/// generator `M` (the in-basis-restricted generator `M`, never
/// materialised).
///
/// Borrowed, not owned: `quspin-types` provides a blanket
/// `LinearOperator` impl for `&T`, so `ExpmOp::from_parts(op, …)` accepts a
/// `MfOp` by value while it keeps borrowing `spec` / `basis`.
pub(crate) struct MfOp<'a> {
    spec: &'a LindbladSpec,
    basis: &'a [Word],
    /// `Word → row` map for `basis`; built once and reused across every
    /// `dot` in the Taylor loop.
    index: FxHashMap<Word, u32>,
}

impl LinearOperator<f64> for MfOp<'_> {
    fn dim(&self) -> usize {
        self.basis.len()
    }

    fn parallel_hint(&self) -> bool {
        // `dot` parallelises internally over basis columns, and we drive the
        // sequential single-vector `apply` path; never let quspin run its
        // persistent-thread pool on top of our rayon parallelism.
        false
    }

    fn dot(&self, overwrite: bool, input: &[f64], output: &mut [f64]) -> Result<(), QuSpinError> {
        if overwrite {
            self.spec
                .spmv_matrix_free(self.basis, &self.index, input, output);
        } else {
            let mut tmp = vec![0.0; output.len()];
            self.spec
                .spmv_matrix_free(self.basis, &self.index, input, &mut tmp);
            for (o, t) in output.iter_mut().zip(tmp.iter()) {
                *o += *t;
            }
        }
        Ok(())
    }

    fn trace(&self) -> f64 {
        // Computed eagerly in `expm_apply_mf`; never reached on the
        // `from_parts` + single-vector `apply` path.
        unreachable!("MfOp::trace not used on the from_parts apply path")
    }

    fn onenorm(&self, _shift: f64) -> f64 {
        unreachable!("MfOp::onenorm not used on the from_parts apply path")
    }

    fn dot_transpose(
        &self,
        _overwrite: bool,
        _input: &[f64],
        _output: &mut [f64],
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "MfOp: dot_transpose not used on the from_parts apply path".into(),
        ))
    }

    fn dot_many(
        &self,
        _overwrite: bool,
        _input: ndarray::ArrayView2<'_, f64>,
        _output: ndarray::ArrayViewMut2<'_, f64>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "MfOp: dot_many not used on the from_parts apply path".into(),
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
            "MfOp: dot_chunk not used on the from_parts apply path".into(),
        ))
    }

    fn dot_transpose_chunk(
        &self,
        _input: &[f64],
        _output: &[<f64 as ExpmComputation>::Atomic],
        _rows: std::ops::Range<usize>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "MfOp: dot_transpose_chunk not used on the from_parts apply path".into(),
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
pub(crate) fn expm_apply_mf(spec: &LindbladSpec, basis: &[Word], dt: f64, coeffs: &[f64]) -> Vec<f64> {
    let n = basis.len();
    if n == 0 {
        return Vec::new();
    }

    // Single matrix-free pass: per column gather `raw = Σ|coeff|` and the
    // diagonal entry `diag` (coeff of the term whose output Word == the input
    // Word). From these we get `trace = Σ diag`, `μ = trace/n`, and the
    // column 1-norm of `M − μ·I`: per column it is
    // `raw − |diag| + |diag − μ|`.
    let per_col: Vec<(f64, f64)> = basis
        .par_iter()
        .map(|p| {
            let terms = spec.action(p);
            let mut raw = 0.0;
            let mut diag = 0.0;
            for (w, c) in &terms {
                raw += c.abs();
                if w == p {
                    diag = *c;
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

    let op = MfOp {
        spec,
        basis,
        index: build_basis_index(basis),
    };
    let expm = quspin_expm::ExpmOp::from_parts(op, dt, mu, s as usize, m_star as usize, 1e-12_f64);

    let mut v = coeffs.to_vec();
    expm.apply(ndarray::ArrayViewMut1::from(v.as_mut_slice()))
        .expect("expm apply");
    v
}
