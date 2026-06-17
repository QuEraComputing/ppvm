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

use crate::expm::CsrCx;
use crate::{LindbladSpec, Word, build_basis_index, expm};
use fxhash::FxHashMap;
use num::Complex;
use quspin_types::{ExpmComputation, LinearOperator, QuSpinError};
use rayon::prelude::*;

/// Borrowed matrix-free view of the in-basis-restricted Lindbladian
/// generator `M` (the same matrix [`LindbladSpec::generator_csr`] would
/// build, never materialised).
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

/// Borrowed CSR-backed view of a genuinely complex generator (the orbit-rep
/// generator built by [`LindbladSpec::generator_csr_orbit_rep`], whose
/// momentum-character phases make the entries complex). There is no
/// matrix-free complex action available, so this drives the materialised
/// `CsrCx` through `quspin-expm`.
pub(crate) struct CsrCxOp<'a> {
    csr: &'a CsrCx,
}

impl LinearOperator<Complex<f64>> for CsrCxOp<'_> {
    fn dim(&self) -> usize {
        self.csr.rows()
    }

    fn parallel_hint(&self) -> bool {
        // `spmv_cx` parallelises internally over rows when worthwhile; we
        // drive the sequential single-vector `apply` path.
        false
    }

    fn dot(
        &self,
        overwrite: bool,
        input: &[Complex<f64>],
        output: &mut [Complex<f64>],
    ) -> Result<(), QuSpinError> {
        if overwrite {
            expm::spmv_cx(self.csr, input, output, usize::MAX);
        } else {
            let mut tmp = vec![Complex::new(0.0, 0.0); output.len()];
            expm::spmv_cx(self.csr, input, &mut tmp, usize::MAX);
            for (o, t) in output.iter_mut().zip(tmp.iter()) {
                *o += *t;
            }
        }
        Ok(())
    }

    fn trace(&self) -> Complex<f64> {
        unreachable!("CsrCxOp::trace not used on the from_parts apply path")
    }

    fn onenorm(&self, _shift: Complex<f64>) -> f64 {
        unreachable!("CsrCxOp::onenorm not used on the from_parts apply path")
    }

    fn dot_transpose(
        &self,
        _overwrite: bool,
        _input: &[Complex<f64>],
        _output: &mut [Complex<f64>],
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CsrCxOp: dot_transpose not used on the from_parts apply path".into(),
        ))
    }

    fn dot_many(
        &self,
        _overwrite: bool,
        _input: ndarray::ArrayView2<'_, Complex<f64>>,
        _output: ndarray::ArrayViewMut2<'_, Complex<f64>>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CsrCxOp: dot_many not used on the from_parts apply path".into(),
        ))
    }

    fn dot_chunk(
        &self,
        _overwrite: bool,
        _input: &[Complex<f64>],
        _output_chunk: &mut [Complex<f64>],
        _row_start: usize,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CsrCxOp: dot_chunk not used on the from_parts apply path".into(),
        ))
    }

    fn dot_transpose_chunk(
        &self,
        _input: &[Complex<f64>],
        _output: &[<Complex<f64> as ExpmComputation>::Atomic],
        _rows: std::ops::Range<usize>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CsrCxOp: dot_transpose_chunk not used on the from_parts apply path".into(),
        ))
    }
}

/// Compute `exp(dt · M) · coeffs` for a genuinely complex CSR generator `M`
/// (the orbit-rep generator) via `quspin-expm`. Returns a fresh
/// `Vec<Complex<f64>>` of length `csr.rows()`.
///
/// Uses the exact zero shift `μ = 0` (the trace shift is only an efficiency
/// optimisation and orbit-rep matrices are small). The Taylor partition
/// `(m*, s)` is picked from `‖dt·M‖₁` via the same selection table as the
/// retired hand-rolled engine in [`crate::expm`].
pub(crate) fn expm_apply_csr_cx(
    csr: &CsrCx,
    dt: f64,
    coeffs: &[Complex<f64>],
) -> Vec<Complex<f64>> {
    let n = csr.rows();
    if n == 0 {
        return Vec::new();
    }

    let mu = Complex::new(0.0, 0.0);
    let onenorm = expm::csr_cx_one_norm(csr);
    let (m_star, s) = expm::select_ms(dt.abs() * onenorm, None);

    let op = CsrCxOp { csr };
    let e = quspin_expm::ExpmOp::from_parts(
        op,
        Complex::new(dt, 0.0),
        mu,
        s as usize,
        m_star as usize,
        1e-12_f64,
    );

    let mut v = coeffs.to_vec();
    e.apply(ndarray::ArrayViewMut1::from(v.as_mut_slice()))
        .expect("expm apply (csr cx)");
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expm::{self, CsrCx, ExpmOpts};
    use num::Complex;
    use sprs::{CsMatI, TriMatI};

    /// Build a small random real CSR + a random `LindbladSpec`-free check is
    /// hard (action depends on a spec), so the real parity test lives in
    /// `lib.rs` where a spec is available. Here we cover the complex
    /// orbit-rep path against the retired engine.
    fn rand_csr_cx(n: usize, density: f64, seed: &mut u64) -> CsrCx {
        let mut next = || {
            // xorshift64
            let mut x = *seed;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            *seed = x;
            (x >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut tri = TriMatI::<Complex<f64>, u32>::new((n, n));
        for i in 0..n {
            for j in 0..n {
                if next() < density {
                    let re = next() * 2.0 - 1.0;
                    let im = next() * 2.0 - 1.0;
                    tri.add_triplet(i, j, Complex::new(re, im));
                }
            }
        }
        tri.to_csr()
    }

    #[test]
    fn parity_complex_csr_against_retired_engine() {
        let mut worst: f64 = 0.0;
        for &n in &[1usize, 3, 7, 12] {
            let mut s2 = 0x1234_5678_9abc_def0u64.wrapping_add(n as u64 * 0x9e37_79b9);
            let csr: CsMatI<Complex<f64>, u32, usize> = rand_csr_cx(n, 0.5, &mut s2);
            let dt = 0.37;
            let mut next = || {
                s2 ^= s2 << 13;
                s2 ^= s2 >> 7;
                s2 ^= s2 << 17;
                (s2 >> 11) as f64 / (1u64 << 53) as f64
            };
            let coeffs: Vec<Complex<f64>> = (0..n)
                .map(|_| Complex::new(next() * 2.0 - 1.0, next() * 2.0 - 1.0))
                .collect();

            let new = expm_apply_csr_cx(&csr, dt, &coeffs);
            let old = expm::expm_multiply_cx(&csr, dt, &coeffs, ExpmOpts::default());

            assert_eq!(new.len(), old.len());
            for (a, b) in new.iter().zip(old.iter()) {
                let d = (a - b).norm();
                worst = worst.max(d);
                assert!(
                    d < 1e-9,
                    "complex CSR parity mismatch (n={n}): new={a} old={b} |Δ|={d:e}"
                );
            }
        }
        eprintln!("parity_complex_csr_against_retired_engine: worst |Δ| = {worst:e}");
    }
}
