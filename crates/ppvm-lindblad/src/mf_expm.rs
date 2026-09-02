// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Matrix-free `exp(dt · L*) · b`, driven by the external `quspin-expm`
//! crate — for both the real (`f64`) adaptive path and the complex,
//! phase-aware orbit-rep path.
//!
//! Instead of materialising the in-basis-restricted generator as a CSR, the
//! per-column generator action is computed ONCE per expm call (via
//! [`build_mf_cols`] / [`build_orbit_rep_cols`]) and reused, CSC-style,
//! across every Krylov/Taylor matvec
//! by [`CscOp`] (a [`quspin_types::LinearOperator`]) fed to
//! [`quspin_expm::ExpmOp::from_parts`]. Each matvec is then a cheap CSC
//! scatter; the Pauli-commutator action is never recomputed per matvec.
//! `from_parts` (rather than `ExpmOp::new`) supplies the diagonal shift `μ`,
//! the partition count `s`, and the truncation order `m*` directly, bypassing
//! quspin's adaptive parameter selection — so the 1-norm *estimator* and
//! `dot_transpose` are never invoked on the single-vector `apply` path; only
//! [`LinearOperator::dot`] runs.
//!
//! `μ`, the trace, and the column 1-norm of `A − μ·I` are computed in the
//! same single action pass as the cache, and turned into an `apply` by the
//! shared [`expm_apply_cached`] tail. The `(m, s)` Taylor partition is
//! picked with the tolerance-matched tables in [`crate::expm`]: a relaxed
//! `tol=1e-6` table when the PC prunes coarsely (`drop_tol ≥ 1e-4`), else the
//! double-precision table (keeping the exact-reference test paths bit-exact).

use crate::scalar::Coeff;
use crate::sector::Sector;
use crate::{LindbladSpec, Word, build_basis_index, expm};
use fxhash::{FxBuildHasher, FxHashMap};
use num::Complex;
use quspin_types::{ExpmComputation, LinearOperator, QuSpinError};
use rayon::prelude::*;
use std::iter::Sum;
use std::ops::{AddAssign, Div, Mul, Sub};

/// CSC columns of a cached in-basis action: `cols[c]` = `(row, coeff)`.
type Cols<T> = Vec<Vec<(u32, T)>>;
/// Per-column `(raw, diag)` for the `μ`/1-norm selection: `raw` bounds
/// `Σ_r |M[r,c]|` from above and `diag = M[c,c]`.
type PerCol<T> = Vec<(f64, T)>;

/// Per-column in-basis action of the real generator `M`, plus the data the
/// `(m, s)`/`μ` selection needs — all from ONE action pass over the basis.
///
/// Returns `(cols, per_col)` where `cols[c]` holds `(row, coeff)` for every
/// action output of `L*(basis[c])` that lands back in `basis` (CSC column
/// `c`), and `per_col[c] = (raw, diag)` with `raw = Σ|coeff|` over ALL action
/// outputs (in- and out-of-basis, an upper bound on the column 1-norm) and
/// `diag` the coefficient of the output Word equal to the input Word. The
/// cache is reused by [`CscOp`] across every Krylov/Taylor matvec.
fn build_mf_cols(
    spec: &LindbladSpec,
    basis: &[Word],
    index: &FxHashMap<Word, u32>,
) -> (Cols<f64>, PerCol<f64>) {
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

/// Per-column **phase-aware** action of the in-basis-restricted orbit-rep
/// generator `M` at momentum `sector`, plus the `(m, s)`/`μ` selection data
/// — from ONE action pass over the basis.
///
/// `cols[c]` holds `(row, χ_k(g_{cnt_q}) · v_q · |orbit_c| / |orbit_row|)`
/// for every action output Pauli `q` of `L*(basis[c])` whose orbit rep
/// `r_q` is in `basis` at index `row`; outputs whose rep is out of basis
/// are dropped. This is the expensive part of the orbit-rep dynamics
/// (`compute_action_terms`, [`Sector::canonicalize_phase`]).
///
/// The character-weighted sum runs over the *output* orbit's distinct
/// members, which makes it the generator in the **summing** convention
/// `ĉ_r = |orbit_r| · c_r`. Coefficients here are in the *averaged*
/// convention (`c_r` = the plain coefficient of the rep word, what
/// `canonicalize_pauli_sum_complex` produces), so each entry carries the
/// similarity factor `|orbit_c| / |orbit_row|` that converts between
/// them. It is 1 exactly when both orbits are free — hence the factor is
/// invisible until an orbit has a non-trivial stabilizer, and cannot be
/// hoisted out as a global `|G|`.
///
/// Unlike [`build_mf_cols`], `per_col[c].0` sums only the retained
/// in-basis entries — the exact column 1-norm of the restricted `M`, not an
/// upper bound: several distinct outputs `q` can share one rep, so the
/// out-of-basis magnitudes are not attributable to a column of `M`. `diag`
/// accumulates for the same reason.
fn build_orbit_rep_cols(
    spec: &LindbladSpec,
    basis: &[Word],
    index: &FxHashMap<Word, u32>,
    sector: Sector<'_>,
) -> (Cols<Complex<f64>>, PerCol<Complex<f64>>) {
    basis
        .par_iter()
        .enumerate()
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
            |(s1, s2, lm), (c, r)| {
                // A rep that cannot carry the sector has coefficient zero
                // identically, so its column is empty.
                let Some(orbit_in) = sector.orbit_size(r) else {
                    return (Vec::new(), (0.0, Complex::new(0.0, 0.0)));
                };
                let terms = spec.compute_action_terms(r, s1, s2, lm);
                let mut out = Vec::with_capacity(terms.len());
                let mut raw = 0.0;
                let mut diag = Complex::new(0.0, 0.0);
                for (q, v) in terms.iter() {
                    let Some((r_q, phase, orbit_out)) = sector.canonicalize_phase(q) else {
                        continue;
                    };
                    if let Some(&row) = index.get(&r_q) {
                        let val = phase * *v * (orbit_in as f64 / orbit_out as f64);
                        raw += val.norm();
                        if row as usize == c {
                            diag += val;
                        }
                        out.push((row, val));
                    }
                }
                (out, (raw, diag))
            },
        )
        .unzip()
}

/// Borrowed CSC-style view of an in-basis-restricted generator `M`, backed
/// by a cached per-column action computed once per expm call
/// ([`build_mf_cols`]). `dot` performs the CSC matvec `y = M·x` against the cache; the
/// remaining `LinearOperator` entry points are unused on the `from_parts` +
/// single-vector `apply` path.
///
/// Borrowed, not owned: `quspin-types` provides a blanket `LinearOperator`
/// impl for `&T`, so `ExpmOp::from_parts(op, ...)` accepts a `CscOp` by
/// value while it keeps borrowing `cols`.
pub(crate) struct CscOp<'a, T> {
    pub(crate) cols: &'a [Vec<(u32, T)>],
    pub(crate) dim: usize,
}

impl<T> LinearOperator<T> for CscOp<'_, T>
where
    T: ExpmComputation
        + Copy
        + PartialEq
        + num::Zero
        + std::ops::AddAssign
        + std::ops::Mul<Output = T>
        + Send
        + Sync,
{
    fn dim(&self) -> usize {
        self.dim
    }

    fn parallel_hint(&self) -> bool {
        // `dot` parallelises internally over column chunks, and we drive the
        // sequential single-vector `apply` path; never let quspin run its
        // persistent-thread pool on top of our rayon parallelism.
        false
    }

    fn dot(&self, overwrite: bool, input: &[T], output: &mut [T]) -> Result<(), QuSpinError> {
        let n = self.dim;
        if n == 0 {
            return Ok(());
        }
        let num_threads = rayon::current_num_threads().max(1);
        let chunk_size = n.div_ceil(num_threads);

        // Parallelise over column chunks; each thread accumulates into a dense
        // local `y` of length `dim`, reading the cached action; the partials
        // are reduced into `output` sequentially at the end.
        let partial_ys: Vec<Vec<T>> = self
            .cols
            .par_chunks(chunk_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| {
                let c_offset = chunk_idx * chunk_size;
                let mut y_local = vec![T::zero(); n];
                for (c_local, col) in chunk.iter().enumerate() {
                    let xc = input[c_offset + c_local];
                    if xc == T::zero() {
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
            output.fill(T::zero());
        }
        for partial in &partial_ys {
            for (oi, &pi) in output.iter_mut().zip(partial.iter()) {
                *oi += pi;
            }
        }
        Ok(())
    }

    fn trace(&self) -> T {
        // Computed eagerly by the callers; never reached on the
        // `from_parts` + single-vector `apply` path.
        unreachable!("CscOp::trace not used on the from_parts apply path")
    }

    fn onenorm(&self, _shift: T) -> <T as ExpmComputation>::Real {
        unreachable!("CscOp::onenorm not used on the from_parts apply path")
    }

    fn dot_transpose(
        &self,
        _overwrite: bool,
        _input: &[T],
        _output: &mut [T],
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CscOp: dot_transpose not used on the from_parts apply path".into(),
        ))
    }

    fn dot_many(
        &self,
        _overwrite: bool,
        _input: ndarray::ArrayView2<'_, T>,
        _output: ndarray::ArrayViewMut2<'_, T>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CscOp: dot_many not used on the from_parts apply path".into(),
        ))
    }

    fn dot_chunk(
        &self,
        _overwrite: bool,
        _input: &[T],
        _output_chunk: &mut [T],
        _row_start: usize,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CscOp: dot_chunk not used on the from_parts apply path".into(),
        ))
    }

    fn dot_transpose_chunk(
        &self,
        _input: &[T],
        _output: &[<T as ExpmComputation>::Atomic],
        _rows: std::ops::Range<usize>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CscOp: dot_transpose_chunk not used on the from_parts apply path".into(),
        ))
    }
}

/// Shared tail of every matrix-free expm: from the cached per-column action
/// derive the diagonal shift `μ = tr(M)/n` and a bound on the column 1-norm
/// of `M − μ·I` (`raw − |diag| + |diag − μ|` per column), pick the Taylor
/// partition via `select` from `‖dt·(M−μI)‖₁`, and hand everything to
/// [`quspin_expm::ExpmOp::from_parts`]. Returns `exp(dt · M) · coeffs`.
///
/// `select` maps `‖dt·(M−μI)‖₁` to `(m*, s, backward-error tol)`; the two
/// call sites differ only in that choice.
fn expm_apply_cached<T>(
    cols: &Cols<T>,
    per_col: &PerCol<T>,
    dt: f64,
    coeffs: &[T],
    select: impl FnOnce(f64) -> (u32, u32, f64),
) -> Vec<T>
where
    T: ExpmComputation<Real = f64>
        + Coeff
        + PartialEq
        + num::Zero
        + AddAssign
        + Mul<Output = T>
        + Sub<Output = T>
        + Div<f64, Output = T>
        + From<f64>
        + Sum,
{
    let n = cols.len();
    let trace: T = per_col.iter().map(|(_, d)| *d).sum();
    let mu = trace / n as f64;
    let onenorm = per_col
        .iter()
        .map(|&(raw, diag)| raw - diag.mag() + (diag - mu).mag())
        .fold(0.0_f64, f64::max);
    let (m_star, s, expm_tol) = select(dt.abs() * onenorm);

    let mut v = coeffs.to_vec();
    let op = CscOp { cols, dim: n };
    let expm =
        quspin_expm::ExpmOp::from_parts(op, T::from(dt), mu, s as usize, m_star as usize, expm_tol);
    expm.apply(ndarray::ArrayViewMut1::from(v.as_mut_slice()))
        .expect("expm apply");
    v
}

/// Compute `exp(dt · M) · coeffs` for the in-basis-restricted generator
/// `M`, matrix-free, via `quspin-expm`. Returns a fresh `Vec<f64>` of length
/// `basis.len()`.
///
/// ONE action pass builds the CSC cache `cols` (reused across every matvec)
/// and, in the same pass, the `(raw, diag)` data the `μ`/1-norm selection
/// needs; [`expm_apply_cached`] does the rest.
pub(crate) fn expm_apply_mf(
    spec: &LindbladSpec,
    basis: &[Word],
    dt: f64,
    coeffs: &[f64],
    drop_tol: f64,
) -> Vec<f64> {
    if basis.is_empty() {
        return Vec::new();
    }
    let index = build_basis_index(basis);
    let (cols, per_col) = build_mf_cols(spec, basis, &index);

    // Pick the Taylor backward-error tolerance to match the basis truncation:
    // when the PC prunes coarsely (drop_tol >= 1e-4) a double-precision exp is
    // ~10 orders more accurate than the state it acts on, so the relaxed
    // (tol=1e-6, still >=100x tighter than the cut) table is used — it admits a
    // lower-degree Taylor polynomial and cuts the SpMV count with no effect on
    // the truncated result. At tight/zero drop_tol we keep double precision so
    // the exact-reference paths (orbit-rep / merged) still agree bit-for-bit.
    expm_apply_cached(&cols, &per_col, dt, coeffs, |t_norm| {
        if drop_tol >= 1e-4 {
            let (m, s) = expm::select_ms_loose(t_norm);
            (m, s, 1e-6)
        } else {
            let (m, s) = expm::select_ms(t_norm);
            (m, s, 1e-12)
        }
    })
}

/// Compute `exp(dt · M) · coeffs` for the in-basis-restricted **orbit-rep**
/// generator `M` at momentum `sector`, via `quspin-expm`. Returns a fresh
/// `Vec<Complex<f64>>` of length `basis.len()`.
///
/// The expensive phase-aware action is computed ONCE here (via
/// [`build_orbit_rep_cols`]) and reused, CSC-style, across every
/// Krylov–Taylor matvec, exactly as on the real path.
pub(crate) fn expm_apply_orbit_rep(
    spec: &LindbladSpec,
    basis: &[Word],
    sector: Sector<'_>,
    dt: f64,
    coeffs: &[Complex<f64>],
) -> Vec<Complex<f64>> {
    if basis.is_empty() {
        return Vec::new();
    }
    let index = build_basis_index(basis);
    let (cols, per_col) = build_orbit_rep_cols(spec, basis, &index, sector);

    expm_apply_cached(&cols, &per_col, dt, coeffs, |t_norm| {
        let (m, s) = expm::select_ms(t_norm);
        (m, s, 1e-12)
    })
}

/// `exp(dt · M) · b` where `M` is the REAL in-basis-restricted generator but
/// the input vector `b` is complex. Because `M` is real,
/// `exp(dt·M)·(re + i·im) = exp(dt·M)·re + i·exp(dt·M)·im`, so we split the
/// complex vector into its real and imaginary parts, run two real
/// matrix-free applies, and recombine. Used by the test-only full-space
/// complex reference step.
#[cfg(test)]
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
