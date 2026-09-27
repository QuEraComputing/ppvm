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

/// Per-column `(raw, diag)` for the `μ`/1-norm selection: `raw` bounds
/// `Σ_r |M[r,c]|` from above and `diag = M[c,c]`.
type PerCol<T> = Vec<(f64, T)>;

/// Scratch buffers for [`LindbladSpec::compute_action_terms`].
type ActionScratch = (Vec<u32>, Vec<u32>, FxHashMap<Word, Complex<f64>>);

/// Consecutive CSC columns stored flat: local column `j` holds
/// `rows[offsets[j]..offsets[j + 1]]` and the matching `vals`.
struct CscBlock<T> {
    offsets: Vec<u32>,
    rows: Vec<u32>,
    vals: Vec<T>,
}

/// Cached in-basis action in CSC form, stored as blocks of `block` columns.
///
/// One exactly-sized allocation triple per block replaces one `Vec` per
/// column: at `|basis| ~ 10^6` the per-column layout reserved every `L*`
/// output (in- and out-of-basis) and left ~10^6 small allocations for the
/// system allocator to retain after the expm call.
pub(crate) struct BlockCsc<T> {
    blocks: Vec<CscBlock<T>>,
    block: usize,
    dim: usize,
}

impl<T> BlockCsc<T> {
    /// Visit the columns `range` in order as `(col, rows, vals)`.
    fn for_each_col(&self, range: std::ops::Range<usize>, mut f: impl FnMut(usize, &[u32], &[T])) {
        let mut c = range.start;
        while c < range.end {
            let b = &self.blocks[c / self.block];
            let base = (c / self.block) * self.block;
            let stop = range.end.min(base + b.offsets.len() - 1);
            for j in (c - base)..(stop - base) {
                let (lo, hi) = (b.offsets[j] as usize, b.offsets[j + 1] as usize);
                f(base + j, &b.rows[lo..hi], &b.vals[lo..hi]);
            }
            c = stop;
        }
    }
}

/// Build the [`BlockCsc`] cache and the per-column `(raw, diag)` data for
/// a `dim`-column generator in one parallel pass. `col(c, scratch, rows,
/// vals)` appends the in-basis entries of column `c` to `rows`/`vals` and
/// returns its `(raw, diag)`.
fn build_block_csc<T, F>(spec: &LindbladSpec, dim: usize, col: F) -> (BlockCsc<T>, PerCol<T>)
where
    T: Copy + Send + Sync,
    F: Fn(usize, &mut ActionScratch, &mut Vec<u32>, &mut Vec<T>) -> (f64, T) + Sync,
{
    // ~16 blocks per thread for load balance, but never so small that the
    // per-block allocations matter.
    let block = dim
        .div_ceil(16 * rayon::current_num_threads().max(1))
        .clamp(64, 4096);
    let (blocks, per_col): (Vec<CscBlock<T>>, Vec<PerCol<T>>) = (0..dim.div_ceil(block))
        .into_par_iter()
        .map_init(
            || {
                let scratch: ActionScratch = (
                    Vec::with_capacity(spec.n_qubits()),
                    Vec::with_capacity(128),
                    FxHashMap::with_capacity_and_hasher(128, FxBuildHasher::default()),
                );
                (scratch, Vec::<u32>::new(), Vec::<T>::new())
            },
            |(scratch, rows, vals), b| {
                let cols = (b * block)..dim.min((b + 1) * block);
                rows.clear();
                vals.clear();
                let mut offsets = Vec::with_capacity(cols.len() + 1);
                let mut per_col = Vec::with_capacity(cols.len());
                offsets.push(0);
                for c in cols {
                    per_col.push(col(c, scratch, rows, vals));
                    offsets.push(u32::try_from(rows.len()).expect("CSC block exceeds u32 entries"));
                }
                // `to_vec` sizes the stored block exactly; the staging
                // buffers are reused for the next block on this thread.
                let blk = CscBlock {
                    offsets,
                    rows: rows.to_vec(),
                    vals: vals.to_vec(),
                };
                (blk, per_col)
            },
        )
        .unzip();
    let per_col = per_col.into_iter().flatten().collect();
    (BlockCsc { blocks, block, dim }, per_col)
}

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
) -> (BlockCsc<f64>, PerCol<f64>) {
    build_block_csc(spec, basis.len(), |c, (s1, s2, lm), rows, vals| {
        let p = &basis[c];
        let terms = spec.compute_action_terms(p, s1, s2, lm);
        let mut raw = 0.0;
        let mut diag = 0.0;
        for (w, v) in terms.iter() {
            raw += v.abs();
            if w == p {
                diag = *v;
            }
            if let Some(&row) = index.get(w) {
                rows.push(row);
                vals.push(*v);
            }
        }
        (raw, diag)
    })
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
) -> (BlockCsc<Complex<f64>>, PerCol<Complex<f64>>) {
    build_block_csc(spec, basis.len(), |c, (s1, s2, lm), rows, vals| {
        let r = &basis[c];
        // A rep that cannot carry the sector has coefficient zero
        // identically, so its column is empty.
        let Some(orbit_in) = sector.orbit_size(r) else {
            return (0.0, Complex::new(0.0, 0.0));
        };
        let terms = spec.compute_action_terms(r, s1, s2, lm);
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
                rows.push(row);
                vals.push(val);
            }
        }
        (raw, diag)
    })
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
    pub(crate) cols: &'a BlockCsc<T>,
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
        self.cols.dim
    }

    fn parallel_hint(&self) -> bool {
        // `dot` parallelises internally over column chunks, and we drive the
        // sequential single-vector `apply` path; never let quspin run its
        // persistent-thread pool on top of our rayon parallelism.
        false
    }

    fn dot(&self, overwrite: bool, input: &[T], output: &mut [T]) -> Result<(), QuSpinError> {
        let n = self.cols.dim;
        if n == 0 {
            return Ok(());
        }
        let num_threads = rayon::current_num_threads().max(1);
        let chunk_size = n.div_ceil(num_threads);

        // Parallelise over column chunks; each thread accumulates into a dense
        // local `y` of length `dim`, reading the cached action; the partials
        // are reduced into `output` sequentially at the end.
        let partial_ys: Vec<Vec<T>> = (0..n.div_ceil(chunk_size))
            .into_par_iter()
            .map(|chunk_idx| {
                let cols = (chunk_idx * chunk_size)..n.min((chunk_idx + 1) * chunk_size);
                let mut y_local = vec![T::zero(); n];
                self.cols.for_each_col(cols, |c, rows, vals| {
                    let xc = input[c];
                    if xc == T::zero() {
                        return;
                    }
                    for (&row, &val) in rows.iter().zip(vals) {
                        y_local[row as usize] += val * xc;
                    }
                });
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
    cols: &BlockCsc<T>,
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
    let n = cols.dim;
    let trace: T = per_col.iter().map(|(_, d)| *d).sum();
    let mu = trace / n as f64;
    let onenorm = per_col
        .iter()
        .map(|&(raw, diag)| raw - diag.mag() + (diag - mu).mag())
        .fold(0.0_f64, f64::max);
    let (m_star, s, expm_tol) = select(dt.abs() * onenorm);

    let mut v = coeffs.to_vec();
    let op = CscOp { cols };
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
