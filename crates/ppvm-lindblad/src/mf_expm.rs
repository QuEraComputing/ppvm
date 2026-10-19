// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Matrix-free `exp(dt · L*) · b`, driven by the external `quspin-expm`
//! crate — for both the real (`f64`) adaptive path and the complex,
//! phase-aware orbit-rep path.
//!
//! The in-basis-restricted generator's per-column action is computed ONCE
//! per expm call (via
//! [`build_mf_cols`] / [`build_orbit_rep_cols`]), stored row-major, and
//! reused across every Krylov/Taylor matvec
//! by [`CsrOp`] (a [`quspin_types::LinearOperator`]) fed to
//! [`quspin_expm::ExpmOp::from_parts`]. Each matvec is then a cheap,
//! deterministic CSR product; the Pauli-commutator action is never
//! recomputed per matvec.
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
type ActionScratch<const C: usize> = (Vec<u32>, Vec<u32>, FxHashMap<Word<C>, Complex<f64>>);

/// One `(row, col, val)` entry of the generator, staged between the
/// column-wise action pass and the row-major layout.
type Entry<T> = (u32, u32, T);

/// Staged entries indexed by row range (one bucket per [`CsrRange`]).
type Buckets<T> = Vec<Vec<Entry<T>>>;

/// Consecutive CSR rows stored flat: local row `i` holds
/// `cols[offsets[i]..offsets[i + 1]]` and the matching `vals`, in
/// ascending column order.
struct CsrRange<T> {
    offsets: Vec<u32>,
    cols: Vec<u32>,
    vals: Vec<T>,
}

/// Cached in-basis action in row-major (CSR) form, stored as ranges of
/// `1 << range_shift` consecutive rows.
///
/// Row-major storage makes the matvec deterministic: each output entry is
/// summed by one thread over its row in ascending column order, so `M·x`
/// is bit-identical for every rayon thread count, with no per-thread
/// partial vectors and no serial reduction. `range_shift` depends on
/// `dim` only.
///
/// One exactly-sized allocation triple per range replaces one `Vec` per
/// row: at `|basis| ~ 10^6` a per-row layout would leave ~10^6 small
/// allocations for the system allocator to retain after the expm call.
pub(crate) struct BlockCsr<T> {
    ranges: Vec<CsrRange<T>>,
    range_shift: u32,
    dim: usize,
}

/// `log2` of the rows per [`CsrRange`]: the smallest power of two that
/// gives at most 256 ranges, and at least 64 rows per range. A function of
/// `dim` alone, so the layout (and every row's summation order) is
/// independent of the thread count. A power of two turns the per-entry
/// range lookup into a shift.
fn range_shift(dim: usize) -> u32 {
    dim.div_ceil(256)
        .max(64)
        .next_power_of_two()
        .trailing_zeros()
}

/// Build the [`BlockCsr`] cache and the per-column `(raw, diag)` data for
/// a `dim`-column generator. `col(c, scratch, rows, vals)` appends the
/// in-basis entries of column `c` to `rows`/`vals` and returns its
/// `(raw, diag)`.
///
/// The action is naturally column-wise, so the build runs in two parallel
/// passes. (1) Over column blocks: compute the action and stage each entry
/// in an exactly-sized bucket of its row range, in column order. (2) Over
/// row ranges: a stable counting sort of the range's buckets (taken in
/// block order) into CSR, dropping each bucket once consumed. The cache is
/// never held in both layouts at once, but staging stores the row as well:
/// 16 B/entry for real `T` and 24 B/entry for complex `T`, against
/// 12 B/entry and 20 B/entry for the final CSR.
fn build_block_csr<const C: usize, T, F>(
    spec: &LindbladSpec<C>,
    dim: usize,
    col: F,
) -> (BlockCsr<T>, PerCol<T>)
where
    T: Coeff,
    F: Fn(usize, &mut ActionScratch<C>, &mut Vec<u32>, &mut Vec<T>) -> (f64, T) + Sync,
{
    // ~16 blocks per thread for load balance, but never so small that the
    // per-block allocations matter. Block boundaries do not affect the
    // result: pass 2 restores ascending column order within every row.
    let block = dim
        .div_ceil(16 * rayon::current_num_threads().max(1))
        .clamp(64, 4096);
    let shift = range_shift(dim);
    let n_ranges = dim.div_ceil(1 << shift);

    // Pass 1: `staged[b][k]` holds block `b`'s entries in row range `k`.
    let (staged, per_col): (Vec<Buckets<T>>, Vec<PerCol<T>>) = (0..dim.div_ceil(block))
        .into_par_iter()
        .map_init(
            || {
                let scratch: ActionScratch<C> = (
                    Vec::with_capacity(spec.n_qubits()),
                    Vec::with_capacity(128),
                    FxHashMap::with_capacity_and_hasher(128, FxBuildHasher::default()),
                );
                (
                    scratch,
                    Vec::<u32>::new(),
                    Vec::<T>::new(),
                    Vec::<usize>::new(),
                    vec![0usize; n_ranges],
                )
            },
            |(scratch, rows, vals, ends, counts), b| {
                let cols = (b * block)..dim.min((b + 1) * block);
                rows.clear();
                vals.clear();
                ends.clear();
                let mut per_col = Vec::with_capacity(cols.len());
                for c in cols.clone() {
                    per_col.push(col(c, scratch, rows, vals));
                    ends.push(rows.len());
                }
                counts.fill(0);
                for &r in rows.iter() {
                    counts[(r >> shift) as usize] += 1;
                }
                let mut buckets: Buckets<T> =
                    counts.iter().map(|&k| Vec::with_capacity(k)).collect();
                let mut lo = 0;
                for (c, &hi) in cols.zip(ends.iter()) {
                    for e in lo..hi {
                        let r = rows[e];
                        buckets[(r >> shift) as usize].push((r, c as u32, vals[e]));
                    }
                    lo = hi;
                }
                (buckets, per_col)
            },
        )
        .unzip();
    let per_col = per_col.into_iter().flatten().collect();

    // Regroup the bucket handles by row range (moves `Vec` headers only).
    let mut by_range: Vec<Buckets<T>> = (0..n_ranges)
        .map(|_| Vec::with_capacity(staged.len()))
        .collect();
    for blk in staged {
        for (k, bucket) in blk.into_iter().enumerate() {
            by_range[k].push(bucket);
        }
    }

    // Pass 2: stable counting sort by row. Buckets are visited in block
    // order and each is in column order, so every row ends up in
    // ascending column order.
    let ranges = by_range
        .into_par_iter()
        .enumerate()
        .map(|(k, buckets)| {
            let row0 = k << shift;
            let nrows = (1usize << shift).min(dim - row0);
            let nnz: usize = buckets.iter().map(Vec::len).sum();
            // The per-range offsets are u32; check once so the counts
            // below cannot wrap.
            u32::try_from(nnz).expect("CSR row range exceeds u32 entries");
            let mut offsets = vec![0u32; nrows + 1];
            for &(r, _, _) in buckets.iter().flatten() {
                offsets[r as usize - row0 + 1] += 1;
            }
            for i in 0..nrows {
                offsets[i + 1] += offsets[i];
            }
            let mut cursor = offsets[..nrows].to_vec();
            let mut cols = vec![0u32; nnz];
            let mut vals = vec![T::zero(); nnz];
            for bucket in buckets {
                for (r, c, v) in bucket {
                    let i = r as usize - row0;
                    let p = cursor[i] as usize;
                    cols[p] = c;
                    vals[p] = v;
                    cursor[i] += 1;
                }
            }
            CsrRange {
                offsets,
                cols,
                vals,
            }
        })
        .collect();
    (
        BlockCsr {
            ranges,
            range_shift: shift,
            dim,
        },
        per_col,
    )
}

/// Per-column in-basis action of the real generator `M`, plus the data the
/// `(m, s)`/`μ` selection needs — all from ONE action pass over the basis.
///
/// Returns `(mat, per_col)`. Column `c` of the row-major cache `mat` holds
/// `M[row, c] = coeff` for every action output of `L*(basis[c])` that lands
/// back in `basis` at index `row`. `per_col[c] = (raw, diag)`, with
/// `raw = Σ|coeff|` over ALL action outputs (in- and out-of-basis, an upper
/// bound on the column 1-norm) and `diag` the coefficient of the output
/// Word equal to the input Word. Both sums run in the canonical term order
/// of [`LindbladSpec::compute_action_terms`]. The cache is reused by
/// [`CsrOp`] across every Krylov/Taylor matvec.
fn build_mf_cols<const C: usize>(
    spec: &LindbladSpec<C>,
    basis: &[Word<C>],
    index: &FxHashMap<Word<C>, u32>,
) -> (BlockCsr<f64>, PerCol<f64>) {
    build_block_csr(spec, basis.len(), |c, (s1, s2, lm), rows, vals| {
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
/// Column `c` of the row-major cache receives
/// `χ_k(g_{cnt_q}) · v_q · |orbit_c| / |orbit_row|` at row `row` for every
/// action output Pauli `q` of `L*(basis[c])` whose orbit rep `r_q` is in
/// `basis` at index `row`; outputs whose rep is out of basis are dropped.
/// Several outputs can share a rep, so a row can hold the same column more
/// than once; those duplicates are emitted, and later summed, in the
/// canonical term order of [`LindbladSpec::compute_action_terms`]. This is the expensive part of the orbit-rep dynamics
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
fn build_orbit_rep_cols<const C: usize>(
    spec: &LindbladSpec<C>,
    basis: &[Word<C>],
    index: &FxHashMap<Word<C>, u32>,
    sector: &Sector<'_>,
) -> (BlockCsr<Complex<f64>>, PerCol<Complex<f64>>) {
    build_block_csr(spec, basis.len(), |c, (s1, s2, lm), rows, vals| {
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

/// Borrowed CSR view of an in-basis-restricted generator `M`, backed by a
/// cached action computed once per expm call ([`build_mf_cols`]). `dot`
/// performs the row-parallel matvec `y = M·x` against the cache; the
/// remaining `LinearOperator` entry points are unused on the `from_parts` +
/// single-vector `apply` path.
///
/// Borrowed, not owned: `quspin-types` provides a blanket `LinearOperator`
/// impl for `&T`, so `ExpmOp::from_parts(op, ...)` accepts a `CsrOp` by
/// value while it keeps borrowing `mat`.
pub(crate) struct CsrOp<'a, T> {
    pub(crate) mat: &'a BlockCsr<T>,
}

impl<T> LinearOperator<T> for CsrOp<'_, T>
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
        self.mat.dim
    }

    fn parallel_hint(&self) -> bool {
        // `dot` parallelises internally over row ranges, and we drive the
        // sequential single-vector `apply` path; never let quspin run its
        // persistent-thread pool on top of our rayon parallelism.
        false
    }

    fn dot(&self, overwrite: bool, input: &[T], output: &mut [T]) -> Result<(), QuSpinError> {
        // Rows in parallel, one range per task; each output entry is summed
        // by one thread in ascending column order, so the result does not
        // depend on the thread count.
        let m = self.mat;
        output
            .par_chunks_mut(1 << m.range_shift)
            .zip(m.ranges.par_iter())
            .for_each(|(y, r)| {
                for (i, yi) in y.iter_mut().enumerate() {
                    let (lo, hi) = (r.offsets[i] as usize, r.offsets[i + 1] as usize);
                    let mut acc = T::zero();
                    for (&c, &v) in r.cols[lo..hi].iter().zip(&r.vals[lo..hi]) {
                        // Skip zero inputs, as the column-wise product did:
                        // saves the multiply-add for sparse `input`, and
                        // keeps a non-finite entry of an unused column out
                        // of `y`.
                        let xc = input[c as usize];
                        if xc != T::zero() {
                            acc += v * xc;
                        }
                    }
                    if overwrite {
                        *yi = acc;
                    } else {
                        *yi += acc;
                    }
                }
            });
        Ok(())
    }

    fn trace(&self) -> T {
        // Computed eagerly by the callers; never reached on the
        // `from_parts` + single-vector `apply` path.
        unreachable!("CsrOp::trace not used on the from_parts apply path")
    }

    fn onenorm(&self, _shift: T) -> <T as ExpmComputation>::Real {
        unreachable!("CsrOp::onenorm not used on the from_parts apply path")
    }

    fn dot_transpose(
        &self,
        _overwrite: bool,
        _input: &[T],
        _output: &mut [T],
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CsrOp: dot_transpose not used on the from_parts apply path".into(),
        ))
    }

    fn dot_many(
        &self,
        _overwrite: bool,
        _input: ndarray::ArrayView2<'_, T>,
        _output: ndarray::ArrayViewMut2<'_, T>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CsrOp: dot_many not used on the from_parts apply path".into(),
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
            "CsrOp: dot_chunk not used on the from_parts apply path".into(),
        ))
    }

    fn dot_transpose_chunk(
        &self,
        _input: &[T],
        _output: &[<T as ExpmComputation>::Atomic],
        _rows: std::ops::Range<usize>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "CsrOp: dot_transpose_chunk not used on the from_parts apply path".into(),
        ))
    }
}

/// Shared tail of every matrix-free expm: from the cached action `mat` and
/// its per-column `(raw, diag)` data derive the diagonal shift `μ = tr(M)/n` and a bound on the column 1-norm
/// of `M − μ·I` (`raw − |diag| + |diag − μ|` per column), pick the Taylor
/// partition via `select` from `‖dt·(M−μI)‖₁`, and hand everything to
/// [`quspin_expm::ExpmOp::from_parts`]. Returns `exp(dt · M) · coeffs`.
///
/// `select` maps `‖dt·(M−μI)‖₁` to `(m*, s, backward-error tol)`; the two
/// call sites differ only in that choice.
fn expm_apply_cached<T>(
    mat: &BlockCsr<T>,
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
    let n = mat.dim;
    let trace: T = per_col.iter().map(|(_, d)| *d).sum();
    let mu = trace / n as f64;
    let onenorm = per_col
        .iter()
        .map(|&(raw, diag)| raw - diag.mag() + (diag - mu).mag())
        .fold(0.0_f64, f64::max);
    let (m_star, s, expm_tol) = select(dt.abs() * onenorm);

    let mut v = coeffs.to_vec();
    let op = CsrOp { mat };
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
/// ONE action pass builds the CSR cache `mat` (reused across every matvec)
/// and, in the same pass, the `(raw, diag)` data the `μ`/1-norm selection
/// needs; [`expm_apply_cached`] does the rest.
pub(crate) fn expm_apply_mf<const C: usize>(
    spec: &LindbladSpec<C>,
    basis: &[Word<C>],
    dt: f64,
    coeffs: &[f64],
    drop_tol: f64,
) -> Vec<f64> {
    if basis.is_empty() {
        return Vec::new();
    }
    let index = build_basis_index(basis);
    let (mat, per_col) = build_mf_cols(spec, basis, &index);

    // Pick the Taylor backward-error tolerance to match the basis truncation:
    // when the PC prunes coarsely (drop_tol >= 1e-4) a double-precision exp is
    // ~10 orders more accurate than the state it acts on, so the relaxed
    // (tol=1e-6, still >=100x tighter than the cut) table is used — it admits a
    // lower-degree Taylor polynomial and cuts the SpMV count with no effect on
    // the truncated result. At tight/zero drop_tol we keep double precision so
    // the exact-reference paths (orbit-rep / merged) still agree bit-for-bit.
    expm_apply_cached(&mat, &per_col, dt, coeffs, |t_norm| {
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
/// [`build_orbit_rep_cols`]) and reused, row-major, across every
/// Krylov–Taylor matvec, exactly as on the real path.
pub(crate) fn expm_apply_orbit_rep<const C: usize>(
    spec: &LindbladSpec<C>,
    basis: &[Word<C>],
    sector: &Sector<'_>,
    dt: f64,
    coeffs: &[Complex<f64>],
) -> Vec<Complex<f64>> {
    if basis.is_empty() {
        return Vec::new();
    }
    let index = build_basis_index(basis);
    let (mat, per_col) = build_orbit_rep_cols(spec, basis, &index, sector);

    expm_apply_cached(&mat, &per_col, dt, coeffs, |t_norm| {
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
