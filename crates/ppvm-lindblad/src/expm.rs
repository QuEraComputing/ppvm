// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Action of the matrix exponential, `y = exp(t·A)·b`, for real sparse `A`.
//!
//! Implements Al-Mohy & Higham (2011), Algorithm 3.2: a Horner-form
//! degree-`m` Taylor polynomial evaluated `s` times (scaling-and-squaring),
//! with `(m, s)` chosen from a precomputed table that minimises total SpMV
//! count subject to a double-precision truncation bound.
//!
//! Storage is `Csr`. The hot path is repeated sparse-matrix ×
//! dense-vector products; SpMV is parallelised over rows with `rayon` when
//! the matrix has more than [`ExpmOpts::parallel_threshold`] nonzeros (below
//! that the task-spawn overhead beats the cache-friendly serial loop).
//! Both branches produce bit-identical output.

use num::Complex;
use rayon::prelude::*;
use sprs::CsMatI;

/// Sparse CSR matrix with `u32` column indices and `usize` indptr.
///
/// `u32` keeps memory traffic in SpMV minimal (4 bytes per nonzero for the
/// column index instead of 8) — important for bandwidth-bound regimes.
/// `usize` indptr accommodates `nnz > 2^32` even when individual indices fit.
pub type Csr = CsMatI<f64, u32, usize>;

/// `θ_m` table from Al-Mohy & Higham (2011), Table A.3, for double
/// precision (unit roundoff `u = 2^{-53}`).
///
/// `θ_m` bounds `‖A‖₁` such that the degree-`m` Taylor polynomial
/// approximates `exp(A)` to within `u`. We pick `(m, s)` with
/// `s ≥ ⌈‖tA‖₁ / θ_m⌉` and minimise `m·s` (total SpMV count).
const THETA: &[(u32, f64)] = &[
    (1, 2.29e-16),
    (2, 2.58e-8),
    (3, 1.39e-5),
    (4, 3.40e-4),
    (5, 2.40e-3),
    (6, 9.07e-3),
    (7, 2.38e-2),
    (8, 5.00e-2),
    (9, 8.96e-2),
    (10, 1.44e-1),
    (11, 2.14e-1),
    (12, 3.00e-1),
    (13, 4.00e-1),
    (14, 5.14e-1),
    (15, 6.41e-1),
    (16, 7.81e-1),
    (17, 9.31e-1),
    (18, 1.09),
    (19, 1.26),
    (20, 1.44),
    (21, 1.62),
    (22, 1.82),
    (23, 2.01),
    (24, 2.22),
    (25, 2.43),
    (26, 2.64),
    (27, 2.86),
    (28, 3.08),
    (29, 3.31),
    (30, 3.54),
];

/// Build a `Csr` (square, `n × n`) from `(row, col, value)`
/// triplets. Duplicate triplets at the same `(row, col)` are summed.
pub fn csr_from_triplets(n: usize, triplets: &[(usize, usize, f64)]) -> Csr {
    let mut counts = vec![0usize; n];
    for &(r, _, _) in triplets {
        debug_assert!(r < n, "row index {r} out of range");
        counts[r] += 1;
    }
    let mut indptr = vec![0usize; n + 1];
    for i in 0..n {
        indptr[i + 1] = indptr[i] + counts[i];
    }
    let nnz = indptr[n];
    let mut indices = vec![0u32; nnz];
    let mut data = vec![0f64; nnz];
    let mut offset = vec![0usize; n];
    for &(r, c, v) in triplets {
        debug_assert!(c < n, "col index {c} out of range");
        let pos = indptr[r] + offset[r];
        indices[pos] = c as u32;
        data[pos] = v;
        offset[r] += 1;
    }
    Csr::new_from_unsorted((n, n), indptr, indices, data)
        .map_err(|(_, _, _, e)| e)
        .expect("invalid CSR structure")
}

/// Matrix 1-norm: `max_j Σ_i |A_{ij}|` (max column sum of absolute values).
/// Used to pick the Taylor parameters `(m, s)`.
pub fn csr_one_norm(m: &Csr) -> f64 {
    if m.cols() == 0 {
        return 0.0;
    }
    let mut col_sums = vec![0f64; m.cols()];
    for (k, &col) in m.indices().iter().enumerate() {
        col_sums[col as usize] += m.data()[k].abs();
    }
    col_sums.into_iter().fold(0f64, f64::max)
}

/// `y ← A · x` (serial).
pub fn spmv_serial(m: &Csr, x: &[f64], y: &mut [f64]) {
    debug_assert_eq!(x.len(), m.cols());
    debug_assert_eq!(y.len(), m.rows());
    let indptr_raw = m.indptr();
    let indptr = indptr_raw.raw_storage();
    let indices = m.indices();
    let data = m.data();
    for (i, yi) in y.iter_mut().enumerate() {
        let mut sum = 0.0;
        for k in indptr[i]..indptr[i + 1] {
            sum += data[k] * x[indices[k] as usize];
        }
        *yi = sum;
    }
}

/// `y ← A · x` (rayon-parallel over rows).
pub fn spmv_parallel(m: &Csr, x: &[f64], y: &mut [f64]) {
    debug_assert_eq!(x.len(), m.cols());
    debug_assert_eq!(y.len(), m.rows());
    let indptr_raw = m.indptr();
    let indptr = indptr_raw.raw_storage();
    let indices = m.indices();
    let data = m.data();
    y.par_iter_mut().enumerate().for_each(|(i, yi)| {
        let mut sum = 0.0;
        for k in indptr[i]..indptr[i + 1] {
            sum += data[k] * x[indices[k] as usize];
        }
        *yi = sum;
    });
}

/// Dispatches to [`spmv_parallel`] when `nnz ≥ parallel_threshold`,
/// else [`spmv_serial`].
#[inline]
pub fn spmv(m: &Csr, x: &[f64], y: &mut [f64], parallel_threshold: usize) {
    if m.nnz() >= parallel_threshold {
        spmv_parallel(m, x, y);
    } else {
        spmv_serial(m, x, y);
    }
}

/// Options for [`expm_multiply`].
#[derive(Clone, Copy, Debug)]
pub struct ExpmOpts {
    /// Stop the inner Horner loop once successive Taylor terms drop below
    /// `tol · ‖F‖`. Default `1e-12`.
    pub tol: f64,
    /// Threshold (in matrix `nnz`) above which SpMV is parallelised. Below
    /// this rayon's task overhead beats cache-friendly serial. Default
    /// `50_000`.
    pub parallel_threshold: usize,
    /// Maximum Krylov-Taylor degree `m_star` considered by `select_ms`.
    /// `None` (default) uses the full table up to `m=30`. Capping at a
    /// smaller value (e.g. `Some(10)`) trades wall (more outer
    /// scaling-and-squaring `s` steps, so more total matvecs) for
    /// instantaneous memory (no `m`-sized stack of vectors held in the
    /// matrix-free SpMV's per-thread accumulators). Useful when `n` is
    /// large enough that `m × n × 8 B` of Krylov scratch dominates RSS.
    pub max_krylov_m: Option<u32>,
}

impl Default for ExpmOpts {
    fn default() -> Self {
        Self {
            tol: 1e-12,
            parallel_threshold: 50_000,
            max_krylov_m: None,
        }
    }
}

/// Pick `(m, s)` minimising `s·m` subject to `s ≥ ⌈t_norm / θ_m⌉, s ≥ 1`.
/// Restricted to `m ≤ 30`; for larger norms `s` simply grows linearly.
/// When `max_m` is set, only entries with `m ≤ max_m` are considered.
fn select_ms(t_norm: f64, max_m: Option<u32>) -> (u32, u32) {
    if t_norm <= 0.0 {
        return (1, 1);
    }
    let mut best_m = 1u32;
    let mut best_s = 1u32;
    let mut best_cost = u64::MAX;
    for &(m, theta) in THETA {
        if let Some(cap) = max_m {
            if m > cap {
                continue;
            }
        }
        let s_f = (t_norm / theta).ceil();
        let s = if s_f >= 1.0 { s_f as u32 } else { 1 };
        let cost = (m as u64) * (s as u64);
        if cost < best_cost {
            best_cost = cost;
            best_m = m;
            best_s = s;
        }
    }
    (best_m, best_s)
}

#[inline]
fn inf_norm(v: &[f64]) -> f64 {
    v.iter().fold(0f64, |acc, &x| acc.max(x.abs()))
}

/// Compute `exp(t · A) · b` for a real square sparse matrix `A`. Allocates
/// a fresh `Vec<f64>` of length `A.rows()`.
///
/// The algorithm is Al-Mohy & Higham (2011) Algorithm 3.2: pick
/// `(m, s)` from the precomputed `θ_m` table to minimise SpMV count,
/// then iteratively apply a degree-`m` Horner-form Taylor polynomial
/// `s` times. Inside each Horner sweep we terminate early once two
/// successive Taylor terms have norm below `opts.tol · ‖F‖`.
pub fn expm_multiply(a: &Csr, t: f64, b: &[f64], opts: ExpmOpts) -> Vec<f64> {
    assert_eq!(
        a.rows(),
        b.len(),
        "matrix size {} mismatches vector length {}",
        a.rows(),
        b.len()
    );
    let n = a.rows();
    if n == 0 {
        return Vec::new();
    }

    let t_a_one_norm = t.abs() * csr_one_norm(a);
    let (m_star, s) = select_ms(t_a_one_norm, opts.max_krylov_m);
    let scale = t / s as f64;

    let mut f = b.to_vec();
    let mut bk = b.to_vec();
    let mut work = vec![0f64; n];

    for _ in 0..s {
        let mut c1 = inf_norm(&bk);
        for i in 1..=m_star {
            // bk ← (scale / i) · A · bk
            spmv(a, &bk, &mut work, opts.parallel_threshold);
            let factor = scale / i as f64;
            for j in 0..n {
                bk[j] = factor * work[j];
            }
            // f ← f + bk
            for j in 0..n {
                f[j] += bk[j];
            }
            // Early termination: two small successive Taylor terms.
            let c2 = inf_norm(&bk);
            let f_norm = inf_norm(&f).max(f64::MIN_POSITIVE);
            if c1 + c2 <= opts.tol * f_norm {
                break;
            }
            c1 = c2;
        }
        // Outer step: `f` becomes the input for the next scaling power.
        bk.copy_from_slice(&f);
    }

    f
}

/// Complex-coefficient SpMV: `y ← A · x` where `A` is real-valued (CSR
/// of `f64`) and `x`, `y` are complex.
///
/// Implemented as two real SpMVs on the real and imaginary parts, then
/// combined. Identical numerics to a single complex SpMV by linearity.
pub fn spmv_complex(m: &Csr, x: &[Complex<f64>], y: &mut [Complex<f64>], parallel_threshold: usize) {
    debug_assert_eq!(x.len(), m.cols());
    debug_assert_eq!(y.len(), m.rows());
    let n = m.cols();
    let mut x_re: Vec<f64> = x.iter().map(|c| c.re).collect();
    let mut x_im: Vec<f64> = x.iter().map(|c| c.im).collect();
    let mut y_re = vec![0.0_f64; n];
    let mut y_im = vec![0.0_f64; n];
    spmv(m, &x_re, &mut y_re, parallel_threshold);
    spmv(m, &x_im, &mut y_im, parallel_threshold);
    for (yi, (re, im)) in y.iter_mut().zip(y_re.iter().zip(y_im.iter())) {
        *yi = Complex::new(*re, *im);
    }
    // Silence unused-mut warnings on Vec<f64> intermediates if they were
    // ever needed for further reuse.
    let _ = (&mut x_re, &mut x_im);
}

/// Compute `exp(t · A) · b` where `A` is a real CSR and `b` is a complex
/// vector. Result is a fresh `Vec<Complex<f64>>` of length `A.rows()`.
///
/// Algorithm: same Al-Mohy & Higham scaling-and-squaring as the real
/// [`expm_multiply`], with the inner SpMV replaced by [`spmv_complex`].
/// Vector arithmetic is in `Complex<f64>` throughout.
pub fn expm_multiply_complex(a: &Csr, t: f64, b: &[Complex<f64>], opts: ExpmOpts) -> Vec<Complex<f64>> {
    assert_eq!(
        a.rows(),
        b.len(),
        "matrix size {} mismatches vector length {}",
        a.rows(),
        b.len()
    );
    let n = a.rows();
    if n == 0 {
        return Vec::new();
    }

    let t_a_one_norm = t.abs() * csr_one_norm(a);
    let (m_star, s) = select_ms(t_a_one_norm, opts.max_krylov_m);
    let scale = t / s as f64;

    let mut f = b.to_vec();
    let mut bk = b.to_vec();
    let mut work = vec![Complex::new(0.0, 0.0); n];

    for _ in 0..s {
        let mut c1 = inf_norm_c(&bk);
        for i in 1..=m_star {
            spmv_complex(a, &bk, &mut work, opts.parallel_threshold);
            let factor = scale / i as f64;
            for j in 0..n {
                bk[j] = work[j] * factor;
            }
            for j in 0..n {
                f[j] += bk[j];
            }
            let c2 = inf_norm_c(&bk);
            let f_norm = inf_norm_c(&f).max(f64::MIN_POSITIVE);
            if c1 + c2 <= opts.tol * f_norm {
                break;
            }
            c1 = c2;
        }
        bk.copy_from_slice(&f);
    }

    f
}

#[inline]
fn inf_norm_c(v: &[Complex<f64>]) -> f64 {
    v.iter().fold(0f64, |acc, c| acc.max(c.norm()))
}

/// Matrix-free variant of [`expm_multiply`]. Same Al-Mohy & Higham
/// algorithm, but takes the matrix 1-norm and SpMV operation as
/// arguments instead of a [`Csr`].
///
/// `spmv_fn(x, y)` must compute `y ← M · x` where `M` is the linear
/// operator being exponentiated. `one_norm` must be the 1-norm of `M`
/// (`max_j Σ_i |M_{ij}|`); the caller computes it once before the call.
///
/// Use this when materialising `M` as a [`Csr`] would dominate memory.
/// Trade: every SpMV recomputes the action instead of streaming through
/// pre-stored values, costing more compute per matvec. Wall time is
/// (per-eval cost × `s × m_star`), vs CSR's (build cost + `s × m_star ×
/// matvec cost`), so matrix-free is wall-competitive only when CSR
/// matvec is memory-bandwidth-bound (i.e. at large `n`).
pub fn expm_multiply_mf<F>(
    n: usize,
    one_norm: f64,
    spmv_fn: F,
    t: f64,
    b: &[f64],
    opts: ExpmOpts,
) -> Vec<f64>
where
    F: Fn(&[f64], &mut [f64]),
{
    assert_eq!(
        n,
        b.len(),
        "matrix dimension {} mismatches vector length {}",
        n,
        b.len()
    );
    if n == 0 {
        return Vec::new();
    }

    let t_a_one_norm = t.abs() * one_norm;
    let (m_star, s) = select_ms(t_a_one_norm, opts.max_krylov_m);
    let scale = t / s as f64;

    let mut f = b.to_vec();
    let mut bk = b.to_vec();
    let mut work = vec![0f64; n];

    for _ in 0..s {
        let mut c1 = inf_norm(&bk);
        for i in 1..=m_star {
            spmv_fn(&bk, &mut work);
            let factor = scale / i as f64;
            for j in 0..n {
                bk[j] = factor * work[j];
            }
            for j in 0..n {
                f[j] += bk[j];
            }
            let c2 = inf_norm(&bk);
            let f_norm = inf_norm(&f).max(f64::MIN_POSITIVE);
            if c1 + c2 <= opts.tol * f_norm {
                break;
            }
            c1 = c2;
        }
        bk.copy_from_slice(&f);
    }

    f
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn diag(vals: &[f64]) -> Csr {
        let trips: Vec<_> = vals.iter().enumerate().map(|(i, &v)| (i, i, v)).collect();
        csr_from_triplets(vals.len(), &trips)
    }

    #[test]
    fn spmv_correctness() {
        // [[1, 0, 2],
        //  [0, 3, 0],
        //  [4, 0, 5]]
        let m = csr_from_triplets(
            3,
            &[
                (0, 0, 1.0),
                (0, 2, 2.0),
                (1, 1, 3.0),
                (2, 0, 4.0),
                (2, 2, 5.0),
            ],
        );
        let x = vec![1.0, 1.0, 1.0];
        let mut y = vec![0.0; 3];
        spmv_serial(&m, &x, &mut y);
        assert_eq!(y, vec![3.0, 3.0, 9.0]);
        let mut y_par = vec![0.0; 3];
        spmv_parallel(&m, &x, &mut y_par);
        assert_eq!(y_par, vec![3.0, 3.0, 9.0]);
    }

    #[test]
    fn one_norm() {
        // 1-norm = max column abs-sum. Above: col0=5, col1=3, col2=7 → 7.
        let m = csr_from_triplets(
            3,
            &[
                (0, 0, 1.0),
                (0, 2, 2.0),
                (1, 1, 3.0),
                (2, 0, 4.0),
                (2, 2, 5.0),
            ],
        );
        assert_abs_diff_eq!(csr_one_norm(&m), 7.0);
    }

    #[test]
    fn expm_zero_t_identity() {
        let m = diag(&[1.5, -2.0, 0.7]);
        let b = vec![1.0, 2.0, 3.0];
        let r = expm_multiply(&m, 0.0, &b, ExpmOpts::default());
        for i in 0..3 {
            assert_abs_diff_eq!(r[i], b[i], epsilon = 1e-14);
        }
    }

    #[test]
    fn expm_diagonal_closed_form() {
        // diag matrix → exp(t·A)·b is elementwise exp(t·a_i) * b_i.
        let m = diag(&[1.5, -2.0, 0.7]);
        let b = vec![1.0, 1.0, 1.0];
        let t = 0.5;
        let r = expm_multiply(&m, t, &b, ExpmOpts::default());
        for (i, &a_i) in [1.5, -2.0, 0.7].iter().enumerate() {
            assert_abs_diff_eq!(r[i], (t * a_i).exp(), epsilon = 1e-12);
        }
    }

    #[test]
    fn expm_skew_rotation() {
        // A = [[0, 1], [-1, 0]]: exp(tA) is rotation by t.
        // exp(t·A)·(1, 0)^T = (cos t, -sin t)^T.
        let m = csr_from_triplets(2, &[(0, 1, 1.0), (1, 0, -1.0)]);
        let b = vec![1.0, 0.0];
        let t = 0.7;
        let r = expm_multiply(&m, t, &b, ExpmOpts::default());
        assert_abs_diff_eq!(r[0], t.cos(), epsilon = 1e-12);
        assert_abs_diff_eq!(r[1], -t.sin(), epsilon = 1e-12);
    }

    #[test]
    fn expm_serial_matches_parallel() {
        // Same algorithm both paths; force one with a high threshold and
        // one with zero. Results must be bit-identical (no floating
        // reordering across rows since each row's reduction is sequential).
        let trips: Vec<_> = (0..50)
            .flat_map(|i| {
                let i: usize = i;
                let prev = if i > 0 { Some((i, i - 1, 0.3)) } else { None };
                let next = if i + 1 < 50 {
                    Some((i, i + 1, 0.5))
                } else {
                    None
                };
                [Some((i, i, -0.4 - 0.01 * i as f64)), prev, next]
                    .into_iter()
                    .flatten()
            })
            .collect();
        let m = csr_from_triplets(50, &trips);
        let b: Vec<f64> = (0..50).map(|i| 1.0 / (i + 1) as f64).collect();
        let r_serial = expm_multiply(
            &m,
            0.3,
            &b,
            ExpmOpts {
                tol: 1e-14,
                parallel_threshold: usize::MAX,
                max_krylov_m: None,
            },
        );
        let r_parallel = expm_multiply(
            &m,
            0.3,
            &b,
            ExpmOpts {
                tol: 1e-14,
                parallel_threshold: 0,
                max_krylov_m: None,
            },
        );
        for (a, b) in r_serial.iter().zip(r_parallel.iter()) {
            assert_eq!(a, b, "serial and parallel SpMV diverged");
        }
    }

    #[test]
    fn ms_selection_sane() {
        // tiny norm → small m, s = 1
        let (m, s) = select_ms(1e-9, None);
        assert!(m <= 5, "expected small m for tiny norm, got m={m}");
        assert_eq!(s, 1);

        // moderate norm → m·s should be ~10-50
        let (m, s) = select_ms(1.0, None);
        assert!((m * s) <= 50, "moderate norm cost too high: m={m} s={s}");

        // large norm → s grows
        let (_m, s) = select_ms(100.0, None);
        assert!(s >= 20, "large norm should require many steps, got s={s}");
    }
}
