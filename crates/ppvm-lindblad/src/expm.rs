// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Action of the matrix exponential, `y = exp(t·A)·b`, for real sparse `A`.
//!
//! Implements Al-Mohy & Higham (2011), Algorithm 3.2: a Horner-form
//! degree-`m` Taylor polynomial evaluated `s` times (scaling-and-squaring),
//! with `(m, s)` chosen from a precomputed table that minimises total SpMV
//! count subject to a double-precision truncation bound.
//!
//! The hot path is repeated sparse-matrix × dense-vector products. SpMV is
//! parallelised over rows with `rayon` when the matrix has more than
//! [`ExpmOpts::parallel_threshold`] nonzeros; below that, the task-spawn
//! overhead beats the cache-friendliness of the serial loop, so we stay
//! single-threaded. Both branches share the same CSR storage and produce
//! bit-identical output.

use rayon::prelude::*;

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

/// Compressed Sparse Row matrix (square, `n × n`).
#[derive(Clone, Debug)]
pub struct CsrMatrix {
    pub n: usize,
    /// Length `n+1`; row `i` spans `[row_ptr[i], row_ptr[i+1])` in
    /// `col_idx` and `values`.
    pub row_ptr: Vec<usize>,
    pub col_idx: Vec<u32>,
    pub values: Vec<f64>,
}

impl CsrMatrix {
    /// Number of structural nonzeros.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Build CSR from a list of `(row, col, value)` triplets. Duplicate
    /// triplets at the same `(row, col)` are summed.
    pub fn from_triplets(n: usize, triplets: &[(usize, usize, f64)]) -> Self {
        let mut counts = vec![0usize; n];
        for &(r, _, _) in triplets {
            debug_assert!(r < n, "row index {r} out of range");
            counts[r] += 1;
        }
        let mut row_ptr = vec![0usize; n + 1];
        for i in 0..n {
            row_ptr[i + 1] = row_ptr[i] + counts[i];
        }
        let mut col_idx = vec![0u32; triplets.len()];
        let mut values = vec![0f64; triplets.len()];
        let mut offset = vec![0usize; n];
        for &(r, c, v) in triplets {
            debug_assert!(c < n, "col index {c} out of range");
            let pos = row_ptr[r] + offset[r];
            col_idx[pos] = c as u32;
            values[pos] = v;
            offset[r] += 1;
        }
        Self {
            n,
            row_ptr,
            col_idx,
            values,
        }
    }

    /// Matrix 1-norm: `max_j Σ_i |A_{ij}|` (max column sum of absolute
    /// values). Used to pick the Taylor parameters `(m, s)`.
    pub fn one_norm(&self) -> f64 {
        if self.n == 0 {
            return 0.0;
        }
        let mut col_sums = vec![0f64; self.n];
        for k in 0..self.values.len() {
            col_sums[self.col_idx[k] as usize] += self.values[k].abs();
        }
        col_sums.into_iter().fold(0f64, f64::max)
    }

    /// `y ← A · x` (serial).
    pub fn spmv_serial(&self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.n);
        debug_assert_eq!(y.len(), self.n);
        for (i, yi) in y.iter_mut().enumerate() {
            let mut sum = 0.0;
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                sum += self.values[k] * x[self.col_idx[k] as usize];
            }
            *yi = sum;
        }
    }

    /// `y ← A · x` (rayon-parallel over rows).
    pub fn spmv_parallel(&self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.n);
        debug_assert_eq!(y.len(), self.n);
        y.par_iter_mut().enumerate().for_each(|(i, yi)| {
            let mut sum = 0.0;
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                sum += self.values[k] * x[self.col_idx[k] as usize];
            }
            *yi = sum;
        });
    }

    /// Dispatches to [`Self::spmv_parallel`] when `nnz ≥ parallel_threshold`,
    /// else [`Self::spmv_serial`].
    #[inline]
    pub fn spmv(&self, x: &[f64], y: &mut [f64], parallel_threshold: usize) {
        if self.nnz() >= parallel_threshold {
            self.spmv_parallel(x, y);
        } else {
            self.spmv_serial(x, y);
        }
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
}

impl Default for ExpmOpts {
    fn default() -> Self {
        Self {
            tol: 1e-12,
            parallel_threshold: 50_000,
        }
    }
}

/// Pick `(m, s)` minimising `s·m` subject to `s ≥ ⌈t_norm / θ_m⌉, s ≥ 1`.
/// Restricted to `m ≤ 30`; for larger norms `s` simply grows linearly.
fn select_ms(t_norm: f64) -> (u32, u32) {
    if t_norm <= 0.0 {
        return (1, 1);
    }
    let mut best_m = 1u32;
    let mut best_s = 1u32;
    let mut best_cost = u64::MAX;
    for &(m, theta) in THETA {
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

/// Compute `exp(t · A) · b` for a real square sparse matrix `A` in CSR
/// form. Allocates a fresh `Vec<f64>` of length `A.n`.
///
/// The algorithm is Al-Mohy & Higham (2011) Algorithm 3.2: pick
/// `(m, s)` from the precomputed `θ_m` table to minimise SpMV count,
/// then iteratively apply a degree-`m` Horner-form Taylor polynomial
/// `s` times. Inside each Horner sweep we terminate early once two
/// successive Taylor terms have norm below `opts.tol · ‖F‖`.
pub fn expm_multiply(a: &CsrMatrix, t: f64, b: &[f64], opts: ExpmOpts) -> Vec<f64> {
    assert_eq!(
        a.n,
        b.len(),
        "matrix size {} mismatches vector length {}",
        a.n,
        b.len()
    );
    let n = a.n;
    if n == 0 {
        return Vec::new();
    }

    let t_a_one_norm = t.abs() * a.one_norm();
    let (m_star, s) = select_ms(t_a_one_norm);
    let scale = t / s as f64;

    let mut f = b.to_vec();
    let mut bk = b.to_vec();
    let mut work = vec![0f64; n];

    for _ in 0..s {
        let mut c1 = inf_norm(&bk);
        for i in 1..=m_star {
            // bk ← (scale / i) · A · bk
            a.spmv(&bk, &mut work, opts.parallel_threshold);
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn diag(vals: &[f64]) -> CsrMatrix {
        let trips: Vec<_> = vals.iter().enumerate().map(|(i, &v)| (i, i, v)).collect();
        CsrMatrix::from_triplets(vals.len(), &trips)
    }

    #[test]
    fn spmv_correctness() {
        // [[1, 0, 2],
        //  [0, 3, 0],
        //  [4, 0, 5]]
        let m = CsrMatrix::from_triplets(
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
        m.spmv_serial(&x, &mut y);
        assert_eq!(y, vec![3.0, 3.0, 9.0]);
        let mut y_par = vec![0.0; 3];
        m.spmv_parallel(&x, &mut y_par);
        assert_eq!(y_par, vec![3.0, 3.0, 9.0]);
    }

    #[test]
    fn one_norm() {
        // 1-norm = max column abs-sum. Above: col0=5, col1=3, col2=7 → 7.
        let m = CsrMatrix::from_triplets(
            3,
            &[
                (0, 0, 1.0),
                (0, 2, 2.0),
                (1, 1, 3.0),
                (2, 0, 4.0),
                (2, 2, 5.0),
            ],
        );
        assert_abs_diff_eq!(m.one_norm(), 7.0);
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
        let m = CsrMatrix::from_triplets(2, &[(0, 1, 1.0), (1, 0, -1.0)]);
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
        let m = CsrMatrix::from_triplets(50, &trips);
        let b: Vec<f64> = (0..50).map(|i| 1.0 / (i + 1) as f64).collect();
        let r_serial = expm_multiply(
            &m,
            0.3,
            &b,
            ExpmOpts {
                tol: 1e-14,
                parallel_threshold: usize::MAX,
            },
        );
        let r_parallel = expm_multiply(
            &m,
            0.3,
            &b,
            ExpmOpts {
                tol: 1e-14,
                parallel_threshold: 0,
            },
        );
        for (a, b) in r_serial.iter().zip(r_parallel.iter()) {
            assert_eq!(a, b, "serial and parallel SpMV diverged");
        }
    }

    #[test]
    fn ms_selection_sane() {
        // tiny norm → small m, s = 1
        let (m, s) = select_ms(1e-9);
        assert!(m <= 5, "expected small m for tiny norm, got m={m}");
        assert_eq!(s, 1);

        // moderate norm → m·s should be ~10-50
        let (m, s) = select_ms(1.0);
        assert!((m * s) <= 50, "moderate norm cost too high: m={m} s={s}");

        // large norm → s grows
        let (_m, s) = select_ms(100.0);
        assert!(s >= 20, "large norm should require many steps, got s={s}");
    }
}
