// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Support pieces for the `quspin-expm`-backed `exp(t·A)·b` engine.
//!
//! The hand-rolled scaling-and-squaring engine that once lived here has been
//! retired in favour of the external `quspin-expm` crate (driven from
//! [`crate::mf_expm`]). What remains is the shared scaffolding the new path
//! still relies on:
//!
//! - the `(m, s)` selection table [`THETA`] / [`select_ms`] from Al-Mohy &
//!   Higham (2011), used to pick the Taylor partition handed to
//!   `quspin-expm`'s `from_parts`;
//! - the complex CSR type [`CsrCx`] and its 1-norm / SpMV
//!   ([`csr_cx_one_norm`], [`spmv_cx`]), used by the orbit-rep evolution
//!   path (whose momentum-character phases make the matrix complex).

use num::Complex;
use rayon::prelude::*;
use sprs::CsMatI;

/// Complex-coefficient CSR matrix. Used by the orbit-rep evolution
/// path, where matrix elements pick up `χ_k(g)` phases from the
/// translation generator and become complex.
pub type CsrCx = CsMatI<Complex<f64>, u32, usize>;

/// `θ_m` table from Al-Mohy & Higham (2011), Table A.3, for double
/// precision (unit roundoff `u = 2^{-53}`).
///
/// `θ_m` bounds `‖A‖₁` such that the degree-`m` Taylor polynomial
/// approximates `exp(A)` to within `u`. We pick `(m, s)` with
/// `s ≥ ⌈‖tA‖₁ / θ_m⌉` and minimise `m·s` (total SpMV count).
pub(crate) const THETA: &[(u32, f64)] = &[
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

/// Pick `(m, s)` minimising `s·m` subject to `s ≥ ⌈t_norm / θ_m⌉, s ≥ 1`.
/// Restricted to `m ≤ 30`; for larger norms `s` simply grows linearly.
/// When `max_m` is set, only entries with `m ≤ max_m` are considered.
pub(crate) fn select_ms(t_norm: f64, max_m: Option<u32>) -> (u32, u32) {
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

/// Matrix 1-norm `max_j Σ_i |A_{ij}|` for a complex CSR.
pub fn csr_cx_one_norm(m: &CsrCx) -> f64 {
    if m.cols() == 0 {
        return 0.0;
    }
    let mut col_sums = vec![0f64; m.cols()];
    for (k, &col) in m.indices().iter().enumerate() {
        col_sums[col as usize] += m.data()[k].norm();
    }
    col_sums.into_iter().fold(0f64, f64::max)
}

/// `y ← A · x` for a complex CSR + complex vector. Serial.
pub fn spmv_cx_serial(m: &CsrCx, x: &[Complex<f64>], y: &mut [Complex<f64>]) {
    debug_assert_eq!(x.len(), m.cols());
    debug_assert_eq!(y.len(), m.rows());
    let indptr_raw = m.indptr();
    let indptr = indptr_raw.raw_storage();
    let indices = m.indices();
    let data = m.data();
    for (i, yi) in y.iter_mut().enumerate() {
        let mut sum = Complex::new(0.0, 0.0);
        for k in indptr[i]..indptr[i + 1] {
            sum += data[k] * x[indices[k] as usize];
        }
        *yi = sum;
    }
}

/// `y ← A · x` for a complex CSR + complex vector. Rayon-parallel over rows.
pub fn spmv_cx_parallel(m: &CsrCx, x: &[Complex<f64>], y: &mut [Complex<f64>]) {
    debug_assert_eq!(x.len(), m.cols());
    debug_assert_eq!(y.len(), m.rows());
    let indptr_raw = m.indptr();
    let indptr = indptr_raw.raw_storage();
    let indices = m.indices();
    let data = m.data();
    y.par_iter_mut().enumerate().for_each(|(i, yi)| {
        let mut sum = Complex::new(0.0, 0.0);
        for k in indptr[i]..indptr[i + 1] {
            sum += data[k] * x[indices[k] as usize];
        }
        *yi = sum;
    });
}

#[inline]
pub fn spmv_cx(m: &CsrCx, x: &[Complex<f64>], y: &mut [Complex<f64>], parallel_threshold: usize) {
    if m.nnz() >= parallel_threshold {
        spmv_cx_parallel(m, x, y);
    } else {
        spmv_cx_serial(m, x, y);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
