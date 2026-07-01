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
//!   `quspin-expm`'s `from_parts`.

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

/// `θ_m` table for a relaxed backward-error tolerance `tol = 1e-6`, computed
/// with the same Al-Mohy & Higham (2011) construction as [`THETA`] (the
/// backward-error series `h_{m+1}(x) = log(e^{-x} T_m(x))`; validated by
/// reproducing the `u = 2^{-53}` table above to ~2 significant figures).
///
/// The predictor-corrector truncates the Pauli basis at `drop_tol` (typically
/// 1e-3), so computing `exp` to double-precision backward error (~1e-16) is
/// ~10 orders more accurate than the state it acts on. Using `tol = 1e-6`
/// (still ~1000x tighter than the truncation) admits a lower-degree Taylor
/// polynomial for the same `‖tA‖`, cutting the SpMV count (e.g. 23 -> 13 at
/// `‖tA‖ ≈ 2`) with no measurable effect on the truncated result.
pub(crate) const THETA_LOOSE: &[(u32, f64)] = &[
    (1, 2.000e-06),
    (2, 2.447e-03),
    (3, 2.863e-02),
    (4, 1.025e-01),
    (5, 2.262e-01),
    (6, 3.911e-01),
    (7, 5.866e-01),
    (8, 8.045e-01),
    (9, 1.039),
    (10, 1.285),
    (11, 1.539),
    (12, 1.801),
    (13, 2.067),
    (14, 2.337),
    (15, 2.610),
    (16, 2.885),
    (17, 3.162),
    (18, 3.441),
    (19, 3.721),
    (20, 4.001),
    (21, 4.282),
    (22, 4.564),
    (23, 4.847),
    (24, 5.129),
    (25, 5.412),
    (26, 5.696),
    (27, 5.979),
    (28, 6.263),
    (29, 6.546),
    (30, 6.830),
];

/// Pick `(m, s)` minimising `s·m` subject to `s ≥ ⌈t_norm / θ_m⌉, s ≥ 1`,
/// using the `θ_m` table `theta`. Restricted to the table's `m` range; for
/// larger norms `s` simply grows linearly. When `max_m` is set, only entries
/// with `m ≤ max_m` are considered.
fn select_ms_with(t_norm: f64, max_m: Option<u32>, theta: &[(u32, f64)]) -> (u32, u32) {
    if t_norm <= 0.0 {
        return (1, 1);
    }
    let mut best_m = 1u32;
    let mut best_s = 1u32;
    let mut best_cost = u64::MAX;
    for &(m, th) in theta {
        if let Some(cap) = max_m {
            if m > cap {
                continue;
            }
        }
        let s_f = (t_norm / th).ceil();
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

/// `(m, s)` selection at double-precision backward error ([`THETA`]).
pub(crate) fn select_ms(t_norm: f64, max_m: Option<u32>) -> (u32, u32) {
    select_ms_with(t_norm, max_m, THETA)
}

/// `(m, s)` selection at the relaxed `tol = 1e-6` backward error
/// ([`THETA_LOOSE`]) — fewer SpMVs, used on the truncated PC expm path.
pub(crate) fn select_ms_loose(t_norm: f64, max_m: Option<u32>) -> (u32, u32) {
    select_ms_with(t_norm, max_m, THETA_LOOSE)
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
