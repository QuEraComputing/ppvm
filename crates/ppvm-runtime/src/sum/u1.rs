// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::traits::*;
use crate::{config::Config, sum::PauliSum};

/// Z-magnetization-conserving two-qubit propagation on a `PauliSum`.
///
/// The current implementation expresses each gate as a composition of
/// the existing `rxx` / `ryy` / `rzz` primitives, which is mathematically
/// equivalent to a fully fused single-pass implementation because the
/// XX, YY, and ZZ generators all pairwise commute. A future revision is
/// free to specialize this impl with a fused single `map_insert_multiple`
/// pass — the public signature is stable.
impl<T: Config> U1Conserving<T> for PauliSum<T>
where
    PauliSum<T>: RotationTwo<T>,
    T::Coeff: Clone,
{
    fn exchange(&mut self, a: usize, b: usize, theta: impl Into<T::Coeff>) {
        let theta: T::Coeff = theta.into();
        self.rxx(a, b, theta.clone());
        self.ryy(a, b, theta);
    }

    fn xyzz(
        &mut self,
        a: usize,
        b: usize,
        theta_xy: impl Into<T::Coeff>,
        theta_zz: impl Into<T::Coeff>,
    ) {
        self.exchange(a, b, theta_xy);
        self.rzz(a, b, theta_zz);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::fxhash::ByteF64;
    use crate::strategy::CoefficientThreshold;
    use std::f64::consts::PI;

    // Use CoefficientThreshold so the propagated maps drop the cos(θ)·0 and
    // sin(θ)·0 ghost entries that rotate_2 leaves behind at θ = 0 or when a
    // sub-branch happens to vanish. Without this, `assert_eq!` would compare
    // maps that differ only in formally-zero terms.
    type Cfg = ByteF64<2, CoefficientThreshold>;

    fn ps(n: usize, term: &str, coeff: f64) -> PauliSum<Cfg> {
        let mut s: PauliSum<Cfg> = PauliSum::builder()
            .n_qubits(n)
            .strategy(CoefficientThreshold(1e-12))
            .build();
        s += (term, coeff);
        s
    }

    /// `exchange(a, b, 0)` is the identity.
    #[test]
    fn exchange_zero_angle_is_identity() {
        for term in ["II", "IZ", "ZI", "ZZ", "XY", "YX", "XX", "YY"] {
            let expect = ps(2, term, 1.0);
            let mut got = ps(2, term, 1.0);
            got.exchange(0, 1, 0.0);
            got.truncate();
            assert_eq!(got, expect, "exchange(0) altered {}", term);
        }
    }

    /// `exchange(a, b, θ)` matches `rxx(θ)` followed by `ryy(θ)`.
    #[test]
    fn exchange_matches_rxx_then_ryy() {
        let thetas = [0.0_f64, 0.123, 0.7, -0.4, PI / 3.0, PI / 2.0];
        let terms = [
            "II", "IZ", "ZI", "ZZ", "XY", "YX", "XX", "YY", "IX", "XI", "IY", "YI", "XZ", "ZX",
            "YZ", "ZY",
        ];
        for theta in thetas {
            for term in terms {
                let mut fused = ps(2, term, 1.0);
                let mut composed = ps(2, term, 1.0);
                fused.exchange(0, 1, theta);
                composed.rxx(0, 1, theta);
                composed.ryy(0, 1, theta);
                fused.truncate();
                composed.truncate();
                assert_eq!(
                    fused, composed,
                    "exchange disagrees on {} at θ={}",
                    term, theta
                );
            }
        }
    }

    /// `exchange` preserves the identity term: `I ↦ I` with the same coefficient.
    #[test]
    fn exchange_preserves_identity() {
        let expect = ps(2, "II", 0.7);
        let mut s = ps(2, "II", 0.7);
        s.exchange(0, 1, 0.42);
        s.truncate();
        assert_eq!(s, expect);
    }

    /// `Z_a + Z_b` commutes with `X_a X_b + Y_a Y_b`, so total magnetization
    /// must propagate trivially through `exchange`.
    #[test]
    fn exchange_preserves_total_z() {
        let mut expect: PauliSum<Cfg> = PauliSum::builder()
            .n_qubits(2)
            .strategy(CoefficientThreshold(1e-12))
            .build();
        expect += ("IZ", 1.0);
        expect += ("ZI", 1.0);
        let mut s = expect.clone();
        s.exchange(0, 1, 0.37);
        s.exchange(0, 1, -1.1);
        s.truncate();
        assert_eq!(s, expect, "total Z was perturbed by exchange");
    }

    /// `xyzz(θ_xy, θ_zz)` matches `exchange(θ_xy)` followed by `rzz(θ_zz)`.
    #[test]
    fn xyzz_matches_exchange_then_rzz() {
        let pairs = [
            (0.0_f64, 0.0_f64),
            (0.3, 0.1),
            (-0.4, 0.7),
            (PI / 5.0, -PI / 7.0),
        ];
        let terms = ["IZ", "ZI", "XY", "YX", "XX", "YY", "ZZ", "IX", "YI"];
        for (theta_xy, theta_zz) in pairs {
            for term in terms {
                let mut combined = ps(2, term, 1.0);
                let mut stepwise = ps(2, term, 1.0);
                combined.xyzz(0, 1, theta_xy, theta_zz);
                stepwise.exchange(0, 1, theta_xy);
                stepwise.rzz(0, 1, theta_zz);
                combined.truncate();
                stepwise.truncate();
                assert_eq!(
                    combined, stepwise,
                    "xyzz disagrees on {} at θ_xy={}, θ_zz={}",
                    term, theta_xy, theta_zz
                );
            }
        }
    }

    /// `xyzz` with both angles zero is the identity.
    #[test]
    fn xyzz_zero_angles_is_identity() {
        for term in ["II", "IZ", "ZI", "XY", "YX", "ZZ"] {
            let expect = ps(2, term, 1.0);
            let mut got = ps(2, term, 1.0);
            got.xyzz(0, 1, 0.0, 0.0);
            got.truncate();
            assert_eq!(got, expect, "xyzz(0,0) altered {}", term);
        }
    }

    /// The magnetization difference and the "current" operator
    /// `Y_a X_b − X_a Y_b` form a closed two-dimensional U(1) sector;
    /// `exchange(θ)` rotates between them with the doubled-angle
    /// `cos(2θ)` / `sin(2θ)` (the ppvm angle convention is
    /// `exp(−i θ/2 G)`, so two stacked rotations contribute a 2θ).
    #[test]
    fn exchange_rotates_magnetization_difference_into_current() {
        let theta = 0.3_f64;
        let mut s: PauliSum<Cfg> = PauliSum::builder()
            .n_qubits(2)
            .strategy(CoefficientThreshold(1e-12))
            .build();
        s += ("ZI", 0.5);
        s += ("IZ", -0.5);
        s.exchange(0, 1, theta);
        s.truncate();

        let c = (2.0 * theta).cos();
        let si = (2.0 * theta).sin();
        let mut want: PauliSum<Cfg> = PauliSum::builder()
            .n_qubits(2)
            .strategy(CoefficientThreshold(1e-12))
            .build();
        want += ("ZI", 0.5 * c);
        want += ("IZ", -0.5 * c);
        want += ("YX", 0.5 * si);
        want += ("XY", -0.5 * si);

        // Compare via terms() under a tolerance — the maps store f64s, so
        // structural equality is too strict here.
        let mut got_terms: Vec<(String, f64)> =
            s.data().iter().map(|(k, v)| (k.to_string(), *v)).collect();
        let mut want_terms: Vec<(String, f64)> = want
            .data()
            .iter()
            .map(|(k, v)| (k.to_string(), *v))
            .collect();
        got_terms.sort_by(|a, b| a.0.cmp(&b.0));
        want_terms.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(got_terms.len(), want_terms.len());
        for ((kg, vg), (kw, vw)) in got_terms.iter().zip(want_terms.iter()) {
            assert_eq!(kg, kw);
            assert!(
                (vg - vw).abs() < 1e-12,
                "coeff for {} differs: got {}, want {}",
                kg,
                vg,
                vw
            );
        }
    }

    // ===========================================================================
    // Truncation-robust conservation hardening tests
    // ===========================================================================
    //
    // These four tests pin down the realistic floating-point-precision guarantee
    // that `exchange` / `xyzz` / `rzz` provide on observables built from `{I, Z}`
    // Pauli strings (e.g. `Σ_i Z_i`, `Σ Z_i Z_j`). The conserved-sector
    // coefficients recover to their starting values modulo per-gate ε (≈ 1e-15),
    // so we compare under a `1e-10` tolerance — three orders of magnitude above
    // accumulated drift for a short circuit, ten orders of magnitude below the
    // conserved-coefficient magnitude (1.0).

    use crate::strategy::{CombinedStrategy, MaxPauliWeight};

    /// Sort-and-compare PauliSums under a per-term tolerance.
    fn assert_close<T>(left: &PauliSum<T>, right: &PauliSum<T>, tol: f64)
    where
        T: crate::config::Config<Coeff = f64>,
        T::Map: for<'a> crate::traits::ACMapIter<
                'a,
                Item = (&'a <T as crate::config::Config>::PauliWordType, &'a f64),
            >,
    {
        let mut left_terms: Vec<(String, f64)> = left
            .data()
            .iter()
            .map(|(k, v)| (k.to_string(), *v))
            .collect();
        let mut right_terms: Vec<(String, f64)> = right
            .data()
            .iter()
            .map(|(k, v)| (k.to_string(), *v))
            .collect();
        left_terms.sort_by(|a, b| a.0.cmp(&b.0));
        right_terms.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(
            left_terms.len(),
            right_terms.len(),
            "term count differs: left={:?}, right={:?}",
            left_terms,
            right_terms
        );
        for ((kl, vl), (kr, vr)) in left_terms.iter().zip(right_terms.iter()) {
            assert_eq!(kl, kr, "key mismatch: {} vs {}", kl, kr);
            assert!(
                (vl - vr).abs() < tol,
                "coeff for {} differs by more than {}: got {}, want {}",
                kl,
                tol,
                vl,
                vr
            );
        }
    }

    /// `Σ Z_i` is preserved across a chain of `exchange` calls when truncation
    /// uses an aggressive `CoefficientThreshold` (well below the conserved
    /// coefficient magnitude). The transient cross terms produced internally
    /// by each gate cancel back to ε before truncation, so the cutoff (0.5)
    /// drops only the ε residues and leaves `Σ Z_i` intact.
    #[test]
    fn exchange_preserves_total_z_under_coefficient_truncation() {
        type C = ByteF64<1, CoefficientThreshold>;
        let mut s: PauliSum<C> = PauliSum::builder()
            .n_qubits(4)
            .strategy(CoefficientThreshold(0.5))
            .build();
        let mut expect: PauliSum<C> = PauliSum::builder()
            .n_qubits(4)
            .strategy(CoefficientThreshold(0.5))
            .build();
        for term in ["ZIII", "IZII", "IIZI", "IIIZ"] {
            s += (term, 1.0);
            expect += (term, 1.0);
        }
        for (a, b) in [(0, 1), (1, 2), (2, 3)] {
            s.exchange(a, b, 0.37);
            s.truncate();
        }
        assert_close(&s, &expect, 1e-10);
    }

    /// Same conservation under `MaxPauliWeight(1)`: the cross terms have weight
    /// 2 and are dropped by the discrete weight check, independent of any
    /// floating-point sensitivity in the threshold comparison.
    #[test]
    fn exchange_preserves_total_z_under_weight_truncation() {
        type C = ByteF64<1, CombinedStrategy<CoefficientThreshold, MaxPauliWeight>>;
        let strat = CombinedStrategy(CoefficientThreshold(1e-12), MaxPauliWeight(1));
        let mut s: PauliSum<C> = PauliSum::builder().n_qubits(4).strategy(strat).build();
        let mut expect: PauliSum<C> = PauliSum::builder().n_qubits(4).strategy(strat).build();
        for term in ["ZIII", "IZII", "IIZI", "IIIZ"] {
            s += (term, 1.0);
            expect += (term, 1.0);
        }
        for (a, b) in [(0, 1), (1, 2), (2, 3)] {
            s.exchange(a, b, 0.41);
            s.truncate();
        }
        assert_close(&s, &expect, 1e-10);
    }

    /// `Σ_{i<j} Z_i Z_j` (sum over *all* pairs, not just nearest neighbours)
    /// commutes with `XX+YY+ZZ` on every edge, so `xyzz` on any edge
    /// preserves the full sum. The nearest-neighbour sum is *not*
    /// conserved: the cross terms that exchange generates from
    /// `Z_aZ_c + Z_bZ_c` cancel only when the matching `Z_aZ_c` and
    /// `Z_bZ_c` pair is present for every "third site" `c`. This
    /// catches any regression in `rzz` semantics under aggressive
    /// truncation and is a tighter check on `xyzz`'s ZZ piece than the
    /// single-Z tests above.
    #[test]
    fn zz_correlation_preserves_under_truncation() {
        type C = ByteF64<1, CombinedStrategy<CoefficientThreshold, MaxPauliWeight>>;
        let strat = CombinedStrategy(CoefficientThreshold(0.5), MaxPauliWeight(2));
        let mut s: PauliSum<C> = PauliSum::builder().n_qubits(4).strategy(strat).build();
        let mut expect: PauliSum<C> = PauliSum::builder().n_qubits(4).strategy(strat).build();
        // All C(4, 2) = 6 unordered pairs.
        for term in ["ZZII", "ZIZI", "ZIIZ", "IZZI", "IZIZ", "IIZZ"] {
            s += (term, 1.0);
            expect += (term, 1.0);
        }
        for (a, b) in [(0, 1), (1, 2), (2, 3)] {
            s.xyzz(a, b, 0.21, 0.07);
            s.truncate();
        }
        assert_close(&s, &expect, 1e-10);
    }

    /// One full Trotter cycle (xyzz on every edge + rz on every site) on
    /// `Σ Z_i`, repeated three times, under combined coefficient + weight
    /// truncation. The Trotter helper itself lives in Python, so we
    /// replicate its gate order manually here.
    #[test]
    fn u1_trotter_total_z_under_combined_truncation() {
        type C = ByteF64<1, CombinedStrategy<CoefficientThreshold, MaxPauliWeight>>;
        let strat = CombinedStrategy(CoefficientThreshold(0.1), MaxPauliWeight(1));
        let mut s: PauliSum<C> = PauliSum::builder().n_qubits(5).strategy(strat).build();
        let mut expect: PauliSum<C> = PauliSum::builder().n_qubits(5).strategy(strat).build();
        for term in ["ZIIII", "IZIII", "IIZII", "IIIZI", "IIIIZ"] {
            s += (term, 1.0);
            expect += (term, 1.0);
        }
        let theta_xy = 0.15;
        let theta_zz = 0.05;
        let h = 0.03;
        for _ in 0..3 {
            for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 4)] {
                s.xyzz(a, b, theta_xy, theta_zz);
                s.truncate();
            }
            for site in 0..5 {
                s.rz(site, h);
                s.truncate();
            }
        }
        assert_close(&s, &expect, 1e-10);
    }
}
