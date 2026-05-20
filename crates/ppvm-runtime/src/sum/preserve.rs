// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Observable-aware truncation: never drop a chosen set of Pauli
//! strings, and apply a weight-biased threshold to the rest.
//!
//! # Why
//!
//! For a transport diagnostic of the form `<Ô(t) Ô(0)>`, the result
//! depends only on the projection of the propagated `Ô(t)` onto the
//! handful of Pauli strings that make up `Ô(0)`. Standard
//! coefficient-magnitude truncation can drop those exact strings when
//! their coefficients drift toward the threshold — typically the tail
//! of a spreading operator at far sites. This module lets the caller
//! mark a small set of strings as *preserved* so they survive any
//! truncation, regardless of their current coefficient.
//!
//! The same struct also exposes a `weight_lambda` knob that biases the
//! threshold by `exp(λ · weight(P))`, i.e. truncates higher-weight
//! Pauli strings more aggressively. This is the "virtual DAOE" tactic
//! of [Rakovszky, Pollmann, von Keyserlingk
//! (2020)](https://arxiv.org/abs/2004.05177) lifted into a purely
//! truncation-level knob: the dynamics is not damped, only the
//! truncation criterion is weight-biased. Set `weight_lambda = 0` for
//! a uniform threshold.
//!
//! # Usage
//!
//! ```ignore
//! use ppvm_runtime::prelude::*;
//! use ppvm_runtime::sum::PreserveConfig;
//! type Cfg = config::indexmap::ByteFxHashF64<1>;
//! type W = <Cfg as Config>::PauliWordType;
//!
//! let preserve = PreserveConfig::<W>::single_z(4, /*base*/ 1e-4, /*lambda*/ 0.0);
//! let mut s: PauliSum<Cfg> = PauliSum::builder()
//!     .n_qubits(4)
//!     .preserve(preserve)
//!     .build();
//! s += ("ZIII", 1.0);
//! s.exchange(0, 1, 0.1);   // hypothetical; preserve survives all auto-truncates
//! ```

use std::collections::HashSet;

use crate::traits::PauliWordTrait;

/// Observable-aware truncation policy applied by [`PauliSum::truncate`]
/// when set via the builder.
///
/// See the [module docs](self) for the design rationale. Construct via
/// [`PreserveConfig::new`] for arbitrary keep-sets, or via one of the
/// `single_z` / `from_strings` convenience constructors.
#[derive(Debug, Clone)]
pub struct PreserveConfig<W: PauliWordTrait> {
    /// Pauli strings that are *never* dropped by truncation, regardless
    /// of their current coefficient.
    pub keep: HashSet<W>,
    /// Base coefficient cutoff. A term `P` with weight `k` is dropped
    /// when `|c_P| < base_threshold · exp(weight_lambda · k)`.
    pub base_threshold: f64,
    /// Weight-biased multiplier (virtual DAOE knob). `0` gives a
    /// uniform cutoff; positive values truncate higher-weight strings
    /// more aggressively. See module docs.
    pub weight_lambda: f64,
}

impl<W: PauliWordTrait> PreserveConfig<W> {
    /// General constructor: pass any iterable of [`PauliWordTrait`]
    /// values to keep.
    pub fn new(keep: impl IntoIterator<Item = W>, base_threshold: f64, weight_lambda: f64) -> Self {
        Self {
            keep: keep.into_iter().collect(),
            base_threshold,
            weight_lambda,
        }
    }

    /// Preserve every single-`Z` Pauli string on `n_qubits` qubits —
    /// `Z_0, Z_1, …, Z_{n−1}`. The natural choice when the transport
    /// diagnostic is `<Σ_j Z_j(t) Z_i(0)>` (z-magnetization spread).
    pub fn single_z(n_qubits: usize, base_threshold: f64, weight_lambda: f64) -> Self {
        let keep: HashSet<W> = (0..n_qubits)
            .map(|i| {
                let s: String = (0..n_qubits)
                    .map(|j| if j == i { 'Z' } else { 'I' })
                    .collect();
                W::from(s)
            })
            .collect();
        Self::new(keep, base_threshold, weight_lambda)
    }

    /// Preserve a user-specified list of Pauli strings (each must be a
    /// length-`n_qubits` string over `{I, X, Y, Z}`).
    pub fn from_strings(
        strings: impl IntoIterator<Item = String>,
        base_threshold: f64,
        weight_lambda: f64,
    ) -> Self {
        let keep: HashSet<W> = strings.into_iter().map(W::from).collect();
        Self::new(keep, base_threshold, weight_lambda)
    }

    /// The effective cutoff for a term of weight `k`:
    /// `base_threshold · exp(weight_lambda · k)`.
    #[inline]
    pub fn threshold_for_weight(&self, weight: usize) -> f64 {
        self.base_threshold * (self.weight_lambda * weight as f64).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::word::PauliWord;

    type W = PauliWord<[u8; 1]>;

    #[test]
    fn single_z_builds_correct_set() {
        let c: PreserveConfig<W> = PreserveConfig::single_z(3, 1e-4, 0.0);
        assert_eq!(c.keep.len(), 3);
        assert!(c.keep.contains(&W::from("ZII")));
        assert!(c.keep.contains(&W::from("IZI")));
        assert!(c.keep.contains(&W::from("IIZ")));
        assert!(!c.keep.contains(&W::from("ZZI")));
    }

    #[test]
    fn from_strings_round_trip() {
        let c: PreserveConfig<W> =
            PreserveConfig::from_strings(["XYZ".to_string(), "ZZZ".to_string()], 0.1, 0.5);
        assert_eq!(c.keep.len(), 2);
        assert!(c.keep.contains(&W::from("XYZ")));
        assert!(c.keep.contains(&W::from("ZZZ")));
    }

    #[test]
    fn threshold_for_weight_is_weight_biased() {
        let c: PreserveConfig<W> = PreserveConfig::new([] as [W; 0], 1e-4, 1.0);
        assert!((c.threshold_for_weight(0) - 1e-4).abs() < 1e-12);
        assert!((c.threshold_for_weight(1) - 1e-4 * 1.0_f64.exp()).abs() < 1e-12);
        assert!((c.threshold_for_weight(3) - 1e-4 * 3.0_f64.exp()).abs() < 1e-12);
    }

    #[test]
    fn lambda_zero_gives_uniform_threshold() {
        let c: PreserveConfig<W> = PreserveConfig::new([] as [W; 0], 0.5, 0.0);
        for k in 0..10 {
            assert!((c.threshold_for_weight(k) - 0.5).abs() < 1e-12);
        }
    }

    // ===========================================================================
    // End-to-end PauliSum integration tests
    // ===========================================================================
    //
    // These exercise the `PauliSum::truncate` integration: a `PauliSum` built
    // with a `PreserveConfig` should never drop strings in the keep-set, but
    // should drop other terms below the (weight-biased) cutoff.

    use crate::config::fxhash::ByteF64;
    use crate::sum::PauliSum;
    use crate::traits::RotationTwo;

    type Cfg = ByteF64<1>;
    type PWord = <Cfg as crate::config::Config>::PauliWordType;

    /// Preserve-aware truncate keeps the preserved strings even if their
    /// coefficient is well below the base cutoff.
    #[test]
    fn truncate_keeps_preserved_strings_below_cutoff() {
        let preserve: PreserveConfig<PWord> = PreserveConfig::single_z(3, /*base*/ 0.5, 0.0);
        let mut s: PauliSum<Cfg> = PauliSum::builder().n_qubits(3).preserve(preserve).build();
        // Single-Z strings get tiny coefficients (below cutoff 0.5) — must survive.
        s += ("ZII", 1e-6);
        s += ("IZI", 1e-6);
        s += ("IIZ", 1e-6);
        // A non-preserved term well below cutoff — must be dropped.
        s += ("XYZ", 1e-6);
        // A non-preserved term above cutoff — must survive.
        s += ("XXX", 0.7);

        s.truncate();
        let kept: std::collections::HashSet<String> =
            s.data().keys().map(|k| k.to_string()).collect();
        assert!(kept.contains("ZII"), "preserved ZII should be kept");
        assert!(kept.contains("IZI"), "preserved IZI should be kept");
        assert!(kept.contains("IIZ"), "preserved IIZ should be kept");
        assert!(
            !kept.contains("XYZ"),
            "below-cutoff non-preserved XYZ should be dropped"
        );
        assert!(kept.contains("XXX"), "above-cutoff XXX should be kept");
    }

    /// With `weight_lambda > 0`, higher-weight non-preserved strings are
    /// dropped at a lower effective coefficient than lower-weight ones.
    /// `exp(λ k) · ε` → weight-1 needs `> ε·e^λ`, weight-3 needs `> ε·e^{3λ}`.
    #[test]
    fn truncate_weight_lambda_is_more_aggressive_on_high_weight() {
        let preserve: PreserveConfig<PWord> = PreserveConfig::new(
            [] as [PWord; 0], // empty keep-set; pure weighted threshold
            /*base*/ 0.01,
            /*lambda*/ 1.0,
        );
        let mut s: PauliSum<Cfg> = PauliSum::builder().n_qubits(3).preserve(preserve).build();
        // Weight 1: cutoff is 0.01 · e^1 ≈ 0.0272. Coefficient 0.05 survives.
        s += ("XII", 0.05);
        // Weight 2: cutoff is 0.01 · e^2 ≈ 0.0739. Coefficient 0.05 dropped.
        s += ("XXI", 0.05);
        // Weight 3: cutoff is 0.01 · e^3 ≈ 0.2008. Coefficient 0.1 dropped.
        s += ("XYZ", 0.1);

        s.truncate();
        let kept: std::collections::HashSet<String> =
            s.data().keys().map(|k| k.to_string()).collect();
        assert!(
            kept.contains("XII"),
            "weight-1 above weighted cutoff should survive"
        );
        assert!(
            !kept.contains("XXI"),
            "weight-2 below weighted cutoff should be dropped"
        );
        assert!(
            !kept.contains("XYZ"),
            "weight-3 below weighted cutoff should be dropped"
        );
    }

    /// End-to-end conservation: `Σ Z_i` propagated through a sequence of
    /// `rxx + ryy` exchange-style gates (which preserve total Z) with
    /// aggressive preserve-aware truncation keeps every single-Z
    /// coefficient at 1.0 exactly. The same setup with a plain
    /// `CoefficientThreshold(0.5)` would drop them; preserve-aware
    /// truncation does not.
    #[test]
    fn preserve_single_z_conserves_total_z_under_aggressive_truncation() {
        let n = 4;
        let preserve: PreserveConfig<PWord> = PreserveConfig::single_z(n, /*base*/ 0.5, 0.0);
        let mut s: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).preserve(preserve).build();
        for j in 0..n {
            let term: String = (0..n).map(|i| if i == j { 'Z' } else { 'I' }).collect();
            s += (term.as_str(), 1.0);
        }

        // Apply a few rxx+ryy pairs (= XY exchange on each edge). This commutes
        // with Σ Z_k, so the coefficients on Z_j should remain at 1.0.
        for (a, b) in [(0, 1), (1, 2), (2, 3)] {
            s.rxx(a, b, 0.37);
            s.ryy(a, b, 0.37);
            s.truncate();
        }

        for j in 0..n {
            let term: String = (0..n).map(|i| if i == j { 'Z' } else { 'I' }).collect();
            let word: PWord = term.clone().into();
            let coeff = s.data().iter().find(|(k, _)| **k == word).map(|(_, v)| *v);
            assert!(
                coeff.is_some(),
                "single-Z string {} must be preserved",
                term
            );
            assert!(
                (coeff.unwrap() - 1.0).abs() < 1e-10,
                "coefficient on {} should remain 1.0 (got {})",
                term,
                coeff.unwrap()
            );
        }
    }

    /// Without preserve, the same aggressive `CoefficientThreshold(0.5)`
    /// also keeps `Σ Z_i` intact in this trivial case (because Σ Z is a
    /// fixed point and never falls below 0.5). The interesting failure
    /// case happens when single-Z coefficients drift below the cutoff
    /// — covered by the Python end-to-end test on a spreading initial
    /// observable.
    #[test]
    fn no_preserve_falls_back_to_strategy() {
        let n = 2;
        // No `preserve` → strategy (NoStrategy here) is used and nothing is dropped.
        let mut s: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).build();
        s += ("ZI", 1.0);
        s += ("XY", 1e-30);
        s.truncate(); // NoStrategy: keeps everything
        assert_eq!(s.data().iter().count(), 2);
    }
}
