// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Observable-aware truncation: a set of Pauli strings that the active
//! truncation [`Strategy`](crate::traits::Strategy) is *not allowed* to
//! drop.
//!
//! # Why
//!
//! For a transport diagnostic of the form `<Ô(t) Ô(0)>`, the answer
//! depends only on the projection of the propagated `Ô(t)` onto the
//! handful of Pauli strings that make up `Ô(0)`. Coefficient-magnitude
//! truncation (or any other truncation strategy) can drop those exact
//! strings when their coefficients drift toward the cutoff — typically
//! the tail of a spreading operator at far sites. The preserve mechanism
//! marks a small set of strings as *kept*, and they survive truncation
//! regardless of what the active strategy would otherwise do.
//!
//! Crucially this is *not* itself a strategy: it composes with **any**
//! [`Strategy`] (coefficient threshold, max weight, both, anything else)
//! via a post-filter inside [`PauliSum::truncate`](crate::sum::PauliSum::truncate)
//! — snapshot, run the strategy verbatim, re-insert any preserved keys
//! that the strategy dropped.
//!
//! # Usage
//!
//! ```ignore
//! use ppvm_runtime::prelude::*;
//! use ppvm_runtime::sum::preserve;
//! type Cfg = config::indexmap::ByteFxHashF64<1>;
//!
//! let mut s: PauliSum<Cfg> = PauliSum::builder()
//!     .n_qubits(4)
//!     // pair the keep-set with any strategy you like; here, default.
//!     .preserve_strings(preserve::single_z(4))
//!     .build();
//! s += ("ZIII", 1.0);
//! s.exchange(0, 1, 0.1);   // hypothetical: ZIII survives all auto-truncates
//! ```

use std::collections::HashSet;

use crate::traits::PauliWordTrait;

/// Build the keep-set of all single-`Z` Pauli strings on `n_qubits` qubits —
/// `Z_0, Z_1, …, Z_{n−1}`. The natural choice when the transport
/// diagnostic is `<Σ_j Z_j(t) Z_i(0)>` (z-magnetization spread).
pub fn single_z<W: PauliWordTrait>(n_qubits: usize) -> HashSet<W> {
    (0..n_qubits)
        .map(|i| {
            let s: String = (0..n_qubits)
                .map(|j| if j == i { 'Z' } else { 'I' })
                .collect();
            W::from(s)
        })
        .collect()
}

/// Build a keep-set from a user-specified list of Pauli strings (each
/// must be a length-`n_qubits` string over `{I, X, Y, Z}`).
pub fn from_strings<W: PauliWordTrait, I>(strings: I) -> HashSet<W>
where
    I: IntoIterator<Item = String>,
{
    strings.into_iter().map(W::from).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::word::PauliWord;

    type W = PauliWord<[u8; 1]>;

    #[test]
    fn single_z_builds_correct_set() {
        let s: HashSet<W> = single_z(3);
        assert_eq!(s.len(), 3);
        assert!(s.contains(&W::from("ZII")));
        assert!(s.contains(&W::from("IZI")));
        assert!(s.contains(&W::from("IIZ")));
        assert!(!s.contains(&W::from("ZZI")));
    }

    #[test]
    fn from_strings_round_trip() {
        let s: HashSet<W> = from_strings(["XYZ".to_string(), "ZZZ".to_string()]);
        assert_eq!(s.len(), 2);
        assert!(s.contains(&W::from("XYZ")));
        assert!(s.contains(&W::from("ZZZ")));
    }

    // ===========================================================================
    // End-to-end PauliSum integration tests
    // ===========================================================================
    //
    // These exercise the `PauliSum::truncate` snapshot-and-restore post-filter:
    // whatever the active strategy decides to drop, preserved strings come back.

    use crate::config::Config;
    use crate::config::fxhash::ByteF64;
    use crate::strategy::CoefficientThreshold;
    use crate::sum::PauliSum;
    use crate::traits::RotationTwo;

    type Cfg = ByteF64<1>;
    type CfgThr = ByteF64<1, CoefficientThreshold>;
    type PWord = <Cfg as Config>::PauliWordType;

    /// The active strategy (`CoefficientThreshold`) drops a tiny coefficient,
    /// but the preserved string is re-inserted.
    #[test]
    fn truncate_restores_preserved_string_dropped_by_strategy() {
        let mut s: PauliSum<CfgThr> = PauliSum::builder()
            .n_qubits(3)
            .strategy(CoefficientThreshold(0.5))
            .preserve_strings(single_z::<PWord>(3))
            .build();
        // Single-Z strings get tiny coefficients (below cutoff 0.5) — must survive.
        s += ("ZII", 1e-6);
        s += ("IZI", 1e-6);
        s += ("IIZ", 1e-6);
        // A non-preserved term well below cutoff — must be dropped.
        s += ("XYZ", 1e-6);
        // A non-preserved term above cutoff — must survive.
        s += ("XXX", 0.7);

        s.truncate();
        let kept: HashSet<String> = s.data().keys().map(|k| k.to_string()).collect();
        assert!(kept.contains("ZII"), "preserved ZII should be kept");
        assert!(kept.contains("IZI"), "preserved IZI should be kept");
        assert!(kept.contains("IIZ"), "preserved IIZ should be kept");
        assert!(
            !kept.contains("XYZ"),
            "below-cutoff non-preserved XYZ should be dropped"
        );
        assert!(kept.contains("XXX"), "above-cutoff XXX should be kept");
    }

    /// End-to-end conservation: `Σ Z_i` propagated through a sequence of
    /// `rxx + ryy` exchange-style gates (which preserve total Z) with
    /// aggressive coefficient truncation keeps every single-Z coefficient
    /// at 1.0 exactly. The same setup without the preserve set would drop
    /// them once their coefficients dipped below the threshold.
    #[test]
    fn preserve_single_z_conserves_total_z_under_aggressive_truncation() {
        let n = 4;
        let mut s: PauliSum<CfgThr> = PauliSum::builder()
            .n_qubits(n)
            .strategy(CoefficientThreshold(0.5))
            .preserve_strings(single_z::<PWord>(n))
            .build();
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

    /// No preserve set → behaviour is identical to the bare strategy.
    #[test]
    fn empty_preserve_falls_back_to_strategy_unchanged() {
        let n = 2;
        let mut s: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).build();
        s += ("ZI", 1.0);
        s += ("XY", 1e-30);
        s.truncate(); // default strategy keeps everything
        assert_eq!(s.data().iter().count(), 2);
    }
}
