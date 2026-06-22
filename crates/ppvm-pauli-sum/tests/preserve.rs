// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the `preserve_strings` snapshot-and-restore
//! post-filter in [`PauliSum::truncate`]: whatever the active strategy
//! decides to drop, preserved strings come back.

use std::collections::HashSet;

use ppvm_pauli_sum::config::Config;
use ppvm_pauli_sum::config::fxhash::ByteF64;
use ppvm_pauli_sum::prelude::*;
use ppvm_pauli_sum::strategy::CoefficientThreshold;
use ppvm_pauli_sum::sum::PauliSum;

type Cfg = ByteF64<1>;
type CfgThr = ByteF64<1, CoefficientThreshold>;
type PWord = <Cfg as Config>::PauliWordType;

fn single_z(n_qubits: usize) -> HashSet<PWord> {
    (0..n_qubits)
        .map(|i| {
            let s: String = (0..n_qubits)
                .map(|j| if j == i { 'Z' } else { 'I' })
                .collect();
            PWord::from(s)
        })
        .collect()
}

/// The active strategy (`CoefficientThreshold`) drops a tiny coefficient,
/// but the preserved string is re-inserted.
#[test]
fn truncate_restores_preserved_string_dropped_by_strategy() {
    let mut s: PauliSum<CfgThr> = PauliSum::builder()
        .n_qubits(3)
        .strategy(CoefficientThreshold(0.5))
        .preserve_strings(single_z(3))
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
        .preserve_strings(single_z(n))
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
