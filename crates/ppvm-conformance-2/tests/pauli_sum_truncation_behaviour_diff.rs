// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Behaviour-parity guard for **when** truncation happens.
//!
//! The old crate's contract is *deferred* truncation: a gate must **never** drop
//! a term on its own — terms leave the support only through an explicit
//! `truncate()` call by the caller. That contract is pinned on the old side by
//! `ppvm-pauli-sum/tests/truncation_semantics.rs`, which documents a
//! since-reverted insertion-pruning optimization: dropping a produced term whose
//! own coefficient is sub-threshold is wrong, because two sub-threshold
//! contributions to the same key can merge into a surviving one
//! (`|a + b| ≥ τ` while `|a|, |b| < τ`) when several gates run between
//! truncations.
//!
//! The `-2` engine briefly broke that contract: `Sum::apply`,
//! `Sum::rekey_bijective` and `Sum::rotate_in_place` each ran
//! `policy.truncate()` internally, so a gate *did* drop terms on its own. Every
//! existing test still passed (they all truncate right after each gate, where
//! the two schedules agree), and the numeric golden masters matched — the
//! divergence only shows up when a caller defers truncation across gates. These
//! tests drive exactly that schedule on **both** crates and diff the result.
//!
//! Config parity (both sides): `[u8; 8]` storage, `f64` coefficients, and the
//! same single-rule truncation (old `CoefficientThreshold`/`MaxPauliWeight`
//! strategies, new the identically-named policies).

use ppvm_conformance_2::assert_close;

use ppvm_pauli_sum::config::fxhash::ByteF64 as OldByteF64;
use ppvm_pauli_sum::strategy::{
    CoefficientThreshold as OldCoeffThreshold, MaxPauliWeight as OldMaxWeight,
};
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_traits::traits::{Clifford as OldClifford, RotationOne as OldRotationOne};

use ppvm_pauli_sum_2::{
    CoefficientThreshold as NewCoeffThreshold, HashMapStore, MaxPauliWeight as NewMaxWeight,
    PauliWord as NewPauliWord, Sum,
};
use ppvm_traits_2::{Clifford as NewClifford, RotationOne as NewRotationOne};

/// Storage-matched key: `[u8; 8]`, exactly the old side's `ByteF64<8, _>`.
type NewKey = NewPauliWord<[u8; 8]>;

type OldThreshSum = OldPauliSum<OldByteF64<8, OldCoeffThreshold>>;
type NewThreshSum = Sum<HashMapStore<NewKey, f64>, NewCoeffThreshold>;

type OldWeightSum = OldPauliSum<OldByteF64<8, OldMaxWeight>>;
type NewWeightSum = Sum<HashMapStore<NewKey, f64>, NewMaxWeight>;

/// Sorted `(pauli_string, coeff)` view of an old sum's support. A macro, not a
/// generic fn: each test uses a different strategy parameter, and the old
/// `iter()` yields borrowed coefficients while the new yields owned ones.
macro_rules! old_support {
    ($sum:expr) => {{
        let mut v: Vec<(String, f64)> = $sum.iter().map(|(k, c)| (k.to_string(), *c)).collect();
        v.sort_by(|a, b| a.0.cmp(&b.0));
        v
    }};
}

/// Sorted `(pauli_string, coeff)` view of a new sum's support.
macro_rules! new_support {
    ($sum:expr) => {{
        let mut v: Vec<(String, f64)> = $sum.iter().map(|(k, c)| (k.to_string(), c)).collect();
        v.sort_by(|a, b| a.0.cmp(&b.0));
        v
    }};
}

/// Assert two supports agree as sorted `(string, coeff)` sets.
#[track_caller]
fn assert_support_eq(old: &[(String, f64)], new: &[(String, f64)], label: &str) {
    let old_keys: Vec<&str> = old.iter().map(|(k, _)| k.as_str()).collect();
    let new_keys: Vec<&str> = new.iter().map(|(k, _)| k.as_str()).collect();
    assert_eq!(
        old_keys, new_keys,
        "[{label}] support keys differ\nold={old:?}\nnew={new:?}"
    );
    for ((_, o), (_, n)) in old.iter().zip(new.iter()) {
        assert_close(*o, *n, 1e-12);
    }
}

// ---------------------------------------------------------------------------
// 1. `rotate_in_place` (rx) must not truncate — the exact scenario the old
//    crate's `truncation_semantics.rs` pins, run differentially.
// ---------------------------------------------------------------------------

/// Two `rx(θ)` on a single `Z` with **no truncate in between**. Each rotation
/// contributes a `Y` branch at `≈ sin(θ) = 0.030 < τ = 0.05`; the accumulated
/// `Y ≈ sin(2θ) ≈ 0.060` is *above* `τ` and must survive the final truncate.
///
/// An engine that truncates inside the gate kills the first (sub-threshold,
/// not-yet-merged) `Y` and never accumulates the second onto it — losing a term
/// the old crate keeps.
#[test]
fn deferred_truncation_across_two_rotations_matches_old() {
    const TAU: f64 = 0.05;
    const THETA: f64 = 0.03;

    let mut old: OldThreshSum = OldPauliSum::builder()
        .n_qubits(1)
        .strategy(OldCoeffThreshold(TAU))
        .build();
    old += ("Z", 1.0);
    old.rx(0, THETA);
    old.rx(0, THETA);
    old.truncate();

    let mut new: NewThreshSum = NewThreshSum::from_terms_with_policy(
        1,
        NewCoeffThreshold { threshold: TAU },
        [(NewKey::from("Z"), 1.0)],
    );
    new.rx(0, THETA);
    new.rx(0, THETA);
    new.truncate();

    let os = old_support!(old);
    let ns = new_support!(new);

    // Setup sanity: the merged `Y` really is above threshold, so this guards a
    // *dropped surviving term*, not just a coefficient nudge.
    assert!(
        os.iter().any(|(k, c)| k == "Y" && c.abs() >= TAU),
        "test setup broken: old should keep an above-threshold Y, got {os:?}"
    );
    assert_support_eq(&os, &ns, "two rx, truncate deferred");
}

/// The same rotations with `truncate()` after **each** gate. Here the two
/// schedules genuinely differ in result (the sub-threshold `Y` is legitimately
/// dropped mid-circuit), and old and new must agree on *that* answer too — this
/// pins that the fix removed the internal truncate without disturbing the
/// explicit one.
#[test]
fn eager_truncation_across_two_rotations_matches_old() {
    const TAU: f64 = 0.05;
    const THETA: f64 = 0.03;

    let mut old: OldThreshSum = OldPauliSum::builder()
        .n_qubits(1)
        .strategy(OldCoeffThreshold(TAU))
        .build();
    old += ("Z", 1.0);
    old.rx(0, THETA);
    old.truncate();
    old.rx(0, THETA);
    old.truncate();

    let mut new: NewThreshSum = NewThreshSum::from_terms_with_policy(
        1,
        NewCoeffThreshold { threshold: TAU },
        [(NewKey::from("Z"), 1.0)],
    );
    new.rx(0, THETA);
    new.truncate();
    new.rx(0, THETA);
    new.truncate();

    let os = old_support!(old);
    let ns = new_support!(new);

    // Setup sanity: this schedule really does lose the Y (otherwise the test
    // would not distinguish the two schedules at all).
    assert!(
        !os.iter().any(|(k, _)| k == "Y"),
        "test setup broken: eager truncation should have dropped the Y, got {os:?}"
    );
    assert_support_eq(&os, &ns, "two rx, truncate after each");
}

// ---------------------------------------------------------------------------
// 2. `rekey_bijective` (Clifford) must not truncate — a CNOT raises a key's
//    Pauli weight past `MaxPauliWeight`, but the term stays until `truncate()`.
// ---------------------------------------------------------------------------

/// Under `MaxPauliWeight(1)`, `cnot(0, 1)` maps `XI ↦ XX` (weight 2). The
/// over-weight term must still be in the support *after the gate*; only the
/// explicit `truncate()` removes it. A second `cnot(0, 1)` (the gate is an
/// involution) maps it back to the weight-1 `XI`, which the old crate keeps —
/// an engine that truncates inside the Clifford re-key has already destroyed it.
#[test]
fn clifford_rekey_does_not_truncate_over_weight_terms() {
    let mut old: OldWeightSum = OldPauliSum::builder()
        .n_qubits(2)
        .strategy(OldMaxWeight(1))
        .build();
    old += ("XI", 1.0);
    old.cnot(0, 1);

    let mut new: NewWeightSum =
        NewWeightSum::from_terms_with_policy(2, NewMaxWeight(1), [(NewKey::from("XI"), 1.0)]);
    new.cnot(0, 1);

    // (a) Immediately after the gate: the over-weight `XX` survives on both.
    let os = old_support!(old);
    let ns = new_support!(new);
    assert!(
        os.iter().any(|(k, _)| k == "XX"),
        "test setup broken: old should still hold the over-weight XX, got {os:?}"
    );
    assert_support_eq(&os, &ns, "after cnot, no truncate");

    // (b) The gate is an involution: a second cnot brings the term back under
    //     the weight cap, so a full round trip is the identity.
    old.cnot(0, 1);
    new.cnot(0, 1);
    let os = old_support!(old);
    let ns = new_support!(new);
    assert_eq!(
        os,
        vec![("XI".to_string(), 1.0)],
        "old lost the term across the cnot round trip"
    );
    assert_support_eq(&os, &ns, "cnot round trip");

    // (c) Explicit truncation still works: re-raise the weight, then truncate.
    old.cnot(0, 1);
    old.truncate();
    new.cnot(0, 1);
    new.truncate();
    let os = old_support!(old);
    let ns = new_support!(new);
    assert!(
        os.is_empty(),
        "explicit truncate should have dropped the over-weight XX, got {os:?}"
    );
    assert_support_eq(&os, &ns, "cnot then explicit truncate");
}
