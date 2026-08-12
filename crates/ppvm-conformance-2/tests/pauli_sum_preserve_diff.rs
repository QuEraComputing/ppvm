// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Behaviour-parity guard for the **preserve set** — old's
//! `PauliSum::preserve_strings`, new's [`Sum::preserving`] keep-set (behavioural
//! contract 5 / architecture feature 11).
//!
//! Old (`ppvm-pauli-sum/src/sum/data.rs:271-306`) makes `truncate()` a
//! snapshot-and-restore post-filter: it records the preserved keys' *pre-truncate*
//! coefficients, runs the configured strategy **verbatim**, and re-inserts any
//! preserved key the strategy dropped — guarded by a membership test so a
//! survivor is never double-added. An empty set short-circuits to the bare
//! strategy with no snapshot scan.
//!
//! These are the old crate's `tests/preserve.rs` cases run **differentially**:
//! same seed, same strategy, same preserve set, both engines, identical support
//! required. The old-side conservation case drives `rxx`/`ryy` (a two-qubit
//! rotation the `-2` engine has not landed yet — `RotationTwo` is a later
//! component), so the propagation case here uses `rx`, which both sides have;
//! it exercises the same mechanism, a preserved key repeatedly pushed below the
//! floor and restored.
//!
//! Config parity (both sides): `[u8; 8]` storage, `f64` coefficients,
//! `CoefficientThreshold`/`MaxPauliWeight` with the same parameters.

use ppvm_conformance_2::assert_close;

use ppvm_pauli_sum::config::fxhash::ByteF64 as OldByteF64;
use ppvm_pauli_sum::strategy::{
    CoefficientThreshold as OldCoeffThreshold, MaxPauliWeight as OldMaxWeight,
};
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_pauli_word::word::PauliWord as OldWordT;
use ppvm_traits::traits::RotationOne as OldRotationOne;

use ppvm_pauli_sum_2::{
    CoefficientThreshold as NewCoeffThreshold, HashMapStore, MaxPauliWeight as NewMaxWeight,
    PauliWord as NewPauliWord, Sum,
};
use ppvm_traits_2::RotationOne as NewRotationOne;

/// Storage-matched key: `[u8; 8]`, exactly the old side's `ByteF64<8, _>`.
type NewKey = NewPauliWord<[u8; 8]>;
type OldKey = OldWordT<[u8; 8]>;

type OldThreshSum = OldPauliSum<OldByteF64<8, OldCoeffThreshold>>;
type NewThreshSum = Sum<HashMapStore<NewKey, f64>, NewCoeffThreshold>;

type OldWeightSum = OldPauliSum<OldByteF64<8, OldMaxWeight>>;
type NewWeightSum = Sum<HashMapStore<NewKey, f64>, NewMaxWeight>;

/// The `n` single-`Z` strings — old's `single_z` helper.
fn single_z_strings(n: usize) -> Vec<String> {
    (0..n)
        .map(|i| (0..n).map(|j| if j == i { 'Z' } else { 'I' }).collect())
        .collect()
}

macro_rules! old_support {
    ($sum:expr) => {{
        let mut v: Vec<(String, f64)> = $sum.iter().map(|(k, c)| (k.to_string(), *c)).collect();
        v.sort_by(|a, b| a.0.cmp(&b.0));
        v
    }};
}

macro_rules! new_support {
    ($sum:expr) => {{
        let mut v: Vec<(String, f64)> = $sum.iter().map(|(k, c)| (k.to_string(), c)).collect();
        v.sort_by(|a, b| a.0.cmp(&b.0));
        v
    }};
}

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
// 1. The strategy drops a tiny coefficient; the preserved string comes back.
//    (Old `tests/preserve.rs::truncate_restores_preserved_string_dropped_by_strategy`.)
// ---------------------------------------------------------------------------

#[test]
fn preserved_keys_survive_the_coefficient_floor_like_old() {
    const TAU: f64 = 0.5;
    let terms = [
        ("ZII", 1e-6), // preserved, below the floor
        ("IZI", 1e-6), // preserved, below the floor
        ("IIZ", 1e-6), // preserved, below the floor
        ("XYZ", 1e-6), // NOT preserved, below the floor → dropped
        ("XXX", 0.7),  // NOT preserved, above the floor → kept
    ];

    let mut old: OldThreshSum = OldPauliSum::builder()
        .n_qubits(3)
        .strategy(OldCoeffThreshold(TAU))
        .preserve_strings(single_z_strings(3).into_iter().map(OldKey::from).collect())
        .build();
    for (s, c) in terms {
        old += (s, c);
    }
    old.truncate();

    let mut new: NewThreshSum = NewThreshSum::with_policy(3, NewCoeffThreshold { threshold: TAU })
        .preserving(single_z_strings(3).into_iter().map(NewKey::from));
    for (s, c) in terms {
        new += (NewKey::from(s), c);
    }
    new.truncate();

    let (o, n) = (old_support!(old), new_support!(new));
    assert_support_eq(&o, &n, "coefficient floor + preserve set");
    let keys: Vec<&str> = n.iter().map(|(k, _)| k.as_str()).collect();
    assert_eq!(keys, vec!["IIZ", "IZI", "XXX", "ZII"]);
}

// ---------------------------------------------------------------------------
// 2. Empty set → the bare strategy, byte for byte (the hot-path short-circuit).
//    (Old `tests/preserve.rs::empty_preserve_falls_back_to_strategy_unchanged`.)
// ---------------------------------------------------------------------------

#[test]
fn empty_preserve_set_is_the_bare_strategy_like_old() {
    const TAU: f64 = 0.5;
    let terms = [("ZI", 1.0), ("XY", 1e-30)];

    let mut old: OldThreshSum = OldPauliSum::builder()
        .n_qubits(2)
        .strategy(OldCoeffThreshold(TAU))
        .build();
    for (s, c) in terms {
        old += (s, c);
    }
    assert!(old.preserve_strings().is_empty());
    old.truncate();

    let mut new: NewThreshSum = NewThreshSum::with_policy(2, NewCoeffThreshold { threshold: TAU });
    for (s, c) in terms {
        new += (NewKey::from(s), c);
    }
    assert!(new.preserved_keys().is_empty());
    new.truncate();

    assert_support_eq(&old_support!(old), &new_support!(new), "empty preserve set");
    assert_eq!(new.len(), 1);
}

// ---------------------------------------------------------------------------
// 3. Repeated propagate + truncate: a preserved key pushed below the floor over
//    and over is restored every time, at its PRE-truncate coefficient.
//    (Old `tests/preserve.rs::preserve_single_z_conserves_total_z_under_
//    aggressive_truncation`, with `rx` standing in for the not-yet-ported
//    `rxx`/`ryy`.)
// ---------------------------------------------------------------------------

#[test]
fn preserved_keys_survive_repeated_propagate_and_truncate_like_old() {
    const N: usize = 4;
    const TAU: f64 = 0.5;
    const THETA: f64 = 0.37;
    const ROUNDS: usize = 12;

    let mut old: OldThreshSum = OldPauliSum::builder()
        .n_qubits(N)
        .strategy(OldCoeffThreshold(TAU))
        .preserve_strings(single_z_strings(N).into_iter().map(OldKey::from).collect())
        .build();
    let mut new: NewThreshSum = NewThreshSum::with_policy(N, NewCoeffThreshold { threshold: TAU })
        .preserving(single_z_strings(N).into_iter().map(NewKey::from));
    for s in single_z_strings(N) {
        old += (s.as_str(), 1.0);
        new += (NewKey::from(s.as_str()), 1.0);
    }

    // 12 rounds: the `Y` branch is below τ from the first truncate, so the
    // single-`Z` coefficient decays as `cos(θ)^k` — it crosses τ = 0.5 at k = 10.
    for _ in 0..ROUNDS {
        for q in 0..N {
            old.rx(q, THETA);
            new.rx(q, THETA);
        }
        old.truncate();
        new.truncate();
    }

    let (o, n) = (old_support!(old), new_support!(new));
    assert_support_eq(&o, &n, "repeated rx + truncate under a preserve set");
    // Every preserved single-Z key is still there on both sides — the whole
    // point of the mechanism (`cos(θ)^12 ≈ 0.42` is below τ = 0.5).
    for s in single_z_strings(N) {
        assert!(
            n.iter().any(|(k, _)| *k == s),
            "preserved key {s} must survive; support = {n:?}"
        );
    }
    // …and without the preserve set the same circuit loses them, so the test is
    // not vacuous.
    let mut bare: NewThreshSum = NewThreshSum::with_policy(N, NewCoeffThreshold { threshold: TAU });
    for s in single_z_strings(N) {
        bare += (NewKey::from(s.as_str()), 1.0);
    }
    for _ in 0..ROUNDS {
        for q in 0..N {
            bare.rx(q, THETA);
        }
        bare.truncate();
    }
    assert!(
        single_z_strings(N)
            .iter()
            .any(|s| !bare.contains_key(&NewKey::from(s.as_str()))),
        "without the preserve set at least one single-Z key must be dropped"
    );
}

// ---------------------------------------------------------------------------
// 4. The mechanism composes with ANY policy — here the weight cap — and a
//    preserved key the policy KEEPS is not double-added (old's `contains_with`
//    guard).
// ---------------------------------------------------------------------------

#[test]
fn preserve_composes_with_the_weight_cap_and_never_double_adds() {
    const MAX_W: usize = 1;
    let terms = [
        ("ZII", 2.0), // preserved, weight 1 → the policy KEEPS it (no restore)
        ("ZZI", 3.0), // preserved, weight 2 → dropped, then restored at 3.0
        ("XYZ", 5.0), // not preserved, weight 3 → dropped for good
    ];

    let mut old: OldWeightSum = OldPauliSum::builder()
        .n_qubits(3)
        .strategy(OldMaxWeight(MAX_W))
        .preserve_strings(["ZII", "ZZI"].into_iter().map(OldKey::from).collect())
        .build();
    for (s, c) in terms {
        old += (s, c);
    }
    old.truncate();

    let mut new: NewWeightSum = NewWeightSum::with_policy(3, NewMaxWeight(MAX_W))
        .preserving(["ZII", "ZZI"].into_iter().map(NewKey::from));
    for (s, c) in terms {
        new += (NewKey::from(s), c);
    }
    new.truncate();

    let (o, n) = (old_support!(old), new_support!(new));
    assert_support_eq(&o, &n, "weight cap + preserve set");
    // The survivor kept its coefficient (a restore through the accumulating
    // `add_term` without the membership guard would have made it 4.0).
    assert_eq!(new.get(&NewKey::from("ZII")), Some(2.0));
    assert_eq!(new.get(&NewKey::from("ZZI")), Some(3.0));
    assert!(!new.contains_key(&NewKey::from("XYZ")));
}

// ---------------------------------------------------------------------------
// 5. A preserved key that was never in the support is NOT conjured into it —
//    the snapshot only records keys that were present, as old's scan does.
// ---------------------------------------------------------------------------

#[test]
fn absent_preserved_key_is_not_inserted_like_old() {
    const TAU: f64 = 0.5;

    let mut old: OldThreshSum = OldPauliSum::builder()
        .n_qubits(2)
        .strategy(OldCoeffThreshold(TAU))
        .preserve_strings(["ZI", "IZ"].into_iter().map(OldKey::from).collect())
        .build();
    old += ("ZI", 1e-9);
    old.truncate();

    let mut new: NewThreshSum = NewThreshSum::with_policy(2, NewCoeffThreshold { threshold: TAU })
        .preserving(["ZI", "IZ"].into_iter().map(NewKey::from));
    new += (NewKey::from("ZI"), 1e-9);
    new.truncate();

    assert_support_eq(
        &old_support!(old),
        &new_support!(new),
        "absent preserved key",
    );
    assert_eq!(new.len(), 1);
    assert!(!new.contains_key(&NewKey::from("IZ")));
}
