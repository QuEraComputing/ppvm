// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential tests for the one property the rest of the suite deliberately
//! filters away: **exact-zero coefficients**.
//!
//! Behavioural contract 2 — old has no `reduce` and no drop-zero logic anywhere,
//! so a term driven to exactly `0.0` stays in the support, and old's exact-map
//! `PartialEq` counts it (`ppvm-pauli-sum/tests/loss.rs::test_reset_channel`
//! asserts `state == (state2 *= 0.0)`). Every other differential test here
//! compares above-a-floor supports, precisely so that rotation-merge residue does
//! not make them flaky — which is exactly why a zero-handling divergence can hide
//! from them. These tests use scenarios whose zeros are *exact* (a `sinθ = 0`
//! branch, a `λ_P = 0` channel eigenvalue, a `× 0.0` scale, an inserted `0.0`),
//! so the whole key set — zeros included — can be compared with no floor at all.

use ppvm_conformance_2::{
    NewKey, NewSum, OldSum, build_new_sum, build_old_sum, new_support, old_support, random_terms,
    seeded_rng,
};

use ppvm_traits::traits::{PauliError as OldPauliError, RotationOne as OldRotationOne};
use ppvm_traits_2::{PauliError as NewPauliError, RotationOne as NewRotationOne};

const SEEDS: [u64; 4] = [1, 42, 777, 31337];
const WIDTHS: [usize; 4] = [1, 3, 5, 8];

/// Compare the FULL supports — no floor, no reduce — key set and coefficients
/// exactly. Zero-coefficient terms must be present on both sides or neither.
#[track_caller]
fn assert_supports_identical_including_zeros(old: &OldSum, new: &NewSum) {
    let os = old_support(old);
    let ns = new_support(new);
    assert_eq!(
        os.len(),
        ns.len(),
        "support size differs (zeros included): old {} vs new {}\nold={os:?}\nnew={ns:?}",
        os.len(),
        ns.len()
    );
    for (o, n) in os.iter().zip(ns.iter()) {
        assert_eq!(o.0, n.0, "key differs: old {} vs new {}", o.0, n.0);
        assert_eq!(
            o.1, n.1,
            "coefficient differs at {}: old {} vs new {}",
            o.0, o.1, n.1
        );
    }
}

// ---------------------------------------------------------------------------
// 1. The identity rotation `R_0` adds a zero-coefficient branch key on BOTH.
// ---------------------------------------------------------------------------

#[test]
fn identity_rotation_adds_the_same_zero_branch() {
    // `sinθ = 0` at θ = 0, so every anticommuting term produces a branch whose
    // coefficient is exactly 0.0. Old merges it (`map_insert` → `add_assign` →
    // `entry().or_insert(0.0)`); the new fast path must too.
    let mut s: OldSum = build_old_sum(1, &[("Z".to_string(), 1.0)]);
    let mut n: NewSum = build_new_sum(1, &[("Z".to_string(), 1.0)]);
    s.rx(0, 0.0);
    n.rx(0, 0.0);
    assert_eq!(s.len(), 2, "old is expected to keep the 0.0 Y branch");
    assert_supports_identical_including_zeros(&s, &n);

    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &w in &WIDTHS {
            let terms = random_terms(&mut rng, w, 12);
            for q in 0..w {
                for axis in 0..3usize {
                    let mut old = build_old_sum(w, &terms);
                    let mut new = build_new_sum(w, &terms);
                    match axis {
                        0 => {
                            old.rx(q, 0.0);
                            new.rx(q, 0.0);
                        }
                        1 => {
                            old.ry(q, 0.0);
                            new.ry(q, 0.0);
                        }
                        _ => {
                            old.rz(q, 0.0);
                            new.rz(q, 0.0);
                        }
                    }
                    assert_supports_identical_including_zeros(&old, &new);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 2. A zero channel eigenvalue keeps the term, at 0.0, on BOTH.
// ---------------------------------------------------------------------------

#[test]
fn zero_eigenvalue_channel_keeps_the_term() {
    // λ_X = 1 − 2(p_Y + p_Z) = 0 for [0.0, 0.25, 0.25].
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &w in &WIDTHS {
            let terms = random_terms(&mut rng, w, 12);
            for q in 0..w {
                let mut old = build_old_sum(w, &terms);
                let mut new = build_new_sum(w, &terms);
                let before = old.len();
                old.pauli_error(q, [0.0, 0.25, 0.25]);
                new.pauli_error(q, [0.0, 0.25, 0.25]);
                assert_eq!(old.len(), before, "old never removes a zeroed term");
                assert_supports_identical_including_zeros(&old, &new);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 3. `*= 0.0` keeps the whole key set on BOTH (the `test_reset_channel` shape).
// ---------------------------------------------------------------------------

#[test]
fn scale_by_zero_keeps_the_whole_key_set() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &w in &WIDTHS {
            let terms = random_terms(&mut rng, w, 15);
            let mut old = build_old_sum(w, &terms);
            let mut new = build_new_sum(w, &terms);
            old *= 0.0;
            new *= 0.0;
            assert_supports_identical_including_zeros(&old, &new);
        }
    }
}

// ---------------------------------------------------------------------------
// 4. Inserting a 0.0 coefficient, and cancelling to 0.0, keeps the key on BOTH.
// ---------------------------------------------------------------------------

#[test]
fn inserted_and_cancelled_zeros_survive_on_both() {
    for &w in &WIDTHS {
        let key = "Z".repeat(w);
        let terms = vec![(key.clone(), 0.0)];
        let old = build_old_sum(w, &terms);
        let new = build_new_sum(w, &terms);
        assert_eq!(old.len(), 1, "old keeps an inserted 0.0");
        assert_supports_identical_including_zeros(&old, &new);

        let cancelling = vec![(key.clone(), 1.5), (key.clone(), -1.5)];
        let old = build_old_sum(w, &cancelling);
        let new = build_new_sum(w, &cancelling);
        assert_eq!(old.len(), 1);
        assert!(new.contains(&NewKey::from(key.as_str()), &0.0));
        assert_supports_identical_including_zeros(&old, &new);
    }
}

// ---------------------------------------------------------------------------
// 5. `contains` is a (key, value) match on both crates.
// ---------------------------------------------------------------------------

#[test]
fn contains_predicate_matches_old() {
    let terms = vec![("XZ".to_string(), 1.0)];
    let old = build_old_sum(2, &terms);
    let new = build_new_sum(2, &terms);
    let old_key: ppvm_pauli_word::word::PauliWord<[u8; 8]> = "XZ".into();
    assert!(old.contains(&old_key, &1.0));
    assert!(!old.contains(&old_key, &2.0));
    assert!(new.contains(&NewKey::from("XZ"), &1.0));
    assert!(!new.contains(&NewKey::from("XZ"), &2.0));
    let _ = new_support(&new);
}
