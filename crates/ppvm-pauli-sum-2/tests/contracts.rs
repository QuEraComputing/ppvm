// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The user-facing behavioural contracts the engine must reproduce from the old
//! crate: construction defaults and the capacity override, the `(key, value)`
//! `contains`, exact/approximate equality, the accumulating operator surface,
//! zero-coefficient survival, and the pattern-trace expectation path.

use approx::assert_relative_eq;
use ppvm_pauli_sum_2::{
    CoefficientThreshold, CombinedPolicy, MaxPauliWeight, NoPolicy, PauliPattern, PauliSum,
    PauliWord, Policy, SiteSet,
};
use ppvm_traits_2::{Clifford, PauliError, RotationOne, Trace};

fn pw(s: &str) -> PauliWord {
    PauliWord::from(s)
}

// ---------------------------------------------------------------------------
// Contract 3 — construction defaults.
// ---------------------------------------------------------------------------

#[test]
fn coefficient_threshold_default_is_1e_12() {
    // Old: `impl Default for CoefficientThreshold { fn default() -> Self { Self(1e-12) } }`.
    assert_eq!(CoefficientThreshold::default().threshold, 1e-12);
    // The disable sentinel / member defaults.
    assert_eq!(MaxPauliWeight::default().max_weight(), usize::MAX);
    let combined = CombinedPolicy::<CoefficientThreshold, MaxPauliWeight>::default();
    assert_eq!(combined.0.threshold, 1e-12);
    assert_eq!(combined.1.max_weight(), usize::MAX);
}

#[test]
fn default_policy_truncates_below_1e_12() {
    // Behavioural form of the same contract: the DEFAULT policy drops 1e-15.
    let mut s: PauliSum<f64, CoefficientThreshold> = PauliSum::<f64, CoefficientThreshold>::new(2);
    s += (pw("XI"), 1e-15);
    s += (pw("IZ"), 1.0);
    s.truncate();
    assert_eq!(s.len(), 1);
    assert!(!s.contains_key(&pw("XI")));
}

// ---------------------------------------------------------------------------
// Contract 4 / feature 6 — capacity hints and the explicit override.
// ---------------------------------------------------------------------------

#[test]
fn capacity_hints_match_old() {
    let n = 7usize;
    assert_eq!(
        Policy::<PauliWord, f64>::capacity(&CoefficientThreshold::default(), n),
        n * 10
    );
    assert_eq!(
        Policy::<PauliWord, f64>::capacity(&MaxPauliWeight::default(), n),
        n * 10
    );
    // CombinedStrategy::capacity == min(members).
    let combined = CombinedPolicy(CoefficientThreshold { threshold: 1e-6 }, MaxPauliWeight(3));
    assert_eq!(
        Policy::<PauliWord, f64>::capacity(&combined, n),
        (n * 10).min(n * 10)
    );
    // NoPolicy reproduces old's 4ⁿ/2 wherever that is allocatable.
    assert_eq!(Policy::<PauliWord, f64>::capacity(&NoPolicy, 4), 1 << 7);
}

#[test]
fn explicit_capacity_override_is_reported() {
    let s: PauliSum<f64, CoefficientThreshold> =
        PauliSum::<f64, CoefficientThreshold>::with_policy(12, CoefficientThreshold::default());
    assert_eq!(s.capacity(), 120);
    assert_eq!(s.n_sites(), 12);
    assert_eq!(s.len(), 0);
    assert!(s.is_empty());

    let s: PauliSum<f64, CoefficientThreshold> =
        PauliSum::<f64, CoefficientThreshold>::with_capacity(
            12,
            CoefficientThreshold::default(),
            144,
        );
    assert_eq!(s.capacity(), 144, "the builder override must win");
}

// ---------------------------------------------------------------------------
// Contract 14 — `contains` is a (key, value) match.
// ---------------------------------------------------------------------------

#[test]
fn contains_matches_the_value_too() {
    let sum: PauliSum = PauliSum::from_terms(2, [(pw("XZ"), 1.0)]);
    assert!(sum.contains(&pw("XZ"), &1.0));
    assert!(!sum.contains(&pw("XZ"), &2.0));
    // The key-only predicate carries a distinct name.
    assert!(sum.contains_key(&pw("XZ")));
}

// ---------------------------------------------------------------------------
// Contract 6 — the accumulating operator surface.
// ---------------------------------------------------------------------------

#[test]
fn add_assign_accumulates_and_scales() {
    let mut s: PauliSum = PauliSum::new(4);
    s += (pw("IIII"), 1.0);
    s += (pw("IIII"), 1.0);
    assert_eq!(s.get(&pw("IIII")), Some(2.0), "`+=` must accumulate");

    // Bare-word `+=` adds coefficient one; the string form is old's spelling.
    s += pw("IYII");
    assert_eq!(s.get(&pw("IYII")), Some(1.0));
    s += ("IIZI", 1.0);
    assert_eq!(s.get(&pw("IIZI")), Some(1.0));

    // `*= 0.0` keeps every key, at 0.0 — it can only mutate, never remove.
    let before: Vec<PauliWord> = s.iter().map(|(k, _)| k).collect();
    s *= 0.0;
    assert_eq!(s.len(), before.len());
    for k in before {
        assert_eq!(s.get(&k), Some(0.0));
    }
}

#[test]
fn collection_surface_preserves_old_extend_and_into_iter_semantics() {
    let mut sum: PauliSum = PauliSum::from_terms(2, [(pw("XI"), 1.0)]);
    sum.extend([(pw("XI"), 7.0), (pw("IZ"), 2.0)]);
    assert_eq!(sum.get(&pw("XI")), Some(7.0), "Extend replaces duplicates");

    let mut terms: Vec<_> = sum
        .into_iter()
        .map(|(key, coeff)| (key.to_string(), coeff))
        .collect();
    terms.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    assert_eq!(terms, vec![("IZ".into(), 2.0), ("XI".into(), 7.0)]);
}

#[test]
#[should_panic(expected = "term key width must match")]
fn width_mismatch_is_a_debug_panic() {
    // Old: `debug_assert_eq!(self.n_qubits(), key.n_qubits())` on every insert.
    let mut s: PauliSum = PauliSum::new(4);
    s += (pw("XX"), 1.0);
}

// ---------------------------------------------------------------------------
// Contract 7 — equality is exact, support-sensitive, width-sensitive.
// ---------------------------------------------------------------------------

#[test]
fn equality_counts_zero_coefficient_terms() {
    let a: PauliSum = PauliSum::from_terms(2, [(pw("XZ"), 1.0)]);
    let b: PauliSum = PauliSum::from_terms(2, [(pw("XZ"), 1.0), (pw("ZX"), 0.0)]);
    assert_ne!(a, b, "an extra 0.0 term is a different support");
    assert!(!approx::AbsDiffEq::abs_diff_eq(&a, &b, f64::EPSILON));

    let c: PauliSum = PauliSum::from_terms(2, [(pw("XZ"), 1.0)]);
    assert_eq!(a, c);
}

#[test]
fn equality_is_width_sensitive_and_ignores_transient_buffers() {
    let a: PauliSum = PauliSum::from_terms(2, [(pw("XI"), 1.0)]);
    let b: PauliSum = PauliSum::from_terms(3, [(pw("XII"), 1.0)]);
    assert_ne!(a, b);

    // A clone taken *after* a gate compares equal: the store's aux/scratch are
    // transient and are not part of the value.
    let mut gated: PauliSum = PauliSum::from_terms(2, [(pw("XI"), 1.0)]);
    gated.h(0);
    let copy = gated.clone();
    assert_eq!(gated, copy);
}

// ---------------------------------------------------------------------------
// Contract 2 — exact-zero coefficients are never dropped.
// ---------------------------------------------------------------------------

#[test]
fn zero_eigenvalue_channel_keeps_the_term() {
    // λ_X = 1 − 2(p_Y + p_Z) = 0 for [0.0, 0.25, 0.25]: the X term stays at 0.0.
    let mut s: PauliSum = PauliSum::from_terms(2, [(pw("XI"), 1.0), (pw("ZI"), 1.0)]);
    let len = s.len();
    s.pauli_error(0, [0.0, 0.25, 0.25]);
    assert_eq!(s.len(), len, "a zeroed term must stay in the support");
    assert_eq!(s.get(&pw("XI")), Some(0.0));
}

#[test]
fn identity_rotation_still_inserts_the_zero_branch() {
    // `rx(q, 0.0)` has sinθ = 0, so the Z-bearing term produces a 0.0 Y branch —
    // which old's `map_insert`/`add_assign` inserts, so the engine must too.
    let mut s: PauliSum = PauliSum::from_terms(1, [(pw("Z"), 1.0)]);
    s.rx(0, 0.0);
    assert_eq!(s.len(), 2, "R_0 must add old's zero-coefficient branch key");
    assert_eq!(s.get(&pw("Z")), Some(1.0));
    assert_eq!(s.get(&pw("Y")), Some(0.0));

    // …and `reduce` is what canonicalizes it away, on the caller's word.
    s.reduce();
    assert_eq!(s.len(), 1);
}

#[test]
fn rotation_branch_cancellation_survives() {
    // Two `rx(θ)` branches onto the same key that cancel exactly: the key stays.
    let mut s: PauliSum = PauliSum::from_terms(1, [(pw("Z"), 1.0), (pw("Y"), 0.0)]);
    s.rx(0, std::f64::consts::FRAC_PI_2);
    assert!(s.contains_key(&pw("Y")));
    assert!(s.contains_key(&pw("Z")));
}

// ---------------------------------------------------------------------------
// Contract 1 — gates never auto-truncate.
// ---------------------------------------------------------------------------

#[test]
fn sub_threshold_branches_merge_when_no_truncate_runs() {
    // `ppvm-pauli-sum/tests/truncation_semantics.rs`: two sub-threshold rx
    // branches on the same key merge into an above-threshold term.
    let tau = 0.05;
    let mut s: PauliSum<f64, CoefficientThreshold> =
        PauliSum::<f64, CoefficientThreshold>::with_policy(
            1,
            CoefficientThreshold { threshold: tau },
        );
    s += (pw("Z"), 1.0);
    s.rx(0, 0.03);
    s.rx(0, 0.03);
    s.truncate();
    let y = s.get(&pw("Y")).expect("the merged Y branch must survive");
    assert!(y.abs() >= tau, "merged branch {y} fell under the floor");
    assert_relative_eq!(y, 0.06_f64.sin(), epsilon = 1e-12);
}

// ---------------------------------------------------------------------------
// Workload 5 — expectation extraction through the `Z?*` pattern trace.
// ---------------------------------------------------------------------------

#[test]
fn ghz_backward_zero_state_trace_is_one() {
    // `ppvm-pauli-sum/tests/ghz.rs::test_ghz_backward`.
    let mut s: PauliSum = PauliSum::from_terms(2, [(pw("ZZ"), 1.0)]);
    s.cnot(0, 1);
    s.h(0);
    assert_eq!(s.trace(&PauliPattern::zero_state()), 1.0);
}

#[test]
fn ghz_forward_matches_the_explicit_operator() {
    // `ppvm-pauli-sum/tests/ghz.rs::test_ghz_forward` — an exact operator match,
    // which is what the ported `PartialEq` is for.
    let mut s: PauliSum = PauliSum::from_terms(
        2,
        [
            (pw("ZI"), 1.0),
            (pw("IZ"), 1.0),
            (pw("ZZ"), 1.0),
            (pw("II"), 1.0),
        ],
    );
    s.h(0);
    s.cnot(0, 1);

    let ghz: PauliSum = PauliSum::from_terms(
        2,
        [
            (pw("XX"), 1.0),
            (pw("YY"), -1.0),
            (pw("ZZ"), 1.0),
            (pw("II"), 1.0),
        ],
    );
    assert_eq!(s, ghz);
}

#[test]
fn rotation_circuit_zero_state_trace_matches_the_frozen_constant() {
    // `ppvm-pauli-sum/tests/cnot.rs::test_cnot`.
    let mut s: PauliSum = PauliSum::from_terms(2, [(pw("ZZ"), 1.0)]);
    for q in 0..2 {
        s.rz(q, 1.1);
        s.ry(q, 2.1);
        s.rz(q, 1.1);
    }
    s.cnot(0, 1);
    s.rx(0, 2.1);
    s.rx(1, 2.1);
    assert_relative_eq!(s.trace(&PauliPattern::zero_state()), 0.18803675917759355);
}

#[test]
fn depolarizing_zero_state_trace_is_the_product_of_factors() {
    // `ppvm-pauli-sum/tests/noise.rs::test_depolarizing_error`: three independent
    // `depolarize1(i, pᵢ)` on ZZZ contract to Πᵢ (1 − 4pᵢ/3). `depolarize1(q, p)`
    // is `pauli_error(q, [p/3; 3])`.
    let mut s: PauliSum = PauliSum::from_terms(3, [(pw("ZZZ"), 1.0)]);
    let ps = [0.1, 0.2, 0.3];
    for (q, p) in ps.iter().enumerate() {
        s.pauli_error(q, [p / 3.0, p / 3.0, p / 3.0]);
    }
    let expected: f64 = ps.iter().map(|p| 1.0 - 4.0 * p / 3.0).product();
    assert!((s.trace(&PauliPattern::zero_state()) - expected).abs() < 1e-10);
}

#[test]
fn pattern_site_sets_classify_words() {
    let zero = PauliPattern::zero_state();
    assert!(zero.matches(&pw("IZZI")));
    assert!(!zero.matches(&pw("IXZI")));

    // An anchored prefix with a wildcard tail.
    let anchored = PauliPattern::new([SiteSet::X, SiteSet::NON_IDENTITY], SiteSet::ANY);
    assert!(anchored.matches(&pw("XYZI")));
    assert!(
        !anchored.matches(&pw("XIZI")),
        "site 1 must be non-identity"
    );
    assert!(!anchored.matches(&pw("YYZI")), "site 0 must be X");
}

/// The trace folds over the *borrowing* scan
/// (`Support::for_each_ref`) rather than the cloning `iter`, so this pins that
/// the scan sees the live support after every buffer-swapping operation — a
/// re-key ping-pongs `primary`/`aux`, and a scan that read the wrong buffer
/// would silently trace a stale sum.
#[test]
fn for_each_ref_sees_the_live_support_after_a_rekey() {
    let mut s: PauliSum = PauliSum::from_terms(3, [(pw("ZZI"), 2.0), (pw("IZZ"), 3.0)]);
    s.cnot(0, 1);
    s.h(2);

    let mut scanned: Vec<(PauliWord, f64)> = Vec::new();
    s.for_each_ref(|k, c| scanned.push((k.clone(), *c)));
    let mut expected: Vec<(PauliWord, f64)> = s.iter().collect();

    let key = |p: &PauliWord| p.to_string();
    scanned.sort_by_key(|(k, _)| key(k));
    expected.sort_by_key(|(k, _)| key(k));
    assert_eq!(scanned, expected);
    assert_eq!(scanned.len(), s.len());

    // …and the trace built on it agrees with an independent filter over `iter`.
    let pattern = PauliPattern::zero_state();
    let want: f64 = s
        .iter()
        .filter(|(k, _)| pattern.matches(k))
        .map(|(_, c)| c)
        .sum();
    assert_eq!(s.trace(&pattern), want);
}
