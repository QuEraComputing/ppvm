// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Integration workload 4 — **truncation boundary and policy-cell semantics**,
//! diffed old-vs-new over the same `(weight profile) × (policy cell)` grid the
//! old crate's `benches/truncation-weight.rs` times and its `tests/strategy.rs`
//! pins.
//!
//! `truncate()` runs after *every* operation in the headline Trotter workload, so
//! its keep-rule is as user-facing as any gate. The rules, from
//! `ppvm-pauli-sum/src/strategy.rs`:
//!
//! * `CoefficientThreshold(τ)` retains `!v.cutoff(τ)` where `f64::cutoff(τ)` is
//!   `|v| < τ` — so a term at **exactly** `|c| == τ` is **KEPT**.
//! * `MaxPauliWeight(w)` retains `weight() <= w` — so `weight == w` is KEPT and
//!   `weight == w + 1` is dropped.
//! * `MaxPauliWeight(usize::MAX)` is the **disable sentinel**: an early return
//!   that skips the retain pass entirely, hence an exact no-op — including on
//!   zero-coefficient terms, which it must not remove.
//! * `CombinedStrategy(S1, S2)` is **two sequential retain passes**, and
//!   `capacity() = min(S1.capacity, S2.capacity)`; a single strategy's
//!   `capacity(n) = n * 10`.
//!
//! The state is the bench's: `n = 128` qubits, 1000 terms, three weight profiles
//! (target weight 3 / 50 / 120, positions spread by stride, Paulis varying with
//! the term index so keys are distinct, coefficients `1/(k+1)`), plus explicitly
//! placed boundary terms at `|c| == τ`, `|c| == τ − 1ulp`, `weight == w` and
//! `weight == w + 1`, so an off-by-one `>` vs `>=` cannot pass.
//!
//! Both sides are storage-matched to `[u8; 16]` (128-qubit capacity) with `f64`
//! coefficients.

use std::collections::BTreeSet;

use ppvm_pauli_sum::config::fxhash::ByteF64 as OldByteF64;
use ppvm_pauli_sum::strategy::{
    CoefficientThreshold as OldCoeffThreshold, CombinedStrategy, MaxPauliWeight as OldMaxWeight,
};
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_traits::traits::Strategy as OldStrategy;

use ppvm_pauli_sum_2::{
    CoefficientThreshold as NewCoeffThreshold, CombinedPolicy, HashMapStore,
    MaxPauliWeight as NewMaxWeight, PauliWord as NewPauliWord, Policy, Sum,
};

/// Qubit count of the truncation workload (old's `benches/truncation-weight.rs`).
const N: usize = 128;
/// Term count.
const TERMS: usize = 1000;
/// The coefficient floor used by the threshold cells.
const TAU: f64 = 1e-12;

type NewKey = NewPauliWord<[u8; 16]>;

/// Build the bench's term list for one weight profile: `TERMS` distinct keys of
/// (approximately) `target_weight` non-identity sites, positions spread by a
/// stride, Pauli letters varying with the term index, coefficients `1/(k+1)`.
fn profile_terms(target_weight: usize) -> Vec<(String, f64)> {
    let stride = (N / target_weight).max(1);
    (0..TERMS)
        .map(|k| {
            let mut w = vec!['I'; N];
            for j in 0..target_weight {
                let pos = (j * stride + k) % N;
                w[pos] = ['X', 'Y', 'Z'][(k + j) % 3];
            }
            // Keys must be distinct: perturb one extra site with the term index
            // in a way that cannot collide with the strided block above.
            let extra = (k * 7 + 3) % N;
            if w[extra] == 'I' {
                w[extra] = ['X', 'Y', 'Z'][k % 3];
            }
            (w.into_iter().collect::<String>(), 1.0 / (k as f64 + 1.0))
        })
        .collect()
}

/// A word with exactly `weight` non-identity sites (packed at the front).
fn word_of_weight(weight: usize, letter: char) -> String {
    (0..N)
        .map(|i| if i < weight { letter } else { 'I' })
        .collect()
}

/// The boundary terms: exactly at and just under the coefficient floor, and
/// exactly at and just over each weight cutoff under test.
fn boundary_terms(cutoffs: &[usize]) -> Vec<(String, f64)> {
    let mut v = vec![
        // |c| == τ  → KEPT (old drops iff |c| < τ).
        (word_of_weight(1, 'X'), TAU),
        // |c| == τ − 1ulp → dropped.
        (word_of_weight(2, 'X'), f64::from_bits(TAU.to_bits() - 1)),
        // |c| == τ + 1ulp → kept.
        (word_of_weight(3, 'X'), f64::from_bits(TAU.to_bits() + 1)),
        // Negative at the boundary: the rule is on the magnitude.
        (word_of_weight(4, 'X'), -TAU),
    ];
    for &w in cutoffs {
        if w == 0 || w >= N {
            continue;
        }
        // weight == w → KEPT; weight == w + 1 → dropped.
        v.push((word_of_weight(w, 'Y'), 1.0));
        v.push((word_of_weight(w + 1, 'Y'), 1.0));
    }
    v
}

/// Build the OLD sum with strategy `$strat` from `$terms`, truncate, and return
/// the surviving key set.
macro_rules! old_survivors {
    ($cfg:ty, $strat:expr, $terms:expr) => {{
        let mut s: OldPauliSum<$cfg> = OldPauliSum::builder()
            .n_qubits(N)
            .strategy($strat)
            .capacity(TERMS * 2)
            .build();
        for (w, c) in $terms {
            s += (w.as_str(), *c);
        }
        s.truncate();
        s.data()
            .iter()
            .map(|(k, _)| k.to_string())
            .collect::<BTreeSet<String>>()
    }};
}

/// Build the NEW sum with policy `$policy` from `$terms`, truncate, and return
/// the surviving key set.
macro_rules! new_survivors {
    ($policy_ty:ty, $policy:expr, $terms:expr) => {{
        let mut s: Sum<HashMapStore<NewKey, f64>, $policy_ty> =
            Sum::with_capacity(N, $policy, TERMS * 2);
        for (w, c) in $terms {
            s += (NewKey::from(w.as_str()), *c);
        }
        s.truncate();
        s.iter()
            .map(|(k, _)| k.to_string())
            .collect::<BTreeSet<String>>()
    }};
}

#[track_caller]
fn assert_same_keys(old: &BTreeSet<String>, new: &BTreeSet<String>, label: &str) {
    let only_old: Vec<&String> = old.difference(new).collect();
    let only_new: Vec<&String> = new.difference(old).collect();
    assert!(
        only_old.is_empty() && only_new.is_empty(),
        "[{label}] surviving key sets differ: {} only in old (e.g. {:?}), \
         {} only in new (e.g. {:?})",
        only_old.len(),
        only_old.first(),
        only_new.len(),
        only_new.first(),
    );
}

/// The `MaxPauliWeight(w)` cells, over every weight profile and every cutoff
/// (including the `usize::MAX` disable sentinel).
#[test]
fn max_weight_cells_keep_the_same_terms_as_old() {
    const CUTOFFS: [usize; 4] = [10, 100, 1000, usize::MAX];
    for target in [3usize, 50, 120] {
        let mut terms = profile_terms(target);
        terms.extend(boundary_terms(&CUTOFFS));
        for w in CUTOFFS {
            let old = old_survivors!(OldByteF64<16, OldMaxWeight>, OldMaxWeight(w), &terms);
            let new = new_survivors!(NewMaxWeight, NewMaxWeight(w), &terms);
            assert_same_keys(&old, &new, &format!("weight profile={target} cutoff={w}"));

            // The keep rule itself, independent of old: `weight <= w`.
            if w != usize::MAX {
                for k in &new {
                    let weight = k.chars().filter(|&c| c != 'I').count();
                    assert!(weight <= w, "kept a weight-{weight} term under cutoff {w}");
                }
            }
        }
    }
}

/// `MaxPauliWeight(usize::MAX)` is an exact no-op: every term survives, including
/// a zero-coefficient one (the sentinel skips the retain pass entirely, so it
/// cannot drop anything).
#[test]
fn max_weight_sentinel_is_an_exact_no_op_on_both() {
    let mut terms = profile_terms(50);
    // A zero-coefficient term and a full-width one, both of which must survive.
    terms.push((word_of_weight(N, 'Z'), 0.0));
    terms.push((word_of_weight(N, 'X'), 1.0));

    let old = old_survivors!(
        OldByteF64<16, OldMaxWeight>,
        OldMaxWeight(usize::MAX),
        &terms
    );
    let new = new_survivors!(NewMaxWeight, NewMaxWeight(usize::MAX), &terms);
    assert_same_keys(&old, &new, "sentinel");

    let all: BTreeSet<String> = terms.iter().map(|(w, _)| w.clone()).collect();
    assert_eq!(new, all, "the sentinel dropped a term");
    assert!(
        new.contains(&word_of_weight(N, 'Z')),
        "the sentinel dropped the zero-coefficient term"
    );
}

/// The `CoefficientThreshold(τ)` cell: the boundary rule is `|c| >= τ` keeps.
#[test]
fn coefficient_threshold_cell_keeps_the_same_terms_as_old() {
    for target in [3usize, 50, 120] {
        let mut terms = profile_terms(target);
        terms.extend(boundary_terms(&[]));

        let old = old_survivors!(
            OldByteF64<16, OldCoeffThreshold>,
            OldCoeffThreshold(TAU),
            &terms
        );
        let new = new_survivors!(
            NewCoeffThreshold,
            NewCoeffThreshold { threshold: TAU },
            &terms
        );
        assert_same_keys(&old, &new, &format!("threshold profile={target}"));

        // The two explicitly placed boundary terms, checked by name.
        assert!(
            new.contains(&word_of_weight(1, 'X')),
            "|c| == τ must be KEPT (old drops iff |c| < τ)"
        );
        assert!(
            !new.contains(&word_of_weight(2, 'X')),
            "|c| == τ − 1ulp must be dropped"
        );
        assert!(
            new.contains(&word_of_weight(4, 'X')),
            "|c| == −τ must be KEPT (the rule is on the magnitude)"
        );
    }
}

/// The `Combined(threshold, weight)` cells — two sequential retain passes.
#[test]
fn combined_policy_cells_keep_the_same_terms_as_old() {
    const CUTOFFS: [usize; 3] = [10, 100, usize::MAX];
    for target in [3usize, 50, 120] {
        let mut terms = profile_terms(target);
        terms.extend(boundary_terms(&CUTOFFS));
        for w in CUTOFFS {
            let old = old_survivors!(
                OldByteF64<16, CombinedStrategy<OldCoeffThreshold, OldMaxWeight>>,
                CombinedStrategy(OldCoeffThreshold(TAU), OldMaxWeight(w)),
                &terms
            );
            let new = new_survivors!(
                CombinedPolicy<NewCoeffThreshold, NewMaxWeight>,
                CombinedPolicy(NewCoeffThreshold { threshold: TAU }, NewMaxWeight(w)),
                &terms
            );
            assert_same_keys(&old, &new, &format!("combined profile={target} cutoff={w}"));
        }
    }
}

/// The capacity hints are observable (they size the maps and are reported by
/// `capacity()`), so they are part of the contract: `n * 10` for a single
/// policy, `min(..)` for the combination.
#[test]
fn capacity_hints_match_old_for_every_cell() {
    for n in [0usize, 1, 4, 12, 128] {
        let old_thresh = OldStrategy::capacity(&OldCoeffThreshold(TAU), n);
        let old_weight = OldStrategy::capacity(&OldMaxWeight(10), n);
        let old_combined = OldStrategy::capacity(
            &CombinedStrategy(OldCoeffThreshold(TAU), OldMaxWeight(10)),
            n,
        );

        let new_thresh = Policy::<NewKey, f64>::capacity(&NewCoeffThreshold { threshold: TAU }, n);
        let new_weight = Policy::<NewKey, f64>::capacity(&NewMaxWeight(10), n);
        let new_combined = Policy::<NewKey, f64>::capacity(
            &CombinedPolicy(NewCoeffThreshold { threshold: TAU }, NewMaxWeight(10)),
            n,
        );

        assert_eq!(new_thresh, old_thresh, "threshold capacity at n={n}");
        assert_eq!(new_weight, old_weight, "weight capacity at n={n}");
        assert_eq!(new_combined, old_combined, "combined capacity at n={n}");
        assert_eq!(new_thresh, n * 10, "single-policy capacity is n*10");
        assert_eq!(
            new_combined,
            new_thresh.min(new_weight),
            "combined capacity is the min"
        );
    }
}
