// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Two user-facing contracts the rest of the suite exercises only statistically:
//!
//! 1. **`overlap` is the bilinear (non-conjugated) trace pairing over the shared
//!    support** — old's `PauliSum::overlap` (`ppvm-pauli-sum/src/sum/trace.rs`).
//!    Its six unit tests are ported here as *differential* assertions: empty
//!    operand → `0`; self-overlap of a single `3.0` term → `9.0`; orthogonal
//!    Paulis → `0`; a dot product over shared keys; partial support; symmetry.
//!    (Old's map-side trace folds over every entry rather than probing, so old's
//!    `overlap` is `O(|A|·|B|)` where the new `Pair::overlap` probes — numerically
//!    identical, and the only allowed difference is speed.)
//!
//! 2. **Iteration order is not a semantic** — the hash backends are unordered on
//!    both sides, so no result may depend on insertion order. Asserted by
//!    building the same term list under several permutations and comparing the
//!    support, the pairing and a propagated circuit.

use ppvm_conformance_2::{
    assert_close, assert_supports_match, build_new_sum, build_old_sum, new_support, old_support,
    random_terms, seeded_rng,
};

use ppvm_traits::traits::Clifford as OldClifford;
use ppvm_traits_2::Clifford as NewClifford;

const TOL: f64 = 1e-12;

/// Build the matched old/new pair from a `(string, coeff)` list.
fn pair(
    n: usize,
    terms: &[(&str, f64)],
) -> (ppvm_conformance_2::OldSum, ppvm_conformance_2::NewSum) {
    let owned: Vec<(String, f64)> = terms.iter().map(|(w, c)| ((*w).to_string(), *c)).collect();
    (build_old_sum(n, &owned), build_new_sum(n, &owned))
}

// ---------------------------------------------------------------------------
// 1. `overlap` — old's six `trace.rs` unit tests, as differential assertions.
// ---------------------------------------------------------------------------

#[test]
fn overlap_with_empty_is_zero_on_both() {
    let (old_a, new_a) = pair(2, &[("XI", 1.0), ("ZZ", 2.0)]);
    let (old_e, new_e) = pair(2, &[]);
    assert_eq!(old_a.overlap(&old_e), 0.0);
    assert_eq!(new_a.overlap(&new_e), 0.0);
    assert_eq!(old_e.overlap(&old_a), 0.0);
    assert_eq!(new_e.overlap(&new_a), 0.0);
    // The empty–empty case too: a pairing of nothing is zero, not a panic.
    assert_eq!(new_e.overlap(&new_e), 0.0);
}

#[test]
fn self_overlap_of_a_single_term_squares_the_coefficient_on_both() {
    let (old_a, new_a) = pair(2, &[("XY", 3.0)]);
    assert_eq!(old_a.overlap(&old_a), 9.0);
    assert_eq!(new_a.overlap(&new_a), 9.0);
}

#[test]
fn orthogonal_paulis_pair_to_zero_on_both() {
    let (old_a, new_a) = pair(2, &[("XI", 1.0)]);
    let (old_b, new_b) = pair(2, &[("ZI", 1.0)]);
    assert_eq!(old_a.overlap(&old_b), 0.0);
    assert_eq!(new_a.overlap(&new_b), 0.0);
}

#[test]
fn overlap_is_the_dot_product_over_shared_keys_on_both() {
    let (old_a, new_a) = pair(2, &[("XI", 1.0), ("ZZ", 2.0)]);
    let (old_b, new_b) = pair(2, &[("XI", 3.0), ("ZZ", 4.0)]);
    let want = 1.0 * 3.0 + 2.0 * 4.0;
    assert_eq!(old_a.overlap(&old_b), want);
    assert_eq!(new_a.overlap(&new_b), want);
}

#[test]
fn keys_present_on_one_side_only_contribute_zero_on_both() {
    let (old_a, new_a) = pair(2, &[("XI", 1.0), ("ZZ", 2.0)]);
    let (old_b, new_b) = pair(2, &[("ZZ", 5.0), ("YY", 7.0)]);
    let want = 2.0 * 5.0;
    assert_eq!(old_a.overlap(&old_b), want);
    assert_eq!(new_a.overlap(&new_b), want);
}

#[test]
fn overlap_is_symmetric_on_both() {
    let (old_a, new_a) = pair(2, &[("XI", 1.0), ("ZZ", 2.0), ("YY", -0.5)]);
    let (old_b, new_b) = pair(2, &[("ZZ", 5.0), ("YY", 7.0), ("IX", 2.0)]);
    assert_eq!(old_a.overlap(&old_b), old_b.overlap(&old_a));
    assert_eq!(new_a.overlap(&new_b), new_b.overlap(&new_a));
    assert_close(old_a.overlap(&old_b), new_a.overlap(&new_b), TOL);
}

/// A zero-coefficient term contributes nothing to the pairing on either side —
/// so `overlap` stays insensitive to the (contractually preserved) zero terms.
#[test]
fn zero_coefficient_terms_do_not_change_the_pairing_on_both() {
    let (old_a, new_a) = pair(2, &[("XI", 1.0), ("ZZ", 2.0)]);
    let (old_z, new_z) = pair(2, &[("XI", 1.0), ("ZZ", 2.0), ("YY", 0.0)]);
    let (old_b, new_b) = pair(2, &[("XI", 3.0), ("ZZ", 4.0), ("YY", 9.0)]);
    assert_eq!(old_a.overlap(&old_b), old_z.overlap(&old_b));
    assert_eq!(new_a.overlap(&new_b), new_z.overlap(&new_b));
}

// ---------------------------------------------------------------------------
// 2. Iteration/insertion order is not a semantic.
// ---------------------------------------------------------------------------

/// Compare two `(key, coeff)` supports for equal key sets and coefficients
/// within `tol`.
///
/// Not exact equality: a permuted insertion order makes duplicate keys
/// accumulate in a different order, and `f64` addition is not associative, so the
/// last ulp legitimately moves. What must not move is the key *set* or the value
/// to any meaningful precision.
#[track_caller]
fn assert_same_terms(a: &[(String, f64)], b: &[(String, f64)], tol: f64, label: &str) {
    assert_eq!(
        a.iter().map(|(k, _)| k).collect::<Vec<_>>(),
        b.iter().map(|(k, _)| k).collect::<Vec<_>>(),
        "[{label}] key set depends on insertion order"
    );
    for ((k, x), (_, y)) in a.iter().zip(b.iter()) {
        assert!(
            (x - y).abs() <= tol.max(x.abs() * 1e-12),
            "[{label}] coefficient at {k} depends on insertion order: {x} vs {y}"
        );
    }
}

/// Building the same terms in a different insertion order must give the same
/// support and the same pairing on both engines — neither backend promises an
/// order, so nothing may depend on one (beyond `f64` reassociation in the last
/// ulp, which is what the tolerance admits).
#[test]
fn insertion_order_does_not_change_the_result_on_either_engine() {
    for &seed in &[1u64, 42, 777] {
        let mut rng = seeded_rng(seed);
        for &n in &[3usize, 6] {
            let terms = random_terms(&mut rng, n, 24);
            let mut reversed = terms.clone();
            reversed.reverse();
            let mut rotated = terms.clone();
            rotated.rotate_left(7);

            let base_old = build_old_sum(n, &terms);
            let base_new = build_new_sum(n, &terms);

            for (label, perm) in [("reversed", &reversed), ("rotated", &rotated)] {
                let o = build_old_sum(n, perm);
                let nw = build_new_sum(n, perm);

                assert_same_terms(
                    &old_support(&base_old),
                    &old_support(&o),
                    TOL,
                    &format!("old/{label}"),
                );
                assert_same_terms(
                    &new_support(&base_new),
                    &new_support(&nw),
                    TOL,
                    &format!("new/{label}"),
                );
                assert_supports_match(&o, &nw, TOL);
                assert_close(
                    base_new.overlap(&base_new),
                    nw.overlap(&nw),
                    TOL.max(base_new.overlap(&base_new).abs() * 1e-12),
                );
            }
        }
    }
}

/// The same, after a short Clifford propagation: a re-key rebuilds the support,
/// so if any gate leaked an order dependence it would show here rather than at
/// construction.
#[test]
fn insertion_order_does_not_change_a_propagated_circuit() {
    let mut rng = seeded_rng(2024);
    let n = 6usize;
    let terms = random_terms(&mut rng, n, 20);
    let mut shuffled = terms.clone();
    shuffled.rotate_left(5);
    shuffled.reverse();

    let mut a_old = build_old_sum(n, &terms);
    let mut a_new = build_new_sum(n, &terms);
    let mut b_old = build_old_sum(n, &shuffled);
    let mut b_new = build_new_sum(n, &shuffled);

    for i in 0..n - 1 {
        a_old.h(i);
        b_old.h(i);
        a_new.h(i);
        b_new.h(i);
        a_old.cnot(i, i + 1);
        b_old.cnot(i, i + 1);
        a_new.cnot(i, i + 1);
        b_new.cnot(i, i + 1);
        a_old.s(i + 1);
        b_old.s(i + 1);
        a_new.s(i + 1);
        b_new.s(i + 1);
    }

    assert_same_terms(
        &old_support(&a_old),
        &old_support(&b_old),
        TOL,
        "old/circuit",
    );
    assert_same_terms(
        &new_support(&a_new),
        &new_support(&b_new),
        TOL,
        "new/circuit",
    );
    assert_supports_match(&a_old, &a_new, TOL);
    assert_supports_match(&b_old, &b_new, TOL);
}
