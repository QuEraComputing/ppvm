// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential + Lean-oracle coverage for **sum + sum** addition
//! (`Sum::add_sum`, `a += &b`).
//!
//! This is the sibling of the `MulAssign<PauliSum>` correction in
//! `pauli_sum_multiply_diff.rs`, and the second (and last) place the `-2` engine
//! deliberately diverges from the old crate.
//!
//! Old ships `impl AddAssign<PauliSum<T>> for PauliSum<T>` and the `&PauliSum`
//! variant (`ppvm-pauli-sum/src/sum/ops.rs:117-140`) as `self.extend(rhs)`.
//! `Extend` for `HashMap`/`IndexMap` **inserts**, i.e. *replaces* the value on a
//! duplicate key, so old's `A += B` overwrites the coefficient of every shared
//! key with `B`'s instead of summing it. That contradicts old's own single-term
//! path `sum += (word, coeff)`, which routes through the accumulating
//! `ACMapAddAssign::add_assign`, and no test in the old crate covers sum-plus-sum
//! — consistent with a latent bug.
//!
//! Lean oracle: free-module addition of two finitely-supported maps is pointwise
//! coefficient **addition** (`lean/PPVM/Algebra/GradedMap.lean` — the `C[K]`
//! module laws, `accumulateTerms_add` for partitioned batches). So the new engine
//! accumulates, and the reference value below is built from old's *own*
//! trustworthy single-term accumulate path — a genuine differential assertion of
//! the CORRECT value, not a self-consistency check.
//!
//! Config parity (both sides): `[u8; 8]` storage, `f64` coefficients, no
//! truncation policy.

use ppvm_conformance_2::assert_close;

use ppvm_pauli_sum::config::fxhash::ByteF64 as OldByteF64;
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;

use ppvm_pauli_sum_2::{HashMapStore, NoPolicy, PauliWord as NewPauliWord, Sum};

type NewKey = NewPauliWord<[u8; 8]>;
type OldSum = OldPauliSum<OldByteF64<8>>;
type NewSum = Sum<HashMapStore<NewKey, f64>, NoPolicy>;

fn build_old(n: usize, terms: &[(&str, f64)]) -> OldSum {
    let mut s: OldSum = OldPauliSum::builder().n_qubits(n).build();
    for (k, c) in terms {
        s += (*k, *c);
    }
    s
}

fn build_new(n: usize, terms: &[(&str, f64)]) -> NewSum {
    let mut s = NewSum::new(n);
    for (k, c) in terms {
        s += (NewKey::from(*k), *c);
    }
    s
}

fn old_support(s: &OldSum) -> Vec<(String, f64)> {
    let mut v: Vec<(String, f64)> = s.iter().map(|(k, c)| (k.to_string(), *c)).collect();
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
}

fn new_support(s: &NewSum) -> Vec<(String, f64)> {
    let mut v: Vec<(String, f64)> = s.iter().map(|(k, c)| (k.to_string(), c)).collect();
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
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

const N: usize = 2;
const A_TERMS: [(&str, f64); 2] = [("ZZ", 1.0), ("XI", 4.0)];
const B_TERMS: [(&str, f64); 2] = [("ZZ", 2.0), ("IY", -1.5)];

// ---------------------------------------------------------------------------
// 1. The divergence itself, pinned on both sides.
// ---------------------------------------------------------------------------

/// **THE SECOND (AND LAST) ALLOWED BEHAVIOUR DIVERGENCE.**
///
/// On a **shared** key old overwrites and new accumulates. The reference is
/// built from old's own single-term `+=`, which is the accumulate path old
/// itself uses everywhere else.
#[test]
fn sum_addition_accumulates_where_old_overwrites() {
    // OLD, through its sum+sum operator: `extend` → right-biased overwrite.
    // NOTE the by-VALUE overload: old's `impl AddAssign<&'a PauliSum<T>>`
    // (ops.rs:128) carries the bound `T::Map: ACMapIter<'a, Item = (W, C)>`,
    // while the shipped map's `ACMapIter` yields borrowed pairs — so the
    // by-reference overload is unsatisfiable for every shipped `Config` and does
    // not compile, exactly like `MulAssign<PauliSum>` (see
    // `pauli_sum_multiply_diff.rs`). Only the by-value form is live on old.
    let mut old = build_old(N, &A_TERMS);
    let old_rhs = build_old(N, &B_TERMS);
    old += old_rhs;
    let old_rhs = build_old(N, &B_TERMS);
    let old_operator = old_support(&old);
    assert_close(
        old_operator
            .iter()
            .find(|(k, _)| k == "ZZ")
            .expect("ZZ present")
            .1,
        2.0,
        1e-12,
    );

    // OLD, through its own single-term accumulate path — the Lean-correct value.
    let mut reference = build_old(N, &A_TERMS);
    for (k, c) in old_rhs.iter().map(|(k, c)| (*k, *c)) {
        reference += (k, c);
    }
    let reference = old_support(&reference);
    assert_close(
        reference.iter().find(|(k, _)| k == "ZZ").expect("ZZ").1,
        3.0,
        1e-12,
    );

    // NEW: the operator *is* the accumulate path.
    let mut new = build_new(N, &A_TERMS);
    let new_rhs = build_new(N, &B_TERMS);
    new += &new_rhs;
    assert_support_eq(&reference, &new_support(&new), "A += B (free-module sum)");

    // …and it is genuinely different from what old's operator produces, so the
    // divergence is observable rather than a coincidentally-equal rewrite.
    assert!(
        old_operator != new_support(&new),
        "the accumulating sum must differ from old's `extend` overwrite"
    );
}

// ---------------------------------------------------------------------------
// 2. Where old is right, new agrees exactly.
// ---------------------------------------------------------------------------

/// On **disjoint** supports `extend` and accumulate coincide, so old's operator
/// and new's must agree key-for-key. This confines the divergence to shared keys.
#[test]
fn disjoint_supports_match_old_exactly() {
    let a = [("ZZ", 1.0), ("XI", 4.0)];
    let b = [("IY", -1.5), ("YY", 0.25)];

    let mut old = build_old(N, &a);
    old += build_old(N, &b);

    let mut new = build_new(N, &a);
    new += &build_new(N, &b);

    assert_support_eq(&old_support(&old), &new_support(&new), "disjoint A += B");
}

/// New ships both `a += b` and `a += &b`, and they agree. (Old ships both
/// spellings too, but only the by-value one is instantiable — see the note in
/// `sum_addition_accumulates_where_old_overwrites`.)
#[test]
fn by_value_and_by_reference_operators_agree() {
    let mut new_ref = build_new(N, &A_TERMS);
    new_ref += &build_new(N, &B_TERMS);
    let mut new_val = build_new(N, &A_TERMS);
    new_val += build_new(N, &B_TERMS);
    assert_eq!(new_support(&new_ref), new_support(&new_val));
}

// ---------------------------------------------------------------------------
// 3. Contract 1/2 carry over to the new operator: no truncation, no zero-drop.
// ---------------------------------------------------------------------------

#[test]
fn sum_addition_neither_truncates_nor_drops_zeros() {
    let mut new = build_new(1, &[("Z", 1.0)]);
    new += &build_new(1, &[("Z", -1.0), ("X", 1e-30)]);
    assert_eq!(new.len(), 2);
    assert_eq!(new.get(&NewKey::from("Z")), Some(0.0));
    assert_eq!(new.get(&NewKey::from("X")), Some(1e-30));

    // Old's single-term accumulate path agrees (it is `entry().and_modify()
    // .or_insert()`, which inserts the zero).
    let mut old = build_old(1, &[("Z", 1.0)]);
    old += ("Z", -1.0);
    old += ("X", 1e-30);
    assert_support_eq(&old_support(&old), &new_support(&new), "zeros survive");
}
