// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Tests for both container backends. They share the fixtures
//! ([`Key`], [`batch`]) and one test spans both backends
//! ([`for_each_ref_agrees_with_iter_on_both_backends`]), so they stay in one
//! module rather than splitting alongside the impls.

use std::collections::HashMap;

use crate::batch::{TermBatch, TermSink};
use crate::graded::{Accumulate, Pair, Retain, Scale, Support};
use crate::hash::{IdentityBuildHasher, Indexable};

/// A minimal [`Indexable`] key, so the hash backend is testable here (the
/// only real one, `PauliWord`, lives downstream). `Hash` is exactly
/// `write_u64(key_hash())`, as the contract requires.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Key(u64);

impl std::hash::Hash for Key {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        state.write_u64(self.key_hash());
    }
}

impl Indexable for Key {
    fn key_hash(&self) -> u64 {
        self.0.wrapping_mul(0x9E37_79B9_7F4A_7C15)
    }
}

fn batch(terms: &[(&str, f64)]) -> TermBatch<String, f64> {
    let mut b = TermBatch::with_capacity(terms.len());
    for (k, c) in terms {
        b.push((*k).to_string(), *c);
    }
    b
}

#[test]
fn vec_accumulate_combines_keys() {
    let mut v: Vec<(String, f64)> = Vec::new();
    v.accumulate_batch(&batch(&[("a", 1.0), ("b", 2.0), ("a", 3.0)]));
    assert_eq!(Support::get(&v, &"a".to_string()), Some(4.0));
    assert_eq!(Support::get(&v, &"b".to_string()), Some(2.0));
    assert_eq!(Support::len(&v), 2);
}

#[test]
fn vec_reduce_drops_zero() {
    let mut v: Vec<(String, f64)> = Vec::new();
    v.accumulate_batch(&batch(&[("a", 1.0), ("a", -1.0), ("b", 2.0)]));
    v.reduce();
    assert_eq!(Support::len(&v), 1);
    assert_eq!(Support::get(&v, &"b".to_string()), Some(2.0));
}

#[test]
fn vec_scale_and_overlap() {
    let mut a: Vec<(String, f64)> = Vec::new();
    a.accumulate_batch(&batch(&[("x", 2.0), ("y", 3.0)]));
    let mut b: Vec<(String, f64)> = Vec::new();
    b.accumulate_batch(&batch(&[("x", 5.0), ("z", 7.0)]));
    a.scale(&2.0);
    // overlap = (2*2)*5 = 20; y and z do not match.
    assert_eq!(Pair::overlap(&a, &b), 20.0);
}

/// The borrowing scan must see exactly what [`Support::iter`] sees — same
/// multiset of `(key, coeff)` pairs — on both backends. It is an *override*
/// of a defaulted method, so nothing but a test keeps the two in step, and a
/// reader that silently visited a stale or partial view would corrupt every
/// trace.
#[test]
fn for_each_ref_agrees_with_iter_on_both_backends() {
    let terms = batch(&[("a", 1.0), ("b", 2.0), ("c", -3.0), ("a", 0.5)]);

    let mut v: Vec<(String, f64)> = Vec::new();
    v.accumulate_batch(&terms);
    let mut seen: Vec<(String, f64)> = Vec::new();
    v.for_each_ref(|k, c| seen.push((k.clone(), *c)));
    seen.sort_by(|a, b| a.0.cmp(&b.0));
    let mut want: Vec<(String, f64)> = Support::iter(&v).collect();
    want.sort_by(|a, b| a.0.cmp(&b.0));
    assert_eq!(seen, want);
    assert_eq!(seen.len(), 3);

    let mut m: HashMap<Key, f64, IdentityBuildHasher> = HashMap::default();
    for (k, c) in [(1u64, 1.0), (2, 2.0), (3, -3.0)] {
        m.insert(Key(k), c);
    }
    let mut seen: Vec<(Key, f64)> = Vec::new();
    m.for_each_ref(|k, c| seen.push((*k, *c)));
    seen.sort_by_key(|(k, _)| k.0);
    let mut want: Vec<(Key, f64)> = Support::iter(&m).collect();
    want.sort_by_key(|(k, _)| k.0);
    assert_eq!(seen, want);
    assert_eq!(seen.len(), 3);
}

#[test]
fn vec_retain_filters() {
    let mut v: Vec<(String, f64)> = Vec::new();
    v.accumulate_batch(&batch(&[("keep", 2.0), ("drop", 0.5)]));
    Retain::retain(&mut v, |_, c| *c >= 1.0);
    assert_eq!(Support::len(&v), 1);
    assert_eq!(Support::get(&v, &"keep".to_string()), Some(2.0));
}
