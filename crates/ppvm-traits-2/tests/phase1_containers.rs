// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Tests for the graded-algebra impls the crate ships **on the containers**
//! (`containers.rs`): the `HashMap<K, C, IdentityBuildHasher>` hash-join backend
//! (only the `Vec` coordinate-list backend has in-module unit tests), the
//! sesquilinear pairing, and the provided `Accumulate::accumulate` sugar. Plus a
//! stub `Columnar`/`KeyColumn` pinning the structure-of-arrays contract.
//!
//! The laws exercised are the machine-checked ones of
//! `lean/PPVM/Algebra/GradedMap.lean`: `accumulate_comm`/`accumulate_assoc`
//! (order-independence of the merge), `reduce_structural` (reduce drops exactly
//! the zero coefficients), `scale_scale`/`scale_accumulate`, `overlap_comm`
//! (the bilinear pairing is symmetric), and `hermitianOverlap_conj_symm` /
//! `hermitianOverlap_self_nonneg` (the sesquilinear one is conjugate-symmetric
//! with a nonnegative diagonal).

use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use num::Complex;

use ppvm_traits_2::batch::{Columnar, KeyBatch, KeyColumn, TermBatch, TermSink};
use ppvm_traits_2::graded::{Accumulate, Pair, Retain, Scale, Support};
use ppvm_traits_2::hash::{IdentityBuildHasher, Indexable};

type C = Complex<f64>;

/// A key whose structural digest is `tag` alone, while equality also compares
/// `body`. Two distinct keys can therefore share a digest on purpose — the case
/// the pass-through `IdentityBuildHasher` must still resolve by `Eq`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Key {
    tag: u64,
    body: u64,
}

impl Key {
    const fn new(tag: u64, body: u64) -> Self {
        Key { tag, body }
    }
}

impl Hash for Key {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Contract (design §"Indexable values"): `Hash` is exactly
        // `state.write_u64(self.key_hash())`.
        state.write_u64(self.key_hash());
    }
}

impl Indexable for Key {
    fn key_hash(&self) -> u64 {
        self.tag
    }
}

type Store<C> = HashMap<Key, C, IdentityBuildHasher>;

fn store<T: Clone>(terms: &[(Key, T)]) -> HashMap<Key, T, IdentityBuildHasher> {
    let mut m: HashMap<Key, T, IdentityBuildHasher> = HashMap::default();
    for (k, c) in terms {
        m.insert(*k, c.clone());
    }
    m
}

fn batch<T: Clone>(terms: &[(Key, T)]) -> TermBatch<Key, T> {
    let mut b = TermBatch::with_capacity(terms.len());
    for (k, c) in terms {
        b.push(*k, c.clone());
    }
    b
}

const A: Key = Key::new(1, 0);
const B: Key = Key::new(2, 0);
const D: Key = Key::new(3, 0);

// ---------------------------------------------------------------------------
// HashMap backend: L1 Accumulate.
// ---------------------------------------------------------------------------

#[test]
fn hashmap_accumulate_merges_onto_existing_keys() {
    let mut m: Store<f64> = HashMap::default();
    assert!(Support::is_empty(&m)); // provided default over `len`

    m.accumulate_batch(&batch(&[(A, 1.0), (B, 2.0), (A, 3.0)]));
    assert_eq!(Support::len(&m), 2);
    assert_eq!(Support::get(&m, &A), Some(4.0));
    assert_eq!(Support::get(&m, &B), Some(2.0));
    assert_eq!(Support::get(&m, &D), None);
    assert!(!Support::is_empty(&m));

    // `iter` exports every (key, coeff) pair (order is unspecified).
    let mut pairs: Vec<(Key, f64)> = Support::iter(&m).collect();
    pairs.sort_by_key(|(k, _)| k.tag);
    assert_eq!(pairs, vec![(A, 4.0), (B, 2.0)]);
}

#[test]
fn hashmap_accumulate_is_order_independent() {
    // `accumulate_comm` / `accumulate_assoc`: the merged support does not depend
    // on the order the terms arrive in, nor on the batch split.
    let mut one: Store<f64> = HashMap::default();
    one.accumulate_batch(&batch(&[(A, 1.0), (B, 2.0), (A, 3.0), (D, -1.0)]));

    let mut other: Store<f64> = HashMap::default();
    other.accumulate_batch(&batch(&[(D, -1.0), (A, 3.0)]));
    other.accumulate_batch(&batch(&[(A, 1.0)]));
    other.accumulate_batch(&batch(&[(B, 2.0)]));

    for k in [A, B, D] {
        assert_eq!(Support::get(&one, &k), Support::get(&other, &k));
    }
    assert_eq!(Support::len(&one), Support::len(&other));
}

#[test]
fn hashmap_accumulate_scalar_sugar_matches_a_batch_of_one() {
    // The design specifies the scalar `accumulate(k, c)` as *provided sugar over
    // a batch of one*; pin that it merges identically (and does not reduce).
    let mut sugar: Store<f64> = HashMap::default();
    sugar.accumulate(A, 1.0);
    sugar.accumulate(B, 2.0);
    sugar.accumulate(A, -1.0);

    let mut batched: Store<f64> = HashMap::default();
    batched.accumulate_batch(&batch(&[(A, 1.0), (B, 2.0), (A, -1.0)]));

    assert_eq!(Support::get(&sugar, &A), Support::get(&batched, &A));
    assert_eq!(Support::get(&sugar, &B), Support::get(&batched, &B));
    // Accumulation alone never reduces: the cancelled key is still supported
    // (with a zero coefficient) until `reduce` runs at finalize.
    assert_eq!(Support::get(&sugar, &A), Some(0.0));
    assert_eq!(Support::len(&sugar), 2);

    // Same on the Vec backend.
    let mut v: Vec<(Key, f64)> = Vec::new();
    v.accumulate(A, 1.0);
    v.accumulate(A, -1.0);
    assert_eq!(Support::len(&v), 1);
    assert_eq!(Support::get(&v, &A), Some(0.0));
}

#[test]
fn hashmap_reduce_drops_exactly_the_zero_coefficients() {
    // `reduce_structural`: reduce is a *finalize* step that drops zeros only.
    let mut m: Store<f64> = HashMap::default();
    m.accumulate_batch(&batch(&[(A, 1.0), (A, -1.0), (B, 2.0), (D, 0.0)]));
    assert_eq!(Support::len(&m), 3); // still there before reduce

    m.reduce();
    assert_eq!(Support::len(&m), 1);
    assert_eq!(Support::get(&m, &B), Some(2.0));
    assert_eq!(Support::get(&m, &A), None);
    assert_eq!(Support::get(&m, &D), None);

    // Idempotent.
    m.reduce();
    assert_eq!(Support::len(&m), 1);
}

// ---------------------------------------------------------------------------
// The pass-through hashing contract: distinct keys may share a digest.
// ---------------------------------------------------------------------------

#[test]
fn hashmap_keeps_distinct_keys_with_an_equal_digest_separate() {
    // `IdentityBuildHasher` hands hashbrown the key's own digest, so a digest
    // collision is resolved by `Eq` alone — two keys with the same `key_hash`
    // but different bodies must remain two terms, not merge into one.
    let k0 = Key::new(7, 0);
    let k1 = Key::new(7, 1);
    assert_eq!(k0.key_hash(), k1.key_hash());
    assert_ne!(k0, k1);

    let mut m: Store<f64> = HashMap::default();
    m.accumulate_batch(&batch(&[(k0, 1.0), (k1, 2.0), (k0, 10.0)]));
    assert_eq!(Support::len(&m), 2);
    assert_eq!(Support::get(&m, &k0), Some(11.0));
    assert_eq!(Support::get(&m, &k1), Some(2.0));
}

// ---------------------------------------------------------------------------
// L2 Scale and L3 Pair on the HashMap backend.
// ---------------------------------------------------------------------------

#[test]
fn hashmap_scale_and_overlap() {
    let mut a: Store<f64> = store(&[(A, 2.0), (B, 3.0)]);
    let b: Store<f64> = store(&[(A, 5.0), (D, 7.0)]);

    // `scale_scale`: scaling twice equals scaling by the product.
    a.scale(&2.0);
    a.scale(&3.0);
    assert_eq!(Support::get(&a, &A), Some(12.0));
    assert_eq!(Support::get(&a, &B), Some(18.0));

    // Only the shared key contributes: 12 * 5 = 60.
    assert_eq!(Pair::overlap(&a, &b), 60.0);
    // `overlap_comm`: the bilinear pairing is symmetric.
    assert_eq!(Pair::overlap(&b, &a), 60.0);
    // Disjoint supports pair to zero.
    let disjoint: Store<f64> = store(&[(D, 1.0)]);
    assert_eq!(Pair::overlap(&store(&[(A, 1.0)]), &disjoint), 0.0);
}

#[test]
fn hashmap_probe_batch_reports_hits_and_misses() {
    let m: Store<f64> = store(&[(A, 2.0), (B, 3.0)]);
    let mut keys: KeyBatch<Key> = KeyBatch::with_capacity(3);
    for k in [A, D, B] {
        keys.push(k);
    }
    let mut out = vec![None; keys.len()];
    Pair::probe_batch(&m, &keys, &mut out);
    assert_eq!(out, vec![Some(2.0), None, Some(3.0)]);
}

#[test]
fn hashmap_retain_filters_without_reducing() {
    let mut m: Store<f64> = store(&[(A, 2.0), (B, 0.5), (D, 0.0)]);
    Retain::retain(&mut m, |_, c| *c >= 1.0);
    assert_eq!(Support::len(&m), 1);
    assert_eq!(Support::get(&m, &A), Some(2.0));

    // The predicate sees the key too.
    let mut m2: Store<f64> = store(&[(A, 2.0), (B, 3.0)]);
    Retain::retain(&mut m2, |k, _| k.tag == B.tag);
    assert_eq!(Support::len(&m2), 1);
    assert_eq!(Support::get(&m2, &B), Some(3.0));
}

// ---------------------------------------------------------------------------
// The sesquilinear pairing, on both backends.
// ---------------------------------------------------------------------------

#[test]
fn hermitian_overlap_is_conjugate_symmetric_with_a_nonnegative_diagonal() {
    let a: Store<C> = store(&[(A, C::new(1.0, 2.0)), (B, C::new(-3.0, 1.0))]);
    let b: Store<C> = store(&[(A, C::new(0.0, 1.0)), (B, C::new(2.0, -1.0))]);

    let ab = Pair::hermitian_overlap(&a, &b);
    let ba = Pair::hermitian_overlap(&b, &a);
    // `hermitianOverlap_conj_symm`: ⟨a, b⟩ = conj(⟨b, a⟩).
    assert_eq!(ab, ba.conj());

    // Explicit value: conj(1+2i)·(i) + conj(-3+i)·(2-i) = ... computed directly.
    let expected = C::new(1.0, -2.0) * C::new(0.0, 1.0) + C::new(-3.0, -1.0) * C::new(2.0, -1.0);
    assert_eq!(ab, expected);

    // `hermitianOverlap_self_nonneg`: ⟨a, a⟩ is real and ≥ 0 (= Σ|a_k|²).
    let aa = Pair::hermitian_overlap(&a, &a);
    assert_eq!(aa.im, 0.0);
    assert!(aa.re >= 0.0);
    assert_eq!(aa.re, 1.0 + 4.0 + 9.0 + 1.0);

    // It is genuinely *not* the bilinear `overlap` over a complex ring.
    assert_ne!(ab, Pair::overlap(&a, &b));
}

#[test]
fn vec_hermitian_overlap_matches_the_hashmap_backend() {
    // The two backends are the same algebra; the sesquilinear pairing must
    // agree term for term.
    let terms_a = [(A, C::new(1.0, 2.0)), (B, C::new(-3.0, 1.0))];
    let terms_b = [(A, C::new(0.0, 1.0)), (D, C::new(2.0, -1.0))];

    let map_a: Store<C> = store(&terms_a);
    let map_b: Store<C> = store(&terms_b);
    let vec_a: Vec<(Key, C)> = terms_a.to_vec();
    let vec_b: Vec<(Key, C)> = terms_b.to_vec();

    assert_eq!(
        Pair::hermitian_overlap(&vec_a, &vec_b),
        Pair::hermitian_overlap(&map_a, &map_b)
    );
    assert_eq!(Pair::overlap(&vec_a, &vec_b), Pair::overlap(&map_a, &map_b));
}

// ---------------------------------------------------------------------------
// Columnar / KeyColumn: the structure-of-arrays contract.
// ---------------------------------------------------------------------------

/// A naive scalar column — the fallback layout the design explicitly permits —
/// used to pin the `KeyColumn` contract (`hash_into` reproduces `key_hash`,
/// `key_eq` confirms a join match, `gather` permutes).
#[derive(Default, Clone)]
struct KeyCol(Vec<Key>);

impl KeyColumn for KeyCol {
    type Key = Key;

    fn len(&self) -> usize {
        self.0.len()
    }
    fn capacity(&self) -> usize {
        self.0.capacity()
    }
    fn with_capacity(n: usize) -> Self {
        KeyCol(Vec::with_capacity(n))
    }
    fn push(&mut self, key: Key) {
        self.0.push(key);
    }
    fn hash_into(&self, out: &mut [u64]) {
        for (slot, k) in out.iter_mut().zip(self.0.iter()) {
            *slot = k.key_hash();
        }
    }
    fn key_eq(&self, i: usize, other: &Key) -> bool {
        self.0[i] == *other
    }
    fn gather(&self, indices: &[u32]) -> Self {
        KeyCol(indices.iter().map(|&i| self.0[i as usize]).collect())
    }
    fn get(&self, i: usize) -> Key {
        self.0[i]
    }
    // The in-place surface a column that is a *store* (the `ColumnStore`
    // backend) needs; on the naive scalar layout they are the `Vec` operations
    // they are named after.
    fn clear(&mut self) {
        self.0.clear();
    }
    fn set(&mut self, i: usize, key: Key) {
        self.0[i] = key;
    }
    fn truncate(&mut self, len: usize) {
        self.0.truncate(len);
    }
}

impl Columnar for Key {
    type Column = KeyCol;
}

#[test]
fn key_column_hash_into_reproduces_key_hash() {
    let mut col = <KeyCol as KeyColumn>::with_capacity(3);
    assert!(col.is_empty()); // provided default over `len`
    for k in [A, B, D] {
        col.push(k);
    }
    assert_eq!(col.len(), 3);
    assert!(!col.is_empty());

    let mut hashes = vec![0u64; col.len()];
    col.hash_into(&mut hashes);
    assert_eq!(hashes, vec![A.key_hash(), B.key_hash(), D.key_hash()]);

    // The same digests a `KeyBatch` precomputes — the two paths must agree.
    let mut kb: KeyBatch<Key> = KeyBatch::with_capacity(3);
    for k in [A, B, D] {
        kb.push(k);
    }
    kb.fill_hashes();
    assert_eq!(kb.hashes(), hashes.as_slice());
}

#[test]
fn key_column_key_eq_and_gather() {
    let mut col = KeyCol::default();
    for k in [A, B, D] {
        col.push(k);
    }

    // Join confirm without materializing the element.
    assert!(col.key_eq(1, &B));
    assert!(!col.key_eq(1, &A));
    // A digest match is not an equality match.
    let shadow = Key::new(B.tag, 99);
    assert_eq!(shadow.key_hash(), B.key_hash());
    assert!(!col.key_eq(1, &shadow));

    // `gather` permutes/selects plane-wise; `get` is the scalar fallback.
    let picked = col.gather(&[2, 0, 2]);
    assert_eq!(picked.len(), 3);
    assert_eq!(picked.get(0), D);
    assert_eq!(picked.get(1), A);
    assert_eq!(picked.get(2), D);
}
