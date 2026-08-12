// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Unit tests pinning the concrete batch leaf types of `ppvm-traits-2`
//! (`KeyBatch` / `TermBatch` / `TermSink` / `TermProducer`).
//!
//! These are the only executable *containers* the crate ships (the graded and
//! gate traits are definitions, impl'd downstream). They back the hash-join
//! contract of the design's §"Batch execution and the hash-join contract", so
//! their column invariants are pinned here:
//!   * `TermBatch` is a structure-of-arrays: the key and coeff columns stay
//!     parallel under `push`/`clear`, and `iter` re-synthesizes `(k, c)` pairs;
//!   * `KeyBatch::fill_hashes` reproduces each key's `Indexable::key_hash`
//!     verbatim into the parallel hash column;
//!   * a `TermProducer` can fan a single input term into a `TermBatch` sink.

use std::hash::{Hash, Hasher};

use ppvm_traits_2::batch::{KeyBatch, TermBatch, TermProducer, TermSink};
use ppvm_traits_2::hash::Indexable;

/// A trivial `Indexable` key: its digest is the wrapped value itself, so
/// `fill_hashes` has a value to reproduce bit for bit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Key(u64);

impl Hash for Key {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Contract: `Hash` writes exactly `key_hash()` as a single u64.
        state.write_u64(self.key_hash());
    }
}

impl Indexable for Key {
    fn key_hash(&self) -> u64 {
        self.0
    }
}

// ---------------------------------------------------------------------------
// TermBatch: structure-of-arrays, parallel columns, SoA iteration.
// ---------------------------------------------------------------------------

#[test]
fn term_batch_keeps_parallel_columns() {
    let mut b: TermBatch<Key, f64> = TermBatch::with_capacity(3);
    assert!(b.is_empty());
    assert_eq!(b.len(), 0);

    // `push` is the `TermSink` append side.
    TermSink::push(&mut b, Key(10), 1.5);
    TermSink::push(&mut b, Key(20), -2.0);
    TermSink::push(&mut b, Key(30), 0.25);

    assert_eq!(b.len(), 3);
    assert!(!b.is_empty());

    // Columns stay parallel and in insertion order.
    assert_eq!(b.keys().keys(), &[Key(10), Key(20), Key(30)]);
    assert_eq!(b.coeffs(), &[1.5, -2.0, 0.25]);

    // `iter` re-synthesizes the (key, coeff) pairs from the two columns.
    let pairs: Vec<(Key, f64)> = b.iter().map(|(k, c)| (*k, *c)).collect();
    assert_eq!(
        pairs,
        vec![(Key(10), 1.5), (Key(20), -2.0), (Key(30), 0.25)]
    );

    // `clear` empties both columns.
    b.clear();
    assert!(b.is_empty());
    assert_eq!(b.keys().keys().len(), 0);
    assert_eq!(b.coeffs().len(), 0);
}

// ---------------------------------------------------------------------------
// KeyBatch::fill_hashes reproduces Indexable::key_hash verbatim.
// ---------------------------------------------------------------------------

#[test]
fn key_batch_fill_hashes_matches_key_hash() {
    let mut kb: KeyBatch<Key> = KeyBatch::with_capacity(4);
    // Hash column is empty until `fill_hashes` runs.
    assert_eq!(kb.hashes().len(), 0);

    let keys = [Key(0), Key(1), Key(u64::MAX), Key(0xdead_beef)];
    for k in keys {
        kb.push(k);
    }
    assert_eq!(kb.len(), keys.len());

    kb.fill_hashes();
    let expected: Vec<u64> = keys.iter().map(Indexable::key_hash).collect();
    assert_eq!(kb.hashes(), expected.as_slice());

    // Idempotent: re-filling does not double the column.
    kb.fill_hashes();
    assert_eq!(kb.hashes(), expected.as_slice());
}

// ---------------------------------------------------------------------------
// TermProducer fans one input term into a TermBatch sink.
// ---------------------------------------------------------------------------

/// A stub producer that emits the input term unchanged plus a scaled copy on a
/// derived key — a stand-in for the rotation branch's two-term fan-out.
struct SplitProducer;

impl TermProducer<Key, f64> for SplitProducer {
    fn produce<S: TermSink<Key, f64>>(&self, key: &Key, coeff: &f64, sink: &mut S) {
        sink.push(*key, *coeff);
        sink.push(Key(key.0 + 1), coeff * 0.5);
    }
}

#[test]
fn term_producer_feeds_sink() {
    let mut sink: TermBatch<Key, f64> = TermBatch::new();
    SplitProducer.produce(&Key(7), &4.0, &mut sink);

    assert_eq!(sink.len(), 2);
    let pairs: Vec<(Key, f64)> = sink.iter().map(|(k, c)| (*k, *c)).collect();
    assert_eq!(pairs, vec![(Key(7), 4.0), (Key(8), 2.0)]);
}
