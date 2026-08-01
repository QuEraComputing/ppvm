// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The key **hashing contract** the pass-through `HashMapStore` storage of
//! `ppvm-pauli-sum-2` relies on — tested as a *contract*, never as a raw-digest
//! diff against the old crate (their finalization folds differ by design).
//!
//! Design (`traits-2-configuration-and-hashing.md` §"Indexable values"): a key's
//! `key_hash()` is the finalized structural digest, and
//!
//!   * `Hash for Self` is exactly `state.write_u64(self.key_hash())` — so the
//!     digest reaches hashbrown untouched through `IdentityBuildHasher`;
//!   * structurally-equal keys return equal digests; and
//!   * the digest is avalanche-quality (low bits + top bits both spread), which we
//!     assert as a low-collision distribution property.

use ppvm_conformance_2::{random_pauli_string, seeded_rng};
use ppvm_pauli_sum_2::{IdentityBuildHasher, PauliWord};
use ppvm_traits_2::Indexable;

use std::collections::HashSet;
use std::hash::BuildHasher;

const SEEDS: [u64; 6] = [1, 7, 42, 123, 2024, 31337];

/// `Hash` must write **exactly** `key_hash()`: pushing the word through the
/// identity build-hasher (which returns verbatim the single `u64` written) yields
/// the finalized digest. If `Hash` wrote anything else — or more than one value —
/// `IdentityHasher::write` would panic or the value would differ.
#[test]
fn hash_writes_exactly_key_hash() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 1..=16 {
            let s = random_pauli_string(&mut rng, n);
            let w: PauliWord = PauliWord::from(s.as_str());
            let via_hasher = IdentityBuildHasher.hash_one(&w);
            assert_eq!(
                via_hasher,
                w.key_hash(),
                "Hash did not write exactly key_hash() for {s}"
            );
        }
    }
}

/// Structurally-equal keys return equal digests — independent of construction
/// path (`from(str)` twice, and a clone) and stable across the lazy cache.
#[test]
fn structurally_equal_keys_equal_digest() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 1..=16 {
            let s = random_pauli_string(&mut rng, n);
            let a: PauliWord = PauliWord::from(s.as_str());
            let b: PauliWord = PauliWord::from(s.as_str());
            assert_eq!(a, b, "equal words {s}");
            assert_eq!(a.key_hash(), b.key_hash(), "equal words, equal digest {s}");
            // Clone preserves the digest (and the cache).
            let c = a.clone();
            assert_eq!(c.key_hash(), a.key_hash(), "clone digest {s}");
        }
    }
}

/// Distinct keys almost never share a digest, and the low 12 bits (hashbrown's
/// bucket index for a moderately-sized map) spread widely — an avalanche /
/// low-collision property test over a large random key population.
#[test]
fn avalanche_low_collision() {
    let mut rng = seeded_rng(12345);
    let mut digests: HashSet<u64> = HashSet::new();
    let mut low_bits: HashSet<u64> = HashSet::new();
    let mut words: HashSet<String> = HashSet::new();
    let n = 16;

    // Draw a large population of distinct 16-qubit words.
    let target = 4000usize;
    let mut attempts = 0usize;
    while words.len() < target && attempts < target * 4 {
        attempts += 1;
        let s = random_pauli_string(&mut rng, n);
        if !words.insert(s.clone()) {
            continue;
        }
        let w: PauliWord = PauliWord::from(s.as_str());
        let h = w.key_hash();
        digests.insert(h);
        low_bits.insert(h & 0xfff); // low 12 bits
    }

    let count = words.len();
    assert!(count > 1000, "too few distinct words drawn: {count}");

    // Full-digest collisions must be vanishingly rare (allow a tiny slack for the
    // birthday bound; in practice this is 0 across all seeds tried).
    let full_collisions = count - digests.len();
    assert!(
        full_collisions <= 2,
        "too many full-digest collisions: {full_collisions} over {count} keys"
    );

    // Low 12 bits: 4096 buckets, ~count keys. A well-mixed hash fills a large
    // fraction; a hash that collapsed the low bits (e.g. no avalanche fold) would
    // fill only a handful.
    assert!(
        low_bits.len() >= 2000,
        "low bits collapsed: only {} of 4096 buckets hit by {count} keys",
        low_bits.len()
    );

    // Single-bit-difference avalanche: two words differing in one site must (with
    // overwhelming probability) differ in many digest bits, not just a few.
    let base: PauliWord<u64> = PauliWord::new(n);
    let h0 = base.key_hash();
    let mut min_flipped = u32::MAX;
    for i in 0..n {
        use ppvm_traits_2::PauliBits;
        let mut w = base.clone();
        w.set_x_bit(i, true);
        let flipped = (w.key_hash() ^ h0).count_ones();
        min_flipped = min_flipped.min(flipped);
    }
    assert!(
        min_flipped >= 8,
        "weak avalanche: a single-site flip changed only {min_flipped} digest bits"
    );
}
