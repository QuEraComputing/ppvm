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

/// An `n`-site word of weight 1–4 — the shape a propagated observable actually
/// holds (`Σ Z_i` fanned out by a few gates), and the population that exposes a
/// weak low-bit avalanche: almost every storage byte is zero, so the hasher sees
/// a nearly constant input and only the few set bits can decorrelate the buckets.
fn low_weight_word(rng: &mut rand::rngs::StdRng, n: usize) -> String {
    use rand::RngExt;
    let weight = rng.random_range(1..5usize).min(n);
    let mut sites = vec!['I'; n];
    for _ in 0..weight {
        let i = rng.random_range(0..n);
        sites[i] = ['X', 'Y', 'Z'][rng.random_range(0..3usize)];
    }
    sites.into_iter().collect()
}

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
            let via_hasher = IdentityBuildHasher.hash_one(w);
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
            let c = a;
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
        let mut w = base;
        w.set_x_bit(i, true);
        let flipped = (w.key_hash() ^ h0).count_ones();
        min_flipped = min_flipped.min(flipped);
    }
    assert!(
        min_flipped >= 8,
        "weak avalanche: a single-site flip changed only {min_flipped} digest bits"
    );
}

/// The **storage-tier fold** must hold on every `[u8; N]` width, not just the
/// default word.
///
/// This is the regression detector for the architecture feature the old crate
/// documents in `examples/trotter_storage_cliff.rs`: raw fxhash correlates in its
/// low bits for short keys, which clusters the bucket index at high map fill and
/// cost old a 4–5× *storage-tier cliff* until `HashFinalize::finalize_hash(h,
/// size_of::<A>())` folded the high bits down. `avalanche_low_collision` above
/// only exercises the default width, so a fold that silently became a no-op for,
/// say, `[u8; 8]` would pass it.
///
/// The low bits matter to BOTH backends and are the *only* thing either one uses
/// to pick a bucket: `HashMapStore` feeds the digest to hashbrown verbatim
/// through `IdentityBuildHasher`, and `ColumnStore` masks it itself
/// (`bucket = hash & mask`). So the property is asserted on the digest's low bits
/// at each tier's *full* width — the fill regime where old cliffed.
///
/// The key population is deliberately **not** uniform-random: it is the
/// *low-weight* population a real propagated observable holds (a Trotter support
/// is `Σ Z_i` fanned out to a few sites, so almost every plane byte is zero).
/// Uniform-random 64-bit planes spread even without the fold — it is exactly the
/// sparse, structured population that exposes fxhash's low-bit correlation, which
/// is why old's cliff showed up in `trotter_storage_cliff.rs` and not in a
/// random-key microbench.
#[test]
fn the_finalize_fold_survives_every_storage_tier() {
    /// Draw `TARGET` distinct `$n`-site words of one storage width and report
    /// `(distinct digests, distinct low-12-bit buckets, population)`. A macro
    /// rather than a generic function: the storage bound lives in a private
    /// module of `ppvm-pauli-word-2`, so the width can only be named at a use
    /// site.
    macro_rules! spread {
        ($storage:ty, $n:expr) => {{
            let mut rng = seeded_rng(0xC0FFEE);
            let mut words: HashSet<String> = HashSet::new();
            let mut digests: HashSet<u64> = HashSet::new();
            let mut buckets: HashSet<u64> = HashSet::new();
            let mut attempts = 0usize;
            while words.len() < TARGET && attempts < TARGET * 40 {
                attempts += 1;
                let s = low_weight_word(&mut rng, $n);
                if !words.insert(s.clone()) {
                    continue;
                }
                let h = PauliWord::<$storage>::from(s.as_str()).key_hash();
                digests.insert(h);
                buckets.insert(h & 0xfff);
            }
            (digests.len(), buckets.len(), words.len())
        }};
    }

    const TARGET: usize = 4000;
    // (tier label, full site width for that storage)
    let tiers: [(&str, (usize, usize, usize)); 4] = [
        ("[u8; 8]", spread!([u8; 8], 64)),
        ("[u8; 16]", spread!([u8; 16], 128)),
        ("[u8; 32]", spread!([u8; 32], 256)),
        ("u64", spread!(u64, 64)),
    ];

    for (label, (digests, buckets, count)) in tiers {
        // Reported so `--nocapture` shows the margin, not just pass/fail.
        println!("[{label}] {count} keys → {digests} digests, {buckets}/4096 buckets");
        assert!(
            count as f64 >= TARGET as f64 * 0.9,
            "[{label}] too few distinct words: {count}"
        );
        // Full-digest collisions stay rare. Note this is NOT the fold's job: the
        // fold `raw ^ (raw >> 32)` is a bijection on `u64`, so it moves no
        // full-digest collision either way (measured: 3 either side on this
        // population — fxhash's own rate on sparse input). It is asserted only as
        // a sanity bound on the underlying hash.
        let full_collisions = count - digests;
        assert!(
            full_collisions * 200 <= count,
            "[{label}] too many full-digest collisions: {full_collisions} over {count} keys"
        );
        // THE detector. 4096 buckets, ~4000 sparse keys: with the fold the digest
        // hits ~2300 of them; with the fold disabled it hits 589 — a 3.9×
        // clustering of hashbrown's / `ColumnStore`'s bucket index, which is old's
        // documented 4–5× storage-tier cliff. The threshold sits between the two
        // measured regimes with margin on both sides.
        assert!(
            buckets >= 1500,
            "[{label}] the finalize fold collapsed: only {buckets} of 4096 buckets \
             hit by {count} keys — this is the storage-tier cliff"
        );
    }
}
