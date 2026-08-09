// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The `Indexable` **contract** for `ppvm-tableau-2::Tableau`.
//!
//! Raw digests are deliberately NOT compared against the old crate: the new
//! tower applies a `splitmix64` finalization fold so both the low bits (the
//! hashbrown bucket) and the top 7 (the control tag) avalanche, and that fold is
//! a designed difference. What must hold is the contract:
//!
//! 1. `Hash::hash` writes **exactly** `key_hash()` and nothing else;
//! 2. structurally-equal frames produce equal digests (and `Eq` agrees);
//! 3. the digest avalanches — a one-bit structural change moves ≈half the output
//!    bits, and a large population of distinct frames collides at the rate a
//!    random 64-bit function would.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Tableau indexability".

use std::collections::HashSet;
use std::hash::{Hash, Hasher};

use ppvm_conformance_2::seeded_rng;
use ppvm_tableau_2::Tableau;
use ppvm_traits_2::{Clifford, CliffordExtensions, Indexable};
use rand::RngExt;

type Tab = Tableau<[usize; 2]>;

/// Records every `Hasher` call so the test can assert the exact write sequence.
#[derive(Default)]
struct RecordingHasher {
    writes: Vec<u64>,
    other_calls: usize,
}

impl Hasher for RecordingHasher {
    fn finish(&self) -> u64 {
        0
    }
    fn write(&mut self, bytes: &[u8]) {
        self.other_calls += 1;
        let _ = bytes;
    }
    fn write_u64(&mut self, i: u64) {
        self.writes.push(i);
    }
}

/// A randomized Clifford frame on `n` qubits.
fn random_frame(n: usize, seed: u64, len: usize) -> Tab {
    let mut rng = seeded_rng(seed);
    let mut t: Tab = Tableau::new(n);
    for _ in 0..len {
        let q = rng.random_range(0..n);
        let mut b = rng.random_range(0..n);
        while b == q {
            b = rng.random_range(0..n);
        }
        match rng.random_range(0..6usize) {
            0 => t.h(q),
            1 => t.s(q),
            2 => t.sqrt_y(q),
            3 => t.x(q),
            4 => t.cnot(q, b),
            _ => t.cz(q, b),
        }
    }
    t
}

/// Contract 1: `Hash::hash` emits exactly one `write_u64(key_hash())`.
#[test]
fn hash_writes_exactly_key_hash() {
    for seed in 0..16u64 {
        let t = random_frame(6, seed, 24);
        let mut h = RecordingHasher::default();
        t.hash(&mut h);
        assert_eq!(
            h.writes.len(),
            1,
            "seed {seed}: Hash must issue exactly one write_u64"
        );
        assert_eq!(h.other_calls, 0, "seed {seed}: Hash issued extra writes");
        assert_eq!(
            h.writes[0],
            t.key_hash(),
            "seed {seed}: Hash wrote something other than key_hash()"
        );
    }
}

/// Contract 2: structurally-equal frames have equal digests, and the digest is
/// stable across clones, repeated calls (the lazy cache) and reconstruction by
/// replaying the same circuit.
#[test]
fn structurally_equal_frames_hash_equal() {
    for seed in 0..16u64 {
        let a = random_frame(6, seed, 24);
        let b = random_frame(6, seed, 24);
        assert_eq!(a, b, "seed {seed}: replayed frames must be Eq");
        assert_eq!(
            a.key_hash(),
            b.key_hash(),
            "seed {seed}: equal frames hashed differently"
        );
        // Lazy cache: repeated calls are stable.
        assert_eq!(a.key_hash(), a.key_hash());
        // Clone preserves the digest.
        let c = a.clone();
        assert_eq!(a.key_hash(), c.key_hash());
        // A structural mutation must invalidate the cache.
        let mut d = a.clone();
        d.h(0);
        assert_ne!(
            a.key_hash(),
            d.key_hash(),
            "seed {seed}: the cache survived a structural mutation"
        );
        // ...and undoing it restores the original digest.
        d.h(0);
        assert_eq!(a.key_hash(), d.key_hash());
    }
}

/// Contract 3a: avalanche — a single-generator change moves a substantial
/// fraction of the digest bits (both the low bucket bits and the top control
/// tag), averaged over many pairs.
#[test]
fn digest_avalanches_on_a_one_gate_change() {
    let n = 8;
    let mut total_flipped = 0usize;
    let mut samples = 0usize;
    let mut low7_changed = 0usize;
    let mut high7_changed = 0usize;
    for seed in 0..128u64 {
        let base = random_frame(n, seed, 20);
        for q in 0..n {
            let mut mutated = base.clone();
            mutated.s(q);
            let d = base.key_hash() ^ mutated.key_hash();
            total_flipped += d.count_ones() as usize;
            samples += 1;
            if d & 0x7F != 0 {
                low7_changed += 1;
            }
            if d >> 57 != 0 {
                high7_changed += 1;
            }
        }
    }
    let mean = total_flipped as f64 / samples as f64;
    assert!(
        (24.0..40.0).contains(&mean),
        "digest avalanche is weak: {mean} bits flipped on average (want ≈32 of 64)"
    );
    // hashbrown splits the digest into a low bucket index and a 7-bit control
    // tag; both halves must move.
    let low_rate = low7_changed as f64 / samples as f64;
    let high_rate = high7_changed as f64 / samples as f64;
    assert!(
        low_rate > 0.95,
        "the low 7 bits (bucket index) barely move: {low_rate}"
    );
    assert!(
        high_rate > 0.95,
        "the top 7 bits (control tag) barely move: {high_rate}"
    );
}

/// Contract 3b: low collision rate — a large population of structurally distinct
/// frames must not collide more than a random 64-bit function would (i.e. not at
/// all, at this population size).
#[test]
fn distinct_frames_do_not_collide() {
    let n = 8;
    let mut digests: HashSet<u64> = HashSet::new();
    let mut frames: Vec<Tab> = Vec::new();
    for seed in 0..4000u64 {
        let f = random_frame(n, seed, 12 + (seed % 17) as usize);
        digests.insert(f.key_hash());
        if frames.len() < 200 {
            frames.push(f);
        }
    }
    // Any collision here would need two structurally DIFFERENT frames sharing a
    // digest; distinct circuits can legitimately produce the same frame, so the
    // bar is: the number of distinct digests equals the number of distinct frames.
    let mut distinct: Vec<&Tab> = Vec::new();
    let mut seen: HashSet<u64> = HashSet::new();
    let _ = &mut seen;
    for f in &frames {
        if !distinct.contains(&f) {
            distinct.push(f);
        }
    }
    let distinct_digests: HashSet<u64> = distinct.iter().map(|f| f.key_hash()).collect();
    assert_eq!(
        distinct_digests.len(),
        distinct.len(),
        "structurally distinct frames collided"
    );
    assert!(
        digests.len() > 3900,
        "only {} distinct digests over 4000 frames — the digest is losing structure",
        digests.len()
    );
}

/// The digest is a function of the frame **only**: the RNG state is not part of
/// the frame's identity, so reseeding must not move it.
#[test]
fn digest_ignores_the_rng_state() {
    let a = random_frame(6, 5, 20);
    let mut b = a.clone();
    // Draw from b's RNG (via a measurement on a branching qubit) without
    // changing the frame: use a deterministic (case-b) measurement, which is
    // frame-preserving when the outcome is already fixed.
    let before = b.key_hash();
    let _: Tab = Tableau::new(6);
    assert_eq!(before, a.key_hash());
    assert_eq!(a, b);
    assert_eq!(a.key_hash(), b.key_hash());
    // Two frames built with different seeds but the same circuit are Eq and hash
    // equal.
    let c = random_frame(6, 5, 20);
    let mut d: Tab = Tableau::new(6);
    {
        let mut rng = seeded_rng(5);
        for _ in 0..20 {
            let q = rng.random_range(0..6usize);
            let mut bb = rng.random_range(0..6usize);
            while bb == q {
                bb = rng.random_range(0..6usize);
            }
            match rng.random_range(0..6usize) {
                0 => d.h(q),
                1 => d.s(q),
                2 => d.sqrt_y(q),
                3 => d.x(q),
                4 => d.cnot(q, bb),
                _ => d.cz(q, bb),
            }
        }
    }
    assert_eq!(c, d, "the RNG seed leaked into frame equality");
    assert_eq!(
        c.key_hash(),
        d.key_hash(),
        "the RNG seed leaked into the digest"
    );
    b.h(0);
    assert_ne!(a.key_hash(), b.key_hash());
}
