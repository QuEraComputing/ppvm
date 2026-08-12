// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Controlled A/B that attributes the ONE op class where the columnar backend
//! loses to the hash-join one: the single-qubit rotation (`column_store/rx`,
//! `pauli_sum/rx` — `ColumnStore` ≈ 1.2× `HashMapStore`, while every other op
//! class is 0.38×–0.76×).
//!
//! **Result** (three runs, this machine, `[u8; 8]` / `f64` / fxhash on both
//! sides): the rotation itself is *not* the cost. Splitting the two passes and
//! then the merge regime:
//!
//! | regime | column/hash |
//! | --- | --- |
//! | commuting (pass 1 only, no branch) | 1.04–1.12 |
//! | anticommuting, steady state (branches merge onto existing rows) | **0.87** |
//! | the microbench's mixed population | 1.01–1.02 |
//! | untruncated growth sweep (branches append NEW rows) | **1.19–1.29** |
//!
//! So the columnar rotation *wins* when branches merge and loses only where they
//! **insert**: appending a row is a push onto every bit plane plus an index
//! insert (and a `reindex` at the load factor), against one hashbrown insert.
//! `column_store/build/from_terms` (1.29×, pure insertion) is the independent
//! corroboration. The criterion `rx` groups are dominated by that insert path
//! because their support grows monotonically.
//!
//! `RotateInPlace` is two passes:
//!
//! 1. **pass 1** — walk the coefficient column in place, materialising the key of
//!    each row (`keys.get(i)`, i.e. an SoA → AoS gather across the bit planes for
//!    the columnar store, versus a plain `&K` for the hash map) and calling the
//!    rotation kernel;
//! 2. **pass 2** — merge each produced branch term through the index
//!    (`find` + accumulate, or append).
//!
//! The single variable held here is **whether pass 1 produces any branch at
//! all**, with pass 1's row count, key width, coefficient type, capacity hint and
//! policy identical on both arms:
//!
//! * `commuting`   — every key carries `I` or `X` on the rotation site, so it
//!   commutes with `X` and the kernel returns `None`: pass 1 runs over the whole
//!   support, pass 2 is empty;
//! * `anticommuting` — every key carries `Y` or `Z` on the rotation site, so every
//!   row produces a branch: pass 1 is the same walk, pass 2 merges `|support|`
//!   terms (all colliding, since the `Y`↔`Z` partner of every key is also in the
//!   support — steady state, no appends, no resize).
//!
//! The difference of the two column/hash ratios is therefore the branch-merge
//! cost alone. Interleaved min-of-N in ONE process (same build, same thermal
//! state), so the layout bias cancels.
//!
//! Run: `cargo run --release -p ppvm-conformance-2 --example column_rx_attrib`

use std::time::{Duration, Instant};

use ppvm_conformance_2::seeded_rng;
use ppvm_pauli_sum_2::{ColumnStore, HashMapStore, NoPolicy, PauliWord, Sum};
use ppvm_traits_2::RotationOne;

use rand::RngExt;
use rand::rngs::StdRng;
use std::collections::BTreeSet;

type Key = PauliWord<[u8; 8]>;
type HashSum = Sum<HashMapStore<Key, f64>, NoPolicy>;
type ColSum = Sum<ColumnStore<Key, f64>, NoPolicy>;

const N: usize = 16;
const TERMS: usize = 1000;
const REPS: usize = 200;
const ROUNDS: usize = 25;

/// `count` distinct Pauli strings on `n` sites whose site-0 letter is drawn from
/// `site0` (the only difference between the two arms).
fn terms(seed: u64, count: usize, n: usize, site0: &[char]) -> Vec<(String, f64)> {
    let mut rng: StdRng = seeded_rng(seed);
    let mut seen: BTreeSet<String> = BTreeSet::new();
    let mut out = Vec::with_capacity(count);
    while out.len() < count {
        let mut w: String = String::with_capacity(n);
        w.push(site0[rng.random_range(0..site0.len())]);
        for _ in 1..n {
            w.push(['I', 'X', 'Y', 'Z'][rng.random_range(0..4usize)]);
        }
        if seen.insert(w.clone()) {
            out.push((w, rng.random_range(0.25..2.0)));
        }
    }
    out
}

fn hash_sum(t: &[(String, f64)]) -> HashSum {
    let mut s = HashSum::with_capacity(N, NoPolicy, N * 10);
    for (w, c) in t {
        s += (Key::from(w.as_str()), *c);
    }
    s
}

fn col_sum(t: &[(String, f64)]) -> ColSum {
    let mut s = ColSum::with_capacity(N, NoPolicy, N * 10);
    for (w, c) in t {
        s += (Key::from(w.as_str()), *c);
    }
    s
}

/// Min-of-`ROUNDS` wall clock for `REPS` back-to-back `rx` calls.
fn time(mut f: impl FnMut()) -> Duration {
    // Warm up: also lets the anticommuting arm reach its steady support (the
    // branch partners are all present, so no arm is timed while appending).
    for _ in 0..REPS {
        f();
    }
    let mut best = Duration::MAX;
    for _ in 0..ROUNDS {
        let t0 = Instant::now();
        for _ in 0..REPS {
            f();
        }
        best = best.min(t0.elapsed());
    }
    best
}

fn main() {
    // `commuting`: site 0 ∈ {I, X} → `rx(0, θ)` produces no branch.
    // `anticommuting`: site 0 ∈ {Y, Z} → every row produces one.
    for (label, site0) in [
        ("commuting   (pass 1 only)", ['I', 'X'].as_slice()),
        (
            "anticommuting (pass 1 + branch merge)",
            ['Y', 'Z'].as_slice(),
        ),
        (
            "mixed (the bench's uniform site 0)",
            ['I', 'X', 'Y', 'Z'].as_slice(),
        ),
    ] {
        let t = terms(8, TERMS, N, site0);
        let mut h = hash_sum(&t);
        let mut c = col_sum(&t);

        // Interleave the two arms so any drift hits both equally.
        let mut hb = Duration::MAX;
        let mut cb = Duration::MAX;
        for _ in 0..2 {
            hb = hb.min(time(|| h.rx(0, 0.1)));
            cb = cb.min(time(|| c.rx(0, 0.1)));
        }
        let (hn, cn) = (h.len(), c.len());
        assert_eq!(hn, cn, "the two backends must hold the same support");
        let per = |d: Duration| d.as_secs_f64() * 1e9 / REPS as f64;
        println!(
            "{label:<40}  support {hn:>5}  hash {:>9.1} ns  column {:>9.1} ns  column/hash {:.3}",
            per(hb),
            per(cb),
            per(cb) / per(hb),
        );
    }

    growth_arm();
}

/// The third regime, and the one the `pauli_sum/rx` bench actually measures: an
/// **untruncated sweep** over every site, where each rotation appends the branch
/// keys that are not yet in the support. Steady state (above) merges into rows
/// that already exist; here the support grows ~2^sites, so the timed work is
/// dominated by *insertion* — a hashbrown insert versus the columnar append
/// (push onto each bit plane, then an index insert and, at the load factor, a
/// full `reindex`).
///
/// Same single-variable discipline: identical seed support, identical gate
/// sequence, the clone kept OUT of the timed region (the two backends clone
/// differently, and that is not what is being attributed).
fn growth_arm() {
    const SITES: usize = 8;
    let t = terms(8, TERMS, N, ['I', 'X', 'Y', 'Z'].as_slice());
    let hseed = hash_sum(&t);
    let cseed = col_sum(&t);

    let mut hb = Duration::MAX;
    let mut cb = Duration::MAX;
    let (mut hn, mut cn) = (0usize, 0usize);
    for _ in 0..ROUNDS {
        let mut h = hseed.clone();
        let t0 = Instant::now();
        for i in 0..SITES {
            h.rx(i, 0.1);
        }
        hb = hb.min(t0.elapsed());
        hn = h.len();

        let mut c = cseed.clone();
        let t0 = Instant::now();
        for i in 0..SITES {
            c.rx(i, 0.1);
        }
        cb = cb.min(t0.elapsed());
        cn = c.len();
    }
    assert_eq!(hn, cn, "the two backends must hold the same support");
    println!(
        "{:<40}  support {hn:>5}  hash {:>9.1} ns  column {:>9.1} ns  column/hash {:.3}",
        format!("growth (untruncated {SITES}-site sweep)"),
        hb.as_secs_f64() * 1e9,
        cb.as_secs_f64() * 1e9,
        cb.as_secs_f64() / hb.as_secs_f64(),
    );
}
