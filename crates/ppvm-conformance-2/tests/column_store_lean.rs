// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Lean-oracle property tests for the **columnar** backend
//! (`Sum<ColumnStore<PauliWord, f64>, _>`).
//!
//! A storage backend is a *representation* change, so the algebra it must model
//! is the one already machine-checked for the graded map — nothing new:
//!
//! * `lean/PPVM/Algebra/GradedMap.lean` — the graded-module laws
//!   (`accumulate_comm` / `accumulate_assoc`, `scale_accumulate` / `scale_scale`,
//!   `reduce_structural`, `pushforward_eq_reset_accumulate`);
//! * `lean/PPVM/Pauli/Symplectic.lean` — `*_bijective`, so the columnar re-key
//!   (which rewrites the key planes **in place**, at the same slot) is still
//!   collision-free.
//!
//! `pauli_sum_lean.rs` already reproduces these on the hash-join backend; this
//! file re-runs them on the columnar one, because the columnar implementations
//! are *different code* for the same laws (a prefix-sum compaction rather than a
//! `HashMap::retain`, an in-place plane rewrite rather than a drain-and-reinsert).
//!
//! It also pins the ONE invariant the layout genuinely introduces, which no
//! existing lemma states because the hash map cannot violate it:
//!
//! > **Column alignment.** The key column, the coefficient column, the digest
//! > column and the bucket index describe the *same* rows: for every `i`,
//! > `get(keys[i]) == Some(coeffs[i])`, and a re-key permutes all columns
//! > consistently. Every mutation — accumulate, `scale`, `scale_by_key`, the
//! > Clifford re-key, the rotation branch merge, `retain`/`truncate`, `reduce`,
//! > the operator product — must preserve it. A missed `reindex()` after a
//! > compaction, or a plane rewritten without its digest, breaks exactly this and
//! > nothing else: the support still *iterates* correctly while every probe
//! > (`get` / `contains` / `overlap` / branch merge) silently misses.
//!
//! That is the SoA-specific hazard the phase brief flags, asserted here as an
//! executable invariant (`assert_aligned`) checked after every mutation class.

use ppvm_conformance_2::{assert_close, random_terms, seeded_rng};
use ppvm_pauli_sum_2::{
    CoefficientThreshold, ColumnPauliSum, ColumnStore, CombinedPolicy, MaxPauliWeight, NoPolicy,
    PauliWord, Policy, Sum,
};
use ppvm_traits_2::{Clifford, PauliBits, PauliError, RotationOne};

use rand::RngExt;
use rand::rngs::StdRng;
use std::collections::BTreeMap;

const SEEDS: [u64; 6] = [1, 2, 7, 42, 777, 31337];
const WIDTHS: [usize; 5] = [1, 2, 3, 5, 8];
const TOL: f64 = 1e-9;

type ColSum<P> = Sum<ColumnStore<PauliWord, f64>, P>;

/// Build a columnar sum on `n` sites from `(string, coeff)` terms, through the
/// accumulating `+=` path (colliding keys merge; zeros are kept).
fn build<P: Policy<PauliWord, f64>>(n: usize, policy: P, terms: &[(String, f64)]) -> ColSum<P> {
    let mut s = ColSum::with_policy(n, policy);
    for (w, c) in terms {
        s += (PauliWord::from(w.as_str()), *c);
    }
    s
}

fn build_plain(n: usize, terms: &[(String, f64)]) -> ColumnPauliSum {
    build(n, NoPolicy, terms)
}

/// The support as a `key → coeff` map (canonical strings).
fn support_map<P: Policy<PauliWord, f64>>(s: &ColSum<P>) -> BTreeMap<String, f64> {
    s.iter().map(|(k, c)| (k.to_string(), c)).collect()
}

#[track_caller]
fn assert_maps_close(a: &BTreeMap<String, f64>, b: &BTreeMap<String, f64>, tol: f64) {
    assert_eq!(a.len(), b.len(), "support size differs\na={a:?}\nb={b:?}");
    for (k, va) in a {
        let vb = b.get(k).unwrap_or_else(|| panic!("key {k} missing in rhs"));
        assert_close(*va, *vb, tol);
    }
}

/// A random scalar in `[-2, -0.3] ∪ [0.3, 2)`.
fn nonzero_scalar(rng: &mut StdRng) -> f64 {
    let mag = rng.random_range(0.3..2.0);
    if rng.random_range(0..2usize) == 0 {
        mag
    } else {
        -mag
    }
}

// ===========================================================================
// The SoA-specific invariant: the columns stay index-aligned under every
// mutation, and the bucket index agrees with them.
// ===========================================================================

/// Assert the columnar store's columns and index describe the same rows.
///
/// The check is deliberately *observational*, using only the public surface:
///
/// 1. every row yielded by `iter()` (the key column paired with the coefficient
///    column) is findable through the index at exactly its own coefficient
///    (`get`) — this fails if a compaction moved rows without a `reindex`, or if
///    a re-key rewrote a key plane without its digest;
/// 2. `len()` equals the number of rows iterated, and the keys are pairwise
///    distinct — the index cannot hold two slots for one key;
/// 3. `contains` agrees with `get` on every row, and a key that is NOT in the
///    support is not found (probed with a key perturbed off the support), so the
///    index has no stale slot pointing at a freed row.
#[track_caller]
fn assert_aligned<P: Policy<PauliWord, f64>>(s: &ColSum<P>, label: &str) {
    let rows: Vec<(PauliWord, f64)> = s.iter().collect();
    assert_eq!(
        rows.len(),
        s.len(),
        "[{label}] len() {} disagrees with the number of iterated rows {}",
        s.len(),
        rows.len()
    );
    let mut seen: BTreeMap<String, f64> = BTreeMap::new();
    for (k, c) in &rows {
        let probed = s.get(k).unwrap_or_else(|| {
            panic!("[{label}] row {k} is in the columns but the index cannot find it")
        });
        assert_eq!(
            probed.to_bits(),
            c.to_bits(),
            "[{label}] key column and coefficient column are misaligned at {k}: \
             iter says {c}, the index resolves {probed}"
        );
        assert!(
            s.contains(k, c),
            "[{label}] contains() disagrees with the columns at {k}"
        );
        assert!(
            seen.insert(k.to_string(), *c).is_none(),
            "[{label}] duplicate row {k}: the index holds two slots for one key"
        );
    }
    // A key outside the support must not resolve — a stale index slot would.
    if let Some((k, _)) = rows.first() {
        let mut absent = k.clone();
        // Flip one site's X bit; skip the check if that lands on another live key.
        let bit = absent.x_bit(0);
        absent.set_x_bit(0, !bit);
        if !seen.contains_key(&absent.to_string()) {
            assert!(
                s.get(&absent).is_none(),
                "[{label}] the index resolved {absent}, which is not in the columns"
            );
        }
    }
}

/// Every mutation class, each followed by the alignment check.
#[test]
fn columns_stay_index_aligned_under_every_mutation() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let terms = random_terms(&mut rng, n, 40);
            let mut s = build(
                n,
                CombinedPolicy(
                    CoefficientThreshold { threshold: 0.5 },
                    MaxPauliWeight(n.max(1)),
                ),
                &terms,
            );
            assert_aligned(&s, "after accumulate");

            // L2 scale: a pure coefficient pass — must not move a row.
            s.scale(&nonzero_scalar(&mut rng));
            assert_aligned(&s, "after scale");

            // `scale_by_key` (the diagonal channel path).
            for q in 0..n {
                s.pauli_error(q, [1e-3, 2e-3, 3e-3]);
            }
            assert_aligned(&s, "after pauli_error");

            // Sign flip by key (`x`/`y`/`z` — in place, no re-key).
            s.x(0);
            s.z(n - 1);
            assert_aligned(&s, "after sign flips");

            // The in-place re-key: the plane rewrite that must permute EVERY
            // column consistently (and refresh the digests + the index).
            for q in 0..n {
                s.h(q);
                assert_aligned(&s, "after h");
            }
            for q in 0..n.saturating_sub(1) {
                s.cnot(q, q + 1);
                assert_aligned(&s, "after cnot");
                s.cz(q, q + 1);
                assert_aligned(&s, "after cz");
            }

            // The rotation branch merge: appends rows AND probes existing ones.
            for q in 0..n {
                s.rx(q, 0.31);
                assert_aligned(&s, "after rx");
                s.rz(q, 0.71);
                assert_aligned(&s, "after rz");
            }

            // The retain compaction (truncate) — the operation that MOVES rows.
            s.truncate();
            assert_aligned(&s, "after truncate");

            // The structural reduce (the prefix-sum compaction proper).
            s.reduce();
            assert_aligned(&s, "after reduce");

            // …and the store still works as a support afterwards: another gate
            // and another probe round trip.
            s.h(0);
            assert_aligned(&s, "after post-compaction re-key");
        }
    }
}

/// The compaction moves rows: a support whose *middle* rows are dropped is the
/// case that breaks a missing `reindex`. Built deliberately so the survivors are
/// interleaved with the dropped rows rather than being a prefix.
#[test]
fn compaction_moves_rows_and_keeps_the_index_consistent() {
    // Width 6 → 4^6 = 4096 distinct words; the base-4 digits of `k` give 64
    // PAIRWISE DISTINCT keys, so every term is its own row (a generator that
    // repeated keys would merge them and, with the coefficients accumulating
    // above the floor, silently stop exercising the compaction at all).
    let n = 6usize;
    for &drop_every in &[2usize, 3, 5] {
        let terms: Vec<(String, f64)> = (0..64usize)
            .map(|k| {
                let s: String = (0..n)
                    .map(|j| ['I', 'X', 'Y', 'Z'][(k >> (2 * j)) & 0b11])
                    .collect();
                // Below the floor on the dropped rows, above it on the survivors.
                let c = if k % drop_every == 0 {
                    1e-9
                } else {
                    1.0 + k as f64
                };
                (s, c)
            })
            .collect();
        let mut s = build(n, CoefficientThreshold { threshold: 1e-6 }, &terms);
        let before = support_map(&s);
        assert_eq!(
            before.len(),
            terms.len(),
            "the generator must produce pairwise-distinct keys"
        );
        assert_aligned(&s, "pre-compaction");
        s.truncate();
        assert_aligned(&s, "post-compaction");
        // The compaction must actually have MOVED rows — otherwise a missing
        // `reindex` after it would be untested (`write == len` short-circuits).
        let dropped = before.len() - s.len();
        assert!(
            dropped > 0 && !s.is_empty(),
            "the cell must drop some rows and keep some: dropped {dropped} of {}",
            before.len()
        );

        // Exactly the sub-threshold rows are gone, and every survivor kept its
        // coefficient bit-for-bit (truncation only removes).
        for (k, c) in &before {
            let key = PauliWord::from(k.as_str());
            if c.abs() >= 1e-6 {
                assert_eq!(
                    s.get(&key).map(f64::to_bits),
                    Some(c.to_bits()),
                    "survivor {k} lost or changed after the compaction"
                );
            } else {
                assert!(s.get(&key).is_none(), "sub-threshold {k} survived");
            }
        }
    }
}

// ===========================================================================
// GradedMap.lean — the graded-module laws, on the columnar backend.
// ===========================================================================

/// `accumulate_comm` / `accumulate_assoc`: accumulation is order-independent.
#[test]
fn accumulate_comm_and_assoc_on_the_columnar_backend() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let terms = random_terms(&mut rng, n, 30);
            let base = build_plain(n, &terms);

            let mut rev = terms.clone();
            rev.reverse();
            assert_maps_close(
                &support_map(&base),
                &support_map(&build_plain(n, &rev)),
                TOL,
            );

            let k = terms.len() / 3;
            let (a, rest) = terms.split_at(k);
            let (b, c) = rest.split_at(k);
            let mut regrouped: Vec<(String, f64)> = Vec::new();
            regrouped.extend_from_slice(c);
            regrouped.extend_from_slice(a);
            regrouped.extend_from_slice(b);
            assert_maps_close(
                &support_map(&base),
                &support_map(&build_plain(n, &regrouped)),
                TOL,
            );
        }
    }
}

/// `reduce_structural`: after `reduce`, a key is in the support iff its
/// coefficient is nonzero — here realized by the prefix-sum compaction.
#[test]
fn reduce_structural_on_the_columnar_backend() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let terms = random_terms(&mut rng, n, 30);
            let mut s = build_plain(n, &terms);
            s.reduce();
            for (_k, c) in s.iter() {
                assert!(c != 0.0, "reduced support contains a zero coefficient");
            }
            assert_aligned(&s, "reduce_structural forward");

            let key = "X".repeat(n);
            let mut with_cancel: Vec<(String, f64)> =
                terms.iter().filter(|(w, _)| *w != key).cloned().collect();
            with_cancel.push((key.clone(), 0.75));
            with_cancel.push((key.clone(), -0.75));
            let mut s2 = build_plain(n, &with_cancel);
            // Caller-driven: the cancellation is live until `reduce` runs.
            assert!(
                s2.contains_key(&PauliWord::from(key.as_str())),
                "the cancelled key must survive until reduce (contract 2)"
            );
            s2.reduce();
            assert!(
                !s2.contains_key(&PauliWord::from(key.as_str())),
                "cancelled key {key} survived reduce"
            );
            assert_aligned(&s2, "reduce_structural backward");
        }
    }
}

/// `scale_accumulate` (`scale` distributes over `+`) and `scale_scale`
/// (`scale s ∘ scale t == scale (s·t)`).
#[test]
fn scale_laws_on_the_columnar_backend() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let a = random_terms(&mut rng, n, 20);
            let b = random_terms(&mut rng, n, 20);
            let s = nonzero_scalar(&mut rng);
            let t = nonzero_scalar(&mut rng);

            // scale(a + b) == scale(a) + scale(b)
            let mut joint: Vec<(String, f64)> = a.clone();
            joint.extend_from_slice(&b);
            let mut lhs = build_plain(n, &joint);
            lhs.scale(&s);

            let mut sa = build_plain(n, &a);
            sa.scale(&s);
            let mut sb = build_plain(n, &b);
            sb.scale(&s);
            let rhs = &sa + &sb;
            assert_maps_close(&support_map(&lhs), &support_map(&rhs), TOL);

            // scale(t)∘scale(s) == scale(s·t)
            let mut twice = build_plain(n, &joint);
            twice.scale(&s);
            twice.scale(&t);
            let mut once = build_plain(n, &joint);
            once.scale(&(s * t));
            assert_maps_close(&support_map(&twice), &support_map(&once), 1e-8);
        }
    }
}

/// The `Pair` overlap is biadditive, symmetric and homogeneous in each slot.
#[test]
fn overlap_is_bilinear_and_symmetric_on_the_columnar_backend() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let a = build_plain(n, &random_terms(&mut rng, n, 20));
            let b = build_plain(n, &random_terms(&mut rng, n, 20));
            let c = build_plain(n, &random_terms(&mut rng, n, 20));
            let s = nonzero_scalar(&mut rng);

            let bar = |x: f64| 1e-9 * x.abs().max(1.0);

            // Symmetry.
            let ab = a.overlap(&b);
            assert_close(ab, b.overlap(&a), bar(ab));

            // Left-additivity: ⟨a + b, c⟩ == ⟨a, c⟩ + ⟨b, c⟩.
            let sum_ab = &a + &b;
            let lhs = sum_ab.overlap(&c);
            let rhs = a.overlap(&c) + b.overlap(&c);
            assert_close(lhs, rhs, bar(lhs));

            // Homogeneity: ⟨s·a, b⟩ == s·⟨a, b⟩.
            let mut sa = a.clone();
            sa.scale(&s);
            assert_close(sa.overlap(&b), s * ab, bar(s * ab));
        }
    }
}

// ===========================================================================
// Symplectic.lean — `*_bijective`: the columnar re-key is collision-free.
// ===========================================================================

/// The in-place plane rewrite is a **bijection on keys**: the support size never
/// changes and no two rows ever collapse onto one. This is sharper on the
/// columnar backend than on the hash map: the hash map would silently *merge* a
/// collision (an accumulate), whereas the columnar re-key writes at the same
/// slot, so a non-injective map would leave two rows with the same key — which
/// `assert_aligned`'s distinctness check catches directly.
#[test]
fn clifford_rekey_is_a_bijection_on_the_columnar_backend() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[2usize, 3, 5, 8] {
            let terms = random_terms(&mut rng, n, 30);
            let mut s = build_plain(n, &terms);
            let before = s.len();
            let keys_before: Vec<String> = support_map(&s).keys().cloned().collect();

            for i in 0..20 {
                let a = rng.random_range(0..n);
                let b = rng.random_range(0..n);
                match i % 5 {
                    0 => s.h(a),
                    1 => s.s(a),
                    2 if a != b => s.cnot(a, b),
                    3 if a != b => s.cz(a, b),
                    _ => s.z(a),
                }
                assert_eq!(s.len(), before, "gate {i} changed the support size");
                assert_eq!(
                    support_map(&s).len(),
                    before,
                    "gate {i} collapsed two keys into one"
                );
                assert_aligned(&s, "after a Clifford re-key");
            }

            // H;H is the identity re-key, keys and all.
            let mut t = build_plain(n, &terms);
            for q in 0..n {
                t.h(q);
                t.h(q);
            }
            let keys_after: Vec<String> = support_map(&t).keys().cloned().collect();
            assert_eq!(keys_before, keys_after, "H;H is not the identity re-key");
        }
    }
}

/// Exhaustive over the finite single-qubit basis: every Clifford is a
/// permutation of `{I, X, Y, Z}` on the columnar backend, and the whole basis
/// carried at once stays four distinct, individually findable rows.
#[test]
fn single_qubit_cliffords_permute_the_basis_exhaustively() {
    const BASIS: [&str; 4] = ["I", "X", "Y", "Z"];
    for gate in ["h", "s", "x", "y", "z"] {
        let terms: Vec<(String, f64)> = BASIS
            .iter()
            .enumerate()
            .map(|(i, p)| (p.to_string(), 1.0 + i as f64))
            .collect();
        let mut s = build_plain(1, &terms);
        match gate {
            "h" => s.h(0),
            "s" => s.s(0),
            "x" => s.x(0),
            "y" => s.y(0),
            _ => s.z(0),
        }
        assert_eq!(s.len(), 4, "{gate} is not a permutation of the basis");
        let m = support_map(&s);
        assert_eq!(m.len(), 4, "{gate} collapsed two basis elements");
        for p in BASIS {
            assert!(m.contains_key(p), "{gate} lost the {p} row");
        }
        assert_aligned(&s, gate);
        // The coefficients are permuted (up to sign), never summed.
        let mags: Vec<f64> = {
            let mut v: Vec<f64> = m.values().map(|c| c.abs()).collect();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            v
        };
        assert_eq!(mags, vec![1.0, 2.0, 3.0, 4.0], "{gate} merged coefficients");
    }
}

/// Exhaustive single-qubit `scale` / `reduce`, mirroring the hash-backend test:
/// `scale(0)` keeps the row (contract 2) and only `reduce` removes it.
#[test]
fn single_qubit_exhaustive_scale_and_reduce_on_the_columnar_backend() {
    for p in ["I", "X", "Y", "Z"] {
        let s = build_plain(1, &[(p.to_string(), 3.0)]);
        let mut zeroed = s.clone();
        zeroed.scale(&0.0);
        assert_eq!(
            zeroed.len(),
            1,
            "scale must not remove the zeroed {p} term (contract 2)"
        );
        assert_aligned(&zeroed, "scale by zero");
        zeroed.reduce();
        assert!(
            zeroed.is_empty(),
            "scale-by-0 then reduce should empty the {p} sum"
        );

        let mut scaled = s.clone();
        scaled.scale(&-2.0);
        assert_close(scaled.get(&PauliWord::from(p)).unwrap(), -6.0, TOL);
    }
}
