// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Lean-oracle property tests for `ppvm-pauli-sum-2::PauliSum<f64>`: the graded
//! `C[K]` module/algebra laws machine-checked in
//! `lean/PPVM/Algebra/GradedMap.lean` and the truncation bounds of
//! `lean/PPVM/Algebra/Truncation.lean`, reproduced as randomized (and, where the
//! single-qubit case is finite, exhaustive) property tests on the sum.
//!
//! Laws reproduced:
//! * `accumulate_comm` / `accumulate_assoc` — accumulation is order-independent;
//! * `reduce_structural` — a key is in support iff its coefficient is nonzero;
//! * `scale_accumulate` — `scale` distributes over accumulation;
//! * `scale_scale` — `scale s ∘ scale t == scale (s·t)`;
//! * `overlap_add_left` / `overlap_add_right` — the `Pair` overlap is biadditive;
//! * `overlap_comm` — and symmetric;
//! * `overlap_smul_left` / `overlap_smul_right` — and `C`-homogeneous in each slot;
//! * `Truncation.l1_bound` — `|Σ_{dropped} c_k e_k| ≤ Σ_{dropped}|c_k|` for
//!   per-key expectations `|e_k| ≤ 1`;
//! * `Truncation.cutoff_mismatch` — `PauliSum` keeps at `|c| == threshold`
//!   (`≥`), whereas the tableau's strict `>` would drop it;
//! * the Clifford re-key is a **support-preserving bijection** (pushforward): it
//!   preserves `len` and maps each key to a distinct key.
//!
//! The phase/conjugation oracles these lean on are already validated for the word
//! (`phased_pauli_word_diff.rs` / `_lean.rs`); here we assert the *sum-level*
//! consequences.

use ppvm_conformance_2::{assert_close, build_new_sum, random_terms, seeded_rng};
use ppvm_pauli_sum_2::{CoefficientThreshold, MaxPauliWeight, PauliSum, PauliWord, Policy};
use ppvm_traits_2::{Clifford, Word};

use rand::RngExt;
use rand::rngs::StdRng;
use std::collections::BTreeMap;

const SEEDS: [u64; 10] = [1, 2, 3, 7, 42, 99, 123, 777, 2024, 31337];
const WIDTHS: [usize; 6] = [1, 2, 3, 5, 8, 12];
const TOL: f64 = 1e-9;

/// The sum's support as a `key → coeff` map (canonical strings), over any policy.
fn support_map<P>(s: &PauliSum<f64, P>) -> BTreeMap<String, f64>
where
    P: Policy<PauliWord, f64>,
{
    s.iter().map(|(k, c)| (k.to_string(), c)).collect()
}

/// Assert two support maps agree within `tol`.
#[track_caller]
fn assert_maps_close(a: &BTreeMap<String, f64>, b: &BTreeMap<String, f64>, tol: f64) {
    assert_eq!(a.len(), b.len(), "support size differs\na={a:?}\nb={b:?}");
    for (k, va) in a {
        let vb = b.get(k).unwrap_or_else(|| panic!("key {k} missing in rhs"));
        assert_close(*va, *vb, tol);
    }
}

/// A random scalar in `[-2, 2)`, avoiding the exact-zero band so `scale` cannot
/// synthesize spurious zeros the laws would then have to reduce away.
fn nonzero_scalar(rng: &mut StdRng) -> f64 {
    let mag = rng.random_range(0.3..2.0);
    if rng.random_range(0..2usize) == 0 {
        mag
    } else {
        -mag
    }
}

// ---------------------------------------------------------------------------
// accumulate_comm / accumulate_assoc — order independence of the module `+`.
// ---------------------------------------------------------------------------

#[test]
fn accumulate_comm_and_assoc() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let terms = random_terms(&mut rng, n, 30);
            let base = build_new_sum(n, &terms);

            // Commutativity: reversed insertion order → identical support.
            let mut rev = terms.clone();
            rev.reverse();
            assert_maps_close(
                &support_map(&base),
                &support_map(&build_new_sum(n, &rev)),
                TOL,
            );

            // Associativity: partition into three blocks and rebuild in a permuted
            // block order; `+` is associative so the merged support is invariant.
            let k = terms.len() / 3;
            let (a, rest) = terms.split_at(k);
            let (b, c) = rest.split_at(k);
            let mut regrouped: Vec<(String, f64)> = Vec::new();
            regrouped.extend_from_slice(c);
            regrouped.extend_from_slice(a);
            regrouped.extend_from_slice(b);
            assert_maps_close(
                &support_map(&base),
                &support_map(&build_new_sum(n, &regrouped)),
                TOL,
            );
        }
    }
}

// ---------------------------------------------------------------------------
// reduce_structural — a key is in support iff its coefficient is nonzero.
// ---------------------------------------------------------------------------

#[test]
fn reduce_structural() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let terms = random_terms(&mut rng, n, 30);
            // `reduce` is caller-driven (construction keeps zeros, as in old), so
            // realize it before asserting the reduced-form invariant.
            let mut s = build_new_sum(n, &terms);
            s.reduce();
            // Forward: every key in support has a nonzero coefficient.
            for (_w, c) in s.iter() {
                assert!(c != 0.0, "reduced support contains a zero coefficient");
            }
            // Backward: a key whose coefficients cancel to 0 is absent.
            let key = "X".repeat(n);
            let mut with_cancel: Vec<(String, f64)> =
                terms.iter().filter(|(w, _)| *w != key).cloned().collect();
            with_cancel.push((key.clone(), 0.75));
            with_cancel.push((key.clone(), -0.75));
            let mut s2 = build_new_sum(n, &with_cancel);
            s2.reduce();
            assert!(
                !s2.contains_key(&PauliWord::from(key.as_str())),
                "cancelled key {key} survived reduce"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// scale_accumulate / scale_scale — the L2 module action.
// ---------------------------------------------------------------------------

#[test]
fn scale_accumulate_and_scale_scale() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let terms = random_terms(&mut rng, n, 25);
            let s = nonzero_scalar(&mut rng);
            let t = nonzero_scalar(&mut rng);

            // scale_accumulate: scaling the whole sum == pre-scaling each term.
            let mut scaled = build_new_sum(n, &terms);
            scaled.scale(&s);
            let prescaled = build_new_sum(
                n,
                &terms
                    .iter()
                    .map(|(w, c)| (w.clone(), c * s))
                    .collect::<Vec<_>>(),
            );
            assert_maps_close(&support_map(&scaled), &support_map(&prescaled), TOL);

            // scale_scale: scale s ∘ scale t == scale (s·t).
            let mut st = build_new_sum(n, &terms);
            st.scale(&s);
            st.scale(&t);
            let mut once = build_new_sum(n, &terms);
            once.scale(&(s * t));
            assert_maps_close(
                &support_map(&st),
                &support_map(&once),
                TOL.max((s * t).abs() * 1e-12),
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Pair overlap — biadditive, symmetric, homogeneous in each slot.
// ---------------------------------------------------------------------------

#[test]
fn overlap_bilinear_symmetric_homogeneous() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let a_terms = random_terms(&mut rng, n, 12);
            let b_terms = random_terms(&mut rng, n, 12);
            let g_terms = random_terms(&mut rng, n, 12);

            let a = build_new_sum(n, &a_terms);
            let b = build_new_sum(n, &b_terms);
            let g = build_new_sum(n, &g_terms);

            // overlap_comm.
            let ab = a.overlap(&b);
            let tol = TOL.max(ab.abs() * 1e-10);
            assert_close(ab, b.overlap(&a), tol);

            // overlap_add_left: ⟨a+b, g⟩ == ⟨a,g⟩ + ⟨b,g⟩.
            let mut ab_terms = a_terms.clone();
            ab_terms.extend_from_slice(&b_terms);
            let a_plus_b = build_new_sum(n, &ab_terms);
            let lhs = a_plus_b.overlap(&g);
            let rhs = a.overlap(&g) + b.overlap(&g);
            assert_close(lhs, rhs, TOL.max(lhs.abs() * 1e-10));

            // overlap_add_right: ⟨a, b+g⟩ == ⟨a,b⟩ + ⟨a,g⟩.
            let mut bg_terms = b_terms.clone();
            bg_terms.extend_from_slice(&g_terms);
            let b_plus_g = build_new_sum(n, &bg_terms);
            let lhs = a.overlap(&b_plus_g);
            let rhs = a.overlap(&b) + a.overlap(&g);
            assert_close(lhs, rhs, TOL.max(lhs.abs() * 1e-10));

            // overlap_smul_left / _right: ⟨s·a, b⟩ == s·⟨a,b⟩ == ⟨a, s·b⟩.
            let s = nonzero_scalar(&mut rng);
            let mut sa = build_new_sum(n, &a_terms);
            sa.scale(&s);
            let mut sb = build_new_sum(n, &b_terms);
            sb.scale(&s);
            let target = s * ab;
            let htol = TOL.max(target.abs() * 1e-10);
            assert_close(sa.overlap(&b), target, htol);
            assert_close(a.overlap(&sb), target, htol);
        }
    }
}

// ---------------------------------------------------------------------------
// Truncation.l1_bound — dropped observable error ≤ dropped ℓ¹ mass.
// ---------------------------------------------------------------------------

#[test]
fn truncation_l1_bound() {
    type WeightSum = PauliSum<f64, MaxPauliWeight>;
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[3usize, 5, 8, 12] {
            let terms = random_terms(&mut rng, n, 40);
            // Cap weight aggressively so a nonempty set is dropped.
            let cap = (n / 2).max(1);
            let build = |pol: MaxPauliWeight| -> WeightSum {
                WeightSum::from_terms_with_policy(
                    n,
                    pol,
                    terms.iter().map(|(w, c)| (PauliWord::from(w.as_str()), *c)),
                )
            };

            let full = build(MaxPauliWeight(usize::MAX));
            let mut kept = build(MaxPauliWeight(cap));
            kept.truncate();

            // Dropped set = full keys absent from the truncated support.
            let kept_map = support_map(&kept);
            let mut dropped: Vec<(String, f64)> = Vec::new();
            for (w, c) in full.iter() {
                let key = w.to_string();
                if !kept_map.contains_key(&key) {
                    dropped.push((key, c));
                }
            }

            // Per-key expectations e_k ∈ [-1, 1]; the incurred observable error is
            // Σ_{dropped} c_k e_k, bounded by the ℓ¹ mass Σ_{dropped}|c_k|.
            let mut error = 0.0f64;
            let mut l1 = 0.0f64;
            for (i, (_w, c)) in dropped.iter().enumerate() {
                let e = ((i as f64 * 0.37 + seed as f64 * 0.11).sin()).clamp(-1.0, 1.0);
                error += c * e;
                l1 += c.abs();
            }
            assert!(
                error.abs() <= l1 + TOL,
                "l1_bound violated: |error| {} > Σ|c| {} (n {n} seed {seed})",
                error.abs(),
                l1
            );

            // Also confirm the truncation actually kept only weight ≤ cap.
            for (w, _c) in kept.iter() {
                assert!(w.weight() <= cap, "kept a term above the weight cap");
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Truncation.cutoff_mismatch — `≥` (PauliSum) keeps at |c| == threshold; the
// tableau's strict `>` would drop it.
// ---------------------------------------------------------------------------

#[test]
fn truncation_cutoff_mismatch() {
    type ThreshSum = PauliSum<f64, CoefficientThreshold>;
    for &threshold in &[0.25f64, 0.5, 1.0, 2.0] {
        let n = 4;
        // One term exactly at the threshold, one strictly above, one strictly below.
        let terms = [
            ("XXII".to_string(), threshold),       // exactly at
            ("ZZII".to_string(), threshold * 2.0), // above
            ("YIII".to_string(), threshold * 0.5), // below
        ];
        let mut s: ThreshSum = ThreshSum::from_terms_with_policy(
            n,
            CoefficientThreshold { threshold },
            terms.iter().map(|(w, c)| (PauliWord::from(w.as_str()), *c)),
        );
        s.truncate();
        let map = support_map(&s);

        // PauliSum keep-rule is `|c| >= threshold`: the boundary term is KEPT.
        assert!(
            map.contains_key("XXII"),
            "PauliSum (>=) dropped the |c| == threshold {threshold} term"
        );
        assert!(map.contains_key("ZZII"), "above-threshold term dropped");
        assert!(!map.contains_key("YIII"), "below-threshold term kept");

        // The tableau's strict `>` keep-rule would DROP the boundary term — this is
        // exactly the Lean `cutoff_mismatch`: at |c| == t, the `>=` keep-rule
        // (`t <= |t|`) and the `>` keep-rule (`t < |t|`) disagree.
        let t = threshold;
        assert_ne!(
            t <= t.abs(),
            t < t.abs(),
            "cutoff_mismatch at threshold {threshold}"
        );
    }
}

// ---------------------------------------------------------------------------
// Clifford re-key is a support-preserving bijection (pushforward).
// ---------------------------------------------------------------------------

#[test]
fn clifford_rekey_is_support_preserving_bijection() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[2usize, 3, 5, 8, 12] {
            let terms = random_terms(&mut rng, n, 30);
            let mut s = build_new_sum(n, &terms);
            let before = s.len();
            let keys_before: Vec<String> = support_map(&s).keys().cloned().collect();

            // Apply a spread of single- and two-qubit Cliffords.
            let gates: Vec<(usize, usize)> = (0..20)
                .map(|_| (rng.random_range(0..n), rng.random_range(0..n)))
                .collect();
            for (i, &(a, b)) in gates.iter().enumerate() {
                match i % 5 {
                    0 => s.h(a),
                    1 => s.s(a),
                    2 if a != b => s.cnot(a, b),
                    3 if a != b => s.cz(a, b),
                    _ => s.z(a),
                }
                // Bijection: len is preserved by every gate...
                assert_eq!(s.len(), before, "gate {i} changed the support size");
                // ...and the keys stay pairwise distinct (a genuine bijection, no
                // two source keys collapsing onto one target).
                let m = support_map(&s);
                assert_eq!(m.len(), before, "gate {i} collapsed two keys into one");
            }

            // A full round-trip identity check: H;H on every qubit returns to the
            // original support (H is an involution), keys and all.
            let mut t = build_new_sum(n, &terms);
            for q in 0..n {
                t.h(q);
                t.h(q);
            }
            let keys_after: Vec<String> = support_map(&t).keys().cloned().collect();
            assert_eq!(keys_before, keys_after, "H;H is not the identity re-key");
        }
    }
}

// ---------------------------------------------------------------------------
// Exhaustive single-qubit reduce_structural + scale over the finite Pauli basis.
// ---------------------------------------------------------------------------

#[test]
fn single_qubit_exhaustive_scale_and_reduce() {
    for p in ["I", "X", "Y", "Z"] {
        // scale by 0 reduces to empty (every coefficient becomes an exact zero,
        // which `scale` leaves in place — `reduce` is the explicit canonicalizer).
        let s = build_new_sum(1, &[(p.to_string(), 3.0)]);
        let mut zeroed = s.clone();
        zeroed.scale(&0.0);
        assert_eq!(
            zeroed.len(),
            1,
            "scale must not remove the zeroed {p} term (contract 2)"
        );
        zeroed.reduce();
        assert!(
            zeroed.is_empty(),
            "scale-by-0 then reduce should empty the {p} sum"
        );

        // scale by a nonzero constant is exact on a single term.
        let mut scaled = s.clone();
        scaled.scale(&-2.0);
        assert_close(scaled.get(&PauliWord::from(p)).unwrap(), -6.0, TOL);
    }
}
