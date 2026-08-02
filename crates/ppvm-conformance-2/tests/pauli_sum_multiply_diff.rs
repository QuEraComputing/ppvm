// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential + Lean-oracle coverage for the **L4 operator product** (the
//! twisted convolution) on `ppvm-pauli-sum-2::Sum`.
//!
//! The product runs over `Complex<f64>`, not `f64`: a Pauli product emits an
//! `iᵏ` phase, so old bounded `Mul`/`MulAssign` on `ComplexCoefficient` and the
//! new crate bounds `multiply_into` on `ImaginaryUnit` — `f64` implements
//! neither, so a real-coefficient sum has **no** product method on either side
//! (behavioural contract 12). Both engines are therefore pinned to
//! `Complex<f64>` here, with `[u8; 8]` storage on the old side.
//!
//! Two very different obligations live in this file:
//!
//! 1. **`A *= P` (single Pauli word) — a strict differential vs old.** Old's
//!    `impl MulAssign<PauliWord> for PauliSum` (`ppvm-pauli-sum/src/sum/ops.rs:95`)
//!    is a *single* `map_add` over a bijection and is believed correct, so the new
//!    [`Sum::mul_word_assign`] must reproduce it key-for-key.
//!
//! 2. **`A * B` (sum × sum) — NOT a diff against old; old is wrong.** Old's
//!    `impl MulAssign<PauliSum<T>> for PauliSum<T>`
//!    (`ppvm-pauli-sum/src/sum/ops.rs:70`) calls `self.map_add(..)` once per rhs
//!    term, and `map_add` *replaces* the support with its image — so it computes
//!    the product **chain** `A·b₀P₀·b₁P₁` instead of the bilinear sum
//!    `A·b₀P₀ + A·b₁P₁`. It is non-bilinear for any rhs with more than one term and
//!    is untested in old. The oracle is therefore the Lean `twistedConv`
//!    (`lean/PPVM/Algebra/Twisted.lean`), which is biadditive by construction; the
//!    tests below assert the **correct** value, built out of old's *trustworthy*
//!    single-word path, and separately pin that old's sum×sum disagrees with it.

use std::collections::BTreeMap;

use num::Complex;
use ppvm_conformance_2::{random_pauli_string, seeded_rng};
use ppvm_pauli_sum::config::fxhash::Byte;
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_pauli_sum_2::{PauliSum as NewPauliSum, PauliWord as NewWord};
use ppvm_pauli_word::word::PauliWord as OldWordT;
use ppvm_traits::traits::{
    Clifford as OldClifford, NoStrategy, PauliError as OldPauliError, RotationOne as OldRotationOne,
};
use rand::RngExt;
use rand::rngs::StdRng;

/// The OLD reference complex sum: `[u8; 8]` storage (64-qubit capacity), FxHash,
/// no truncation strategy — the product must not truncate on either side.
type OldCSum = OldPauliSum<Byte<8, Complex<f64>, NoStrategy>>;
/// The OLD word matching `OldCSum`'s storage.
type OldCWord = OldWordT<[u8; 8]>;
/// The NEW complex sum under test (`NoPolicy`).
type NewCSum = NewPauliSum<Complex<f64>>;

const SEEDS: [u64; 6] = [1, 7, 42, 99, 2024, 31337];
const WIDTHS: [usize; 4] = [2, 3, 4, 5];
/// Coefficient tolerance. The product is a handful of multiply-adds per key, so
/// `1e-12` is far above the achievable error and far below any real difference.
const TOL: f64 = 1e-12;

// ---------------------------------------------------------------------------
// Harness
// ---------------------------------------------------------------------------

fn random_complex_terms(rng: &mut StdRng, n: usize, count: usize) -> Vec<(String, Complex<f64>)> {
    (0..count)
        .map(|_| {
            let w = random_pauli_string(rng, n);
            let re = rng.random_range(-2.0..2.0f64);
            let im = rng.random_range(-2.0..2.0f64);
            (w, Complex::new(re, im))
        })
        .collect()
}

fn build_old(n: usize, terms: &[(String, Complex<f64>)]) -> OldCSum {
    let mut s: OldCSum = OldPauliSum::builder().n_qubits(n).build();
    for (w, c) in terms {
        s += (w.as_str(), *c);
    }
    s
}

fn build_new(n: usize, terms: &[(String, Complex<f64>)]) -> NewCSum {
    NewCSum::from_terms(
        n,
        terms.iter().map(|(w, c)| (NewWord::from(w.as_str()), *c)),
    )
}

fn old_support(s: &OldCSum) -> BTreeMap<String, Complex<f64>> {
    s.data().iter().map(|(k, c)| (k.to_string(), *c)).collect()
}

fn new_support(s: &NewCSum) -> BTreeMap<String, Complex<f64>> {
    s.iter().map(|(k, c)| (k.to_string(), c)).collect()
}

#[track_caller]
fn assert_maps_match(
    expected: &BTreeMap<String, Complex<f64>>,
    got: &BTreeMap<String, Complex<f64>>,
    what: &str,
) {
    assert_eq!(
        expected.keys().collect::<Vec<_>>(),
        got.keys().collect::<Vec<_>>(),
        "{what}: key sets differ"
    );
    for (k, e) in expected {
        let g = got[k];
        assert!((e - g).norm() <= TOL, "{what}: {k} -> {e} vs {g}");
    }
}

/// Whether two supports agree — the boolean form of [`assert_maps_match`], for
/// the test that asserts old *disagrees*.
fn maps_match(
    a: &BTreeMap<String, Complex<f64>>,
    b: &BTreeMap<String, Complex<f64>>,
    tol: f64,
) -> bool {
    a.len() == b.len()
        && a.iter()
            .all(|(k, v)| b.get(k).is_some_and(|w| (v - w).norm() <= tol))
}

/// The **Lean-correct** product `A · B`, assembled entirely out of OLD
/// primitives that are believed correct: for each rhs monomial `(q, b)`, run old's
/// bijective `A *= q` (`ops.rs:95`) and accumulate `b · (A·q)` into a map.
///
/// This is exactly `twistedConv` (`lean/PPVM/Algebra/Twisted.lean`) unrolled over
/// the rhs support — biadditive in `B` by construction — so it is a differential
/// reference that does *not* inherit old's sum×sum bug.
fn old_reference_product(
    a: &OldCSum,
    rhs_terms: &[(String, Complex<f64>)],
) -> BTreeMap<String, Complex<f64>> {
    let mut acc: BTreeMap<String, Complex<f64>> = BTreeMap::new();
    for (qstr, b) in rhs_terms {
        let mut t = a.clone();
        t *= OldCWord::from(qstr.as_str());
        for (k, c) in t.data().iter() {
            *acc.entry(k.to_string()).or_insert(Complex::new(0.0, 0.0)) += c * b;
        }
    }
    acc
}

// ---------------------------------------------------------------------------
// 1. `A *= P` — strict differential against old's bijective single-word path
// ---------------------------------------------------------------------------

#[test]
fn mul_by_single_word_matches_old() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            for count in [1usize, 5, 20] {
                let terms = random_complex_terms(&mut rng, n, count);
                let q = random_pauli_string(&mut rng, n);

                let mut old = build_old(n, &terms);
                old *= OldCWord::from(q.as_str());

                let mut new = build_new(n, &terms);
                new.mul_word_assign(&NewWord::from(q.as_str()));

                assert_maps_match(
                    &old_support(&old),
                    &new_support(&new),
                    &format!("A *= {q} (seed {seed}, n {n}, {count} terms)"),
                );
            }
        }
    }
}

/// The single-word product is a **bijection**: it neither grows nor shrinks the
/// support, and (like old's `map_add`, which replaces the support wholesale) it
/// truncates nothing and drops nothing — including a term whose coefficient is far
/// below any plausible threshold.
#[test]
fn mul_by_single_word_preserves_support_size() {
    let n = 3;
    let terms = vec![
        ("XYZ".to_string(), Complex::new(1.5, -0.5)),
        ("IZI".to_string(), Complex::new(1e-30, 0.0)),
        ("ZZZ".to_string(), Complex::new(-2.0, 0.25)),
    ];
    let mut old = build_old(n, &terms);
    let mut new = build_new(n, &terms);

    old *= OldCWord::from("YZX");
    new.mul_word_assign(&NewWord::from("YZX"));

    assert_eq!(old.len(), terms.len(), "old: bijection preserves |support|");
    assert_eq!(new.len(), terms.len(), "new: bijection preserves |support|");
    assert_maps_match(&old_support(&old), &new_support(&new), "A *= YZX");
}

// ---------------------------------------------------------------------------
// 2. `A * B` — the documented divergence from old
// ---------------------------------------------------------------------------

/// **THE ONE ALLOWED BEHAVIOUR DIVERGENCE.**
///
/// `A · (b₀P₀ + b₁P₁)` must be the bilinear sum `A·b₀P₀ + A·b₁P₁`. Old
/// (`ppvm-pauli-sum/src/sum/ops.rs:70`) instead folds each rhs monomial back into
/// `self` via `map_add`, which *replaces* the support, so it computes the product
/// **chain** `A·b₀P₀·b₁P₁`.
///
/// Lean oracle: `twistedConv` in `lean/PPVM/Algebra/Twisted.lean` (the outer
/// product with the `i^{phaseExpN p q}` twist), whose monomial case is `tmul` and
/// which is biadditive in each argument by construction — the property old
/// violates. The reference value below is built from old's *own* believed-correct
/// single-word path, so this is a genuine differential assertion of the CORRECT
/// value, not a self-consistency check.
#[test]
fn sum_product_is_bilinear_where_old_computes_a_chain() {
    let n = 2;
    let a_terms = vec![
        ("XZ".to_string(), Complex::new(1.0, 0.0)),
        ("IY".to_string(), Complex::new(-0.5, 0.25)),
    ];
    let b_terms = vec![
        ("ZI".to_string(), Complex::new(2.0, 0.0)),
        ("IX".to_string(), Complex::new(0.0, 3.0)),
    ];

    let a_old = build_old(n, &a_terms);
    let a_new = build_new(n, &a_terms);
    let b_new = build_new(n, &b_terms);

    // The Lean-correct value, assembled from old's bijective single-word path.
    let reference = old_reference_product(&a_old, &b_terms);

    // NEW: accumulates every monomial product into a fresh accumulator.
    let got = new_support(&a_new.multiply(&b_new));
    assert_maps_match(&reference, &got, "A·(b₀P₀ + b₁P₁) (Lean twistedConv)");

    // OLD cannot be run here — and that is the second half of the finding. Its
    // `impl MulAssign<PauliSum<T>> for PauliSum<T>` (ops.rs:70) carries the bound
    //     PhasedPauliWord<..>: for<'a> From<&'a T::PauliWordType>
    // while `PhasedPauliWord`'s only word conversion is `impl<W: PauliWordTrait>
    // From<W>` (`ppvm-pauli-word/src/phase/data.rs:114`) and `PauliWordTrait` is
    // implemented **only** for `PauliWord` itself, never for a reference
    // (`ppvm-pauli-word/src/word/data.rs:100`). The bound is therefore
    // unsatisfiable for every shipped `Config`: `old *= other_sum` does not
    // compile, so old's sum×sum is not merely untested, it is UNREACHABLE dead
    // code. (Only `Mul`/`MulAssign<PauliWord>`, exercised above, are live.) That
    // is why the reference is assembled from the single-word path instead, and why
    // there is no old-vs-new numeric diff for this shape.
    //
    // Chain-vs-sum is still observable in the new engine, so pin that the correct
    // (bilinear) answer is NOT old's intended `A·b₀P₀·b₁P₁` fold:
    let chain = new_support(&a_new.multiply(&b_new_first()).multiply(&b_new_second()));
    assert!(
        !maps_match(&reference, &chain, 1e-9),
        "the bilinear product must differ from the product chain ops.rs:70 computes"
    );
}

/// `b₀P₀` alone — the first rhs monomial of
/// [`sum_product_is_bilinear_where_old_computes_a_chain`].
fn b_new_first() -> NewCSum {
    build_new(2, &[("ZI".to_string(), Complex::new(2.0, 0.0))])
}

/// `b₁P₁` alone — the second rhs monomial.
fn b_new_second() -> NewCSum {
    build_new(2, &[("IX".to_string(), Complex::new(0.0, 3.0))])
}

/// Bilinearity in the rhs on random operands: `A·(B + C) == A·B + A·C`, with the
/// reference again built from old's single-word path. (`twistedConv` biadditivity,
/// `lean/PPVM/Algebra/Twisted.lean`.)
#[test]
fn sum_product_matches_the_monomial_expansion() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let a_terms = random_complex_terms(&mut rng, n, 6);
            // Distinct rhs keys: a repeated key would be merged by `build_*` but
            // double-counted by the unrolled reference.
            let mut b_terms = random_complex_terms(&mut rng, n, 4);
            b_terms.sort_by(|x, y| x.0.cmp(&y.0));
            b_terms.dedup_by(|x, y| x.0 == y.0);

            let a_old = build_old(n, &a_terms);
            let a_new = build_new(n, &a_terms);
            let b_new = build_new(n, &b_terms);

            let reference = old_reference_product(&a_old, &b_terms);
            let mut got = a_new.multiply(&b_new);
            // The reference map keeps no exact-zero keys only if none arise; the
            // product deliberately does keep them, so canonicalize both sides.
            got.reduce();
            let reference: BTreeMap<_, _> = reference
                .into_iter()
                .filter(|(_, c)| c.norm() > 0.0)
                .collect();

            assert_maps_match(
                &reference,
                &new_support(&got),
                &format!("A·B monomial expansion (seed {seed}, n {n})"),
            );
        }
    }
}

// ---------------------------------------------------------------------------
// 3. Lean-oracle algebra laws (no old counterpart exists — that absence is
//    itself a finding; `Twisted.lean` is the oracle)
// ---------------------------------------------------------------------------

/// `tmul_assoc` / `gtmul_assoc`: the twisted product is associative over any
/// commutative ring with `i⁴ = 1`.
#[test]
fn product_is_associative() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let a = build_new(n, &random_complex_terms(&mut rng, n, 5));
            let b = build_new(n, &random_complex_terms(&mut rng, n, 4));
            let c = build_new(n, &random_complex_terms(&mut rng, n, 3));

            let mut left = a.multiply(&b).multiply(&c);
            let mut right = a.multiply(&b.multiply(&c));
            left.reduce();
            right.reduce();
            assert_maps_match(
                &new_support(&left),
                &new_support(&right),
                &format!("(A·B)·C vs A·(B·C) (seed {seed}, n {n})"),
            );
        }
    }
}

/// `one_tmul` / `tmul_one`: the all-identity word is the unit.
#[test]
fn identity_is_the_unit() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let a = build_new(n, &random_complex_terms(&mut rng, n, 6));
            let id = build_new(n, &[("I".repeat(n), Complex::new(1.0, 0.0))]);
            assert_maps_match(&new_support(&a), &new_support(&a.multiply(&id)), "A·I");
            assert_maps_match(&new_support(&a), &new_support(&id.multiply(&a)), "I·A");
        }
    }
}

/// `twistedConv_apply_id`: `(A·B)[I] == ⟨A, B⟩` — the L4↔L3 tie that licenses
/// calling `Pair::overlap` the Hilbert–Schmidt pairing.
#[test]
fn identity_coefficient_is_the_overlap() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let a = build_new(n, &random_complex_terms(&mut rng, n, 8));
            let b = build_new(n, &random_complex_terms(&mut rng, n, 8));
            let product = a.multiply(&b);
            let id = product
                .get(&NewWord::from("I".repeat(n).as_str()))
                .unwrap_or(Complex::new(0.0, 0.0));
            let ov = a.overlap(&b);
            assert!(
                (id - ov).norm() <= TOL,
                "(A·B)[I] {id} != overlap {ov} (seed {seed}, n {n})"
            );
        }
    }
}

/// `phaseExpN_self`: `P·P = +I`, so a single-term sum squares to `c²·I`.
#[test]
fn word_squares_to_identity() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let w = random_pauli_string(&mut rng, n);
            let c = Complex::new(1.5, -0.5);
            let p = build_new(n, &[(w.clone(), c)]);
            let sq = p.multiply(&p);
            assert_eq!(sq.len(), 1, "{w}² support");
            let got = sq.get(&NewWord::from("I".repeat(n).as_str())).unwrap();
            assert!((got - c * c).norm() <= TOL, "{w}²: {got} vs {}", c * c);
        }
    }
}

/// The product must **not** truncate and must **not** drop exact-zero
/// cancellations (behavioural contracts 1 and 2).
#[test]
fn product_neither_truncates_nor_drops_zeros() {
    let n = 1;
    let a = build_new(n, &[("X".to_string(), Complex::new(1.0, 0.0))]);
    let mut acc = NewCSum::new(n);
    a.multiply_into(
        &build_new(n, &[("X".to_string(), Complex::new(1.0, 0.0))]),
        &mut acc,
    );
    a.multiply_into(
        &build_new(n, &[("X".to_string(), Complex::new(-1.0, 0.0))]),
        &mut acc,
    );
    assert_eq!(acc.len(), 1, "the cancelled identity key must survive");
    assert_eq!(acc.get(&NewWord::from("I")), Some(Complex::new(0.0, 0.0)));

    // A tiny coefficient survives too: no policy runs inside the product.
    let tiny = build_new(n, &[("Z".to_string(), Complex::new(1e-30, 0.0))]);
    let prod = tiny.multiply(&build_new(n, &[("Z".to_string(), Complex::new(1.0, 0.0))]));
    assert_eq!(prod.len(), 1);
    assert!(prod.get(&NewWord::from("I")).unwrap().norm() > 0.0);
}

// ---------------------------------------------------------------------------
// 4. Integration workload 6 — "observable-product-and-variance", at real scale
//
// The shapes above run on hand-sized supports. The baseline workload is the
// *observable* one: build `O` by propagating the noisy-TFIM Trotter circuit until
// it carries a realistic support, then (i) right-multiply by a single Pauli word
// and (ii) square it for the ⟨O²⟩ variance estimate.
// ---------------------------------------------------------------------------

/// The OLD `f64` sum used only to GROW a realistic support. `[u8; 8]` storage and
/// the headline `CoefficientThreshold(1e-6)` floor, matching
/// `tests/pauli_sum_integration_diff.rs`.
type OldGrowSum = OldPauliSum<
    ppvm_pauli_sum::config::fxhash::ByteF64<8, ppvm_pauli_sum::strategy::CoefficientThreshold>,
>;

/// Propagate `Σᵢ Zᵢ` through the noisy-TFIM Trotter circuit on the OLD engine and
/// return its final key set (sorted). Grown on the *old* side deliberately: the
/// keys the product is then benchmarked/diffed over come from an engine that is
/// not the one under test.
fn grown_keys(n: usize, steps: usize) -> Vec<String> {
    let dt = 0.1_f64;
    let theta_x = dt;
    // `J = 1.0` (the example's support-growing form) rather than the headline
    // `1/8`: the point here is a large, realistic support for the product.
    let theta_zz = dt;
    let noise = [1e-4 / 4.0; 3];

    let mut state: OldGrowSum = OldPauliSum::builder()
        .n_qubits(n)
        .capacity(n.pow(2))
        .strategy(ppvm_pauli_sum::strategy::CoefficientThreshold(1e-6))
        .build();
    for i in 0..n {
        let s: String = (0..n).map(|j| if j == i { 'Z' } else { 'I' }).collect();
        state += (s.as_str(), 1.0);
    }
    for _ in 0..steps {
        for i in 0..n {
            state.pauli_error(i, noise);
            state.truncate();
            state.rx(i, theta_x);
            state.truncate();
        }
        for i in 0..n - 1 {
            state.pauli_error(i + 1, noise);
            state.truncate();
            state.pauli_error(i, noise);
            state.truncate();
            state.cnot(i, i + 1);
            state.rz(i + 1, theta_zz);
            state.cnot(i, i + 1);
            state.truncate();
        }
    }
    let mut keys: Vec<String> = state.data().keys().map(|k| k.to_string()).collect();
    keys.sort();
    keys
}

/// Deterministic complex coefficients over a key list — no RNG, so old and new
/// see byte-identical input.
fn deterministic_complex(keys: &[String]) -> Vec<(String, Complex<f64>)> {
    keys.iter()
        .enumerate()
        .map(|(i, k)| {
            let t = (i % 17) as f64;
            (k.clone(), Complex::new(0.5 + t / 32.0, -0.25 + t / 64.0))
        })
        .collect()
}

/// **Workload 6 (i) at scale.** `O *= P` over a realistic multi-hundred-term
/// observable, strictly diffed against old's bijective `MulAssign<PauliWord>`:
/// identical key sets, coefficients within `1e-12`, support size preserved
/// exactly (the re-key is a bijection and drops nothing — injectivity of
/// `p ↦ p·q` is machine-checked as `mulWord_right_injective` in
/// `lean/PPVM/Pauli/Word.lean`, which is what licenses the plain-`insert`
/// `RekeyBijective` path where a collision would silently drop a term).
#[test]
fn mul_by_single_word_matches_old_on_a_grown_observable() {
    let n = 10;
    let keys = grown_keys(n, 8);
    assert!(
        keys.len() >= 200,
        "the grown observable is too small to be representative: {} terms",
        keys.len()
    );
    let terms = deterministic_complex(&keys);

    let mut old = build_old(n, &terms);
    let mut new = build_new(n, &terms);
    assert_eq!(old.len(), terms.len());
    assert_eq!(new.len(), terms.len());

    // A word that is non-identity on every site, so no term takes a trivial path.
    let word: String = (0..n).map(|i| ['X', 'Y', 'Z', 'Y'][i % 4]).collect();
    old *= OldCWord::from(word.as_str());
    new.mul_word_assign(&NewWord::from(word.as_str()));

    assert_eq!(old.len(), terms.len(), "old: bijection preserves |support|");
    assert_eq!(new.len(), terms.len(), "new: bijection preserves |support|");
    assert_maps_match(
        &old_support(&old),
        &new_support(&new),
        &format!("grown O *= {word}"),
    );
}

/// **Workload 6 (ii).** The ⟨O²⟩ variance shape on the baseline's stated scale
/// (3–5 qubits, 10–50 terms), checked against the Lean laws rather than old
/// (whose sum×sum is unreachable — see the module docs):
/// bilinearity, associativity, `A·I = A`, and `(A·B)[I] = ⟨A, B⟩`.
#[test]
fn observable_squaring_satisfies_the_lean_laws_at_workload_scale() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 3..=5usize {
            for &count in &[10usize, 25, 50] {
                let a_terms = random_complex_terms(&mut rng, n, count);
                let b_terms = random_complex_terms(&mut rng, n, count);
                let d_terms = random_complex_terms(&mut rng, n, count);
                let a = build_new(n, &a_terms);
                let b = build_new(n, &b_terms);
                let d = build_new(n, &d_terms);
                let id = build_new(n, &[("I".repeat(n), Complex::new(1.0, 0.0))]);

                // A·I == A.
                assert_maps_match(
                    &new_support(&a),
                    &new_support(&a.multiply(&id)),
                    &format!("A·I (seed {seed}, n {n}, {count} terms)"),
                );

                // Bilinearity in the rhs, via two `multiply_into` calls.
                let mut b_plus_d = NewCSum::new(n);
                a.multiply_into(&b, &mut b_plus_d);
                a.multiply_into(&d, &mut b_plus_d);
                // `B + D` as one sum: `from_terms` runs `accumulate_batch`, which
                // genuinely *adds* colliding keys (unlike old's `AddAssign<
                // PauliSum>`, which routes through `Extend` and overwrites).
                let summed = build_new(n, &[b_terms.clone(), d_terms.clone()].concat());
                let mut lhs = a.multiply(&summed);
                let mut rhs = b_plus_d;
                lhs.reduce();
                rhs.reduce();
                assert_maps_match(
                    &new_support(&lhs),
                    &new_support(&rhs),
                    &format!("A·(B+D) (seed {seed}, n {n}, {count} terms)"),
                );

                // Associativity.
                let mut left = a.multiply(&b).multiply(&d);
                let mut right = a.multiply(&b.multiply(&d));
                left.reduce();
                right.reduce();
                assert_maps_match(
                    &new_support(&left),
                    &new_support(&right),
                    &format!("(A·B)·D (seed {seed}, n {n}, {count} terms)"),
                );

                // ⟨O²⟩: the identity coefficient of `A·A` is `⟨A, A⟩`.
                let sq = a.multiply(&a);
                let got = sq
                    .get(&NewWord::from("I".repeat(n).as_str()))
                    .unwrap_or(Complex::new(0.0, 0.0));
                let ov = a.overlap(&a);
                assert!(
                    (got - ov).norm() <= TOL,
                    "⟨O²⟩ {got} != overlap {ov} (seed {seed}, n {n}, {count} terms)"
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 5. Behavioural contract 1, extended to the NEW operation
//
// "Gates never auto-truncate — truncation is caller-driven only." The tests
// above run the product under `NoPolicy`, which cannot demonstrate this: a
// policy that would drop nothing proves nothing. These install an ACTIVE
// threshold policy that *would* delete most of the product's support and assert
// the product still leaves every term in place until the caller asks.
// ---------------------------------------------------------------------------

/// A complex sum carrying a real, non-trivial coefficient floor.
type NewThreshSum = ppvm_pauli_sum_2::Sum<
    ppvm_pauli_sum_2::HashMapStore<NewWord, Complex<f64>>,
    ppvm_pauli_sum_2::CoefficientThreshold,
>;

fn build_thresh(n: usize, terms: &[(String, Complex<f64>)], tau: f64) -> NewThreshSum {
    NewThreshSum::from_terms_with_policy(
        n,
        ppvm_pauli_sum_2::CoefficientThreshold { threshold: tau },
        terms.iter().map(|(w, c)| (NewWord::from(w.as_str()), *c)),
    )
}

/// `multiply` / `multiply_into` / `mul_word_assign` must **not** invoke the
/// policy, even when the policy would delete nearly everything they produce —
/// `Sum::truncate` stays the single place `Policy::truncate` runs (behavioural
/// contract 1, which old upholds because no gate in `ppvm-pauli-sum` ever calls
/// `Strategy::truncate`).
///
/// The operands' coefficients are ~1e-4 — comfortably *above* the 1e-6 floor, so
/// the operands themselves are not at risk — while every product coefficient is
/// ~1e-8, two orders of magnitude *below* it. If any product path truncated, the
/// result would come back empty.
#[test]
fn the_product_never_truncates_even_under_an_active_policy() {
    let n = 3;
    let tau = 1e-6;
    let a_terms: Vec<(String, Complex<f64>)> = ["XYZ", "ZIX", "IZI", "XXI"]
        .iter()
        .map(|s| (s.to_string(), Complex::new(1e-4, 0.0)))
        .collect();
    let b_terms: Vec<(String, Complex<f64>)> = ["YZX", "IIZ"]
        .iter()
        .map(|s| (s.to_string(), Complex::new(1e-4, 0.0)))
        .collect();

    let a = build_thresh(n, &a_terms, tau);
    let b = build_thresh(n, &b_terms, tau);

    // Sanity, both directions: the floor leaves the *operands* alone (so an empty
    // result could only come from truncating the product), and it really would
    // delete terms at the product's magnitude if it ever ran.
    {
        let mut operands = a.clone();
        operands.truncate();
        assert_eq!(operands.len(), a.len(), "the floor must spare the operands");

        let mut probe = build_thresh(n, &[("III".to_string(), Complex::new(1e-8, 0.0))], tau);
        probe.truncate();
        assert_eq!(
            probe.len(),
            0,
            "the policy must be able to delete a product-magnitude term"
        );
    }

    // (i) allocating product
    let p = a.multiply(&b);
    assert_eq!(p.len(), a.len() * b.len(), "A·B kept every produced term");
    assert!(
        p.iter().all(|(_, c)| c.norm() < tau),
        "every product coefficient is below the floor, so a truncating product \
         would have emptied the support"
    );

    // (ii) accumulate-into form
    let mut acc =
        NewThreshSum::with_policy(n, ppvm_pauli_sum_2::CoefficientThreshold { threshold: tau });
    a.multiply_into(&b, &mut acc);
    assert_eq!(
        acc.len(),
        a.len() * b.len(),
        "multiply_into kept every term"
    );

    // (iii) in-place product
    let mut ip = a.clone();
    ip *= &b;
    assert_eq!(ip.len(), a.len() * b.len(), "A *= B kept every term");

    // (iv) single-word re-key
    let mut w = a.clone();
    w.mul_word_assign(&NewWord::from("YZX"));
    assert_eq!(w.len(), a.len(), "mul_word_assign kept every term");

    // …and the caller CAN still truncate afterwards — the policy is installed and
    // live, it simply never fires on its own.
    let mut after = p.clone();
    after.truncate();
    assert_eq!(
        after.len(),
        0,
        "explicit truncate() still applies the policy"
    );
}
