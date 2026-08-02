// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Lean-oracle property tests for the **L4 operator product** — checked against a
//! *dense `2ⁿ × 2ⁿ` matrix* model, not against the engine itself.
//!
//! `tests/pauli_sum_multiply_diff.rs` diffs the product against the old crate's
//! believed-correct single-word path. That is a differential check, but it shares
//! a phase convention with the thing under test: if `key_mul`'s `iᵏ` and the old
//! `PhasedPauliWord` product were *both* wrong the same way, both would agree.
//!
//! This file closes that by re-deriving the oracle from first principles. The
//! governing Lean spec is `twistedConv` in `lean/PPVM/Algebra/Twisted.lean`:
//!
//! ```text
//! twistedConv A B = A.sum fun p a => B.sum fun q b =>
//!     single (mulWord p q) (a * b * iPow i (phaseExpN p q))
//! ```
//!
//! and `lean/PPVM/Pauli/Matrix.lean` interprets a word as a genuine matrix via
//! `toOperator` (with `PPVM.PauliMatrix.trace_toOperator_mul` for the trace
//! statement). Composing the two, the property that pins **both** the key product
//! and its phase is the ring homomorphism
//!
//! ```text
//! mat(A · B) == mat(A) * mat(B)
//! ```
//!
//! where `mat(Σ cₖ k) = Σ cₖ · (Pauli Kronecker product of k)`. `mat` is injective
//! (the Paulis are an orthogonal basis of the matrix algebra), so equality of the
//! matrices is equality of the sums — and nothing about the engine's own phase
//! encoding enters the right-hand side. The Kronecker *ordering* is irrelevant as
//! long as it is used consistently: a different factor order is a relabelling of
//! the tensor factors, an algebra automorphism.
//!
//! Coverage per the component brief: **exhaustive** for the finite single-qubit
//! case (every one of the 16×16 support pairs, and every one of the 64 monomial
//! associativity triples), **randomized** for `n` qubits.

use std::collections::BTreeMap;

use num::Complex;
use ppvm_conformance_2::{PAULIS, random_pauli_string, seeded_rng};
use ppvm_pauli_sum_2::{PauliSum as NewPauliSum, PauliWord as NewWord};
use rand::RngExt;
use rand::rngs::StdRng;

type CSum = NewPauliSum<Complex<f64>>;
type C = Complex<f64>;

const TOL: f64 = 1e-12;
const SEEDS: [u64; 4] = [1, 7, 2024, 31337];

fn c(re: f64, im: f64) -> C {
    Complex::new(re, im)
}

fn zero() -> C {
    Complex::new(0.0, 0.0)
}

// ---------------------------------------------------------------------------
// The dense matrix model (the oracle) — `PPVM.PauliMatrix.toOperator`.
// ---------------------------------------------------------------------------

/// A dense complex matrix, row-major.
#[derive(Clone, PartialEq, Debug)]
struct Mat {
    dim: usize,
    a: Vec<C>,
}

impl Mat {
    fn zeros(dim: usize) -> Self {
        Self {
            dim,
            a: vec![zero(); dim * dim],
        }
    }

    fn at(&self, r: usize, col: usize) -> C {
        self.a[r * self.dim + col]
    }

    /// The Kronecker product `self ⊗ rhs`.
    fn kron(&self, rhs: &Mat) -> Mat {
        let dim = self.dim * rhs.dim;
        let mut out = Mat::zeros(dim);
        for i in 0..self.dim {
            for j in 0..self.dim {
                for k in 0..rhs.dim {
                    for l in 0..rhs.dim {
                        out.a[(i * rhs.dim + k) * dim + (j * rhs.dim + l)] =
                            self.at(i, j) * rhs.at(k, l);
                    }
                }
            }
        }
        out
    }

    fn matmul(&self, rhs: &Mat) -> Mat {
        assert_eq!(self.dim, rhs.dim);
        let mut out = Mat::zeros(self.dim);
        for i in 0..self.dim {
            for k in 0..self.dim {
                let aik = self.at(i, k);
                if aik == zero() {
                    continue;
                }
                for j in 0..self.dim {
                    out.a[i * self.dim + j] += aik * rhs.at(k, j);
                }
            }
        }
        out
    }

    fn add_scaled(&mut self, s: C, rhs: &Mat) {
        assert_eq!(self.dim, rhs.dim);
        for (o, r) in self.a.iter_mut().zip(rhs.a.iter()) {
            *o += s * r;
        }
    }

    fn max_abs_diff(&self, rhs: &Mat) -> f64 {
        assert_eq!(self.dim, rhs.dim);
        self.a
            .iter()
            .zip(rhs.a.iter())
            .map(|(x, y)| (x - y).norm())
            .fold(0.0, f64::max)
    }
}

/// The four single-qubit Paulis as `2 × 2` matrices (`PPVM.PauliMatrix`).
fn pauli_matrix(p: char) -> Mat {
    let m = match p {
        'I' => [c(1.0, 0.0), zero(), zero(), c(1.0, 0.0)],
        'X' => [zero(), c(1.0, 0.0), c(1.0, 0.0), zero()],
        // Y = [[0, -i], [i, 0]].
        'Y' => [zero(), c(0.0, -1.0), c(0.0, 1.0), zero()],
        'Z' => [c(1.0, 0.0), zero(), zero(), c(-1.0, 0.0)],
        other => panic!("not a Pauli letter: {other}"),
    };
    Mat {
        dim: 2,
        a: m.to_vec(),
    }
}

/// `toOperator` of a whole word: the Kronecker product over its sites.
fn word_matrix(s: &str) -> Mat {
    let mut it = s.chars();
    let mut m = pauli_matrix(it.next().expect("empty word"));
    for p in it {
        m = m.kron(&pauli_matrix(p));
    }
    m
}

/// `mat(Σ cₖ k)` — the dense image of a whole sum.
fn sum_matrix(s: &CSum) -> Mat {
    let mut out = Mat::zeros(1 << s.n_sites());
    for (k, v) in s.iter() {
        out.add_scaled(v, &word_matrix(&k.to_string()));
    }
    out
}

// ---------------------------------------------------------------------------
// Harness
// ---------------------------------------------------------------------------

fn build(n: usize, terms: &[(String, C)]) -> CSum {
    CSum::from_terms(
        n,
        terms.iter().map(|(w, v)| (NewWord::from(w.as_str()), *v)),
    )
}

/// `mat(A · B) == mat(A) · mat(B)` — the homomorphism property that *is*
/// `twistedConv` (`lean/PPVM/Algebra/Twisted.lean`) read through
/// `PPVM.PauliMatrix.toOperator`.
#[track_caller]
fn assert_product_matches_matrices(a: &CSum, b: &CSum, what: &str) {
    let product = a.multiply(b);
    let lhs = sum_matrix(&product);
    let rhs = sum_matrix(a).matmul(&sum_matrix(b));
    let d = lhs.max_abs_diff(&rhs);
    assert!(d <= TOL, "{what}: mat(A·B) != mat(A)·mat(B), max |Δ| = {d}");
}

fn support(s: &CSum) -> BTreeMap<String, C> {
    s.iter().map(|(k, v)| (k.to_string(), v)).collect()
}

#[track_caller]
fn assert_same_support(a: &CSum, b: &CSum, what: &str) {
    let (x, y) = (support(a), support(b));
    assert_eq!(
        x.keys().collect::<Vec<_>>(),
        y.keys().collect::<Vec<_>>(),
        "{what}: key sets differ"
    );
    for (k, v) in &x {
        let w = y[k];
        assert!((v - w).norm() <= TOL, "{what}: {k} -> {v} vs {w}");
    }
}

fn random_terms(rng: &mut StdRng, n: usize, count: usize) -> Vec<(String, C)> {
    (0..count)
        .map(|_| {
            let w = random_pauli_string(rng, n);
            (
                w,
                Complex::new(
                    rng.random_range(-2.0..2.0f64),
                    rng.random_range(-2.0..2.0f64),
                ),
            )
        })
        .collect()
}

// ---------------------------------------------------------------------------
// 1. EXHAUSTIVE single-qubit coverage
// ---------------------------------------------------------------------------

/// Every one of the 16 monomial pairs `(P, Q)` reproduces the Pauli
/// multiplication table exactly — key **and** `iᵏ` phase (`tmul`, the monomial
/// case of `twistedConv`).
#[test]
fn every_single_qubit_monomial_pair_matches_the_matrix_product() {
    for p in PAULIS {
        for q in PAULIS {
            let a = build(1, &[(p.to_string(), c(1.5, -0.5))]);
            let b = build(1, &[(q.to_string(), c(-0.25, 2.0))]);
            assert_product_matches_matrices(&a, &b, &format!("{p}·{q}"));
            // A monomial product is a single monomial: exactly one key out.
            assert_eq!(a.multiply(&b).len(), 1, "{p}·{q} support size");
        }
    }
}

/// **Exhaustive over the finite single-qubit case**: all 2⁴ = 16 possible
/// supports for `A` crossed with all 16 for `B` (256 pairs, including the empty
/// sum), each checked against the dense matrix product.
#[test]
fn every_single_qubit_support_pair_matches_the_matrix_product() {
    let coeff = |i: usize| c(0.5 + i as f64, -1.25 + 0.75 * i as f64);
    for mask_a in 0..16usize {
        for mask_b in 0..16usize {
            let ta: Vec<(String, C)> = (0..4)
                .filter(|i| mask_a >> i & 1 == 1)
                .map(|i| (PAULIS[i].to_string(), coeff(i)))
                .collect();
            let tb: Vec<(String, C)> = (0..4)
                .filter(|i| mask_b >> i & 1 == 1)
                .map(|i| (PAULIS[i].to_string(), coeff(3 - i)))
                .collect();
            let a = build(1, &ta);
            let b = build(1, &tb);
            assert_product_matches_matrices(&a, &b, &format!("mask {mask_a:04b}·{mask_b:04b}"));
        }
    }
}

/// **Exhaustive associativity** over all 4³ = 64 single-qubit monomial triples —
/// the Rust witness for `tmul_assoc` / `gtmul_assoc`
/// (`lean/PPVM/Algebra/Twisted.lean`), whose hypothesis `i⁴ = 1` holds on
/// `Complex<f64>`. The whole-map lift these monomial laws do *not* by themselves
/// give is `twistedConv_assoc` in the same file (it needs biadditivity, below).
#[test]
fn every_single_qubit_monomial_triple_is_associative() {
    for p in PAULIS {
        for q in PAULIS {
            for r in PAULIS {
                let a = build(1, &[(p.to_string(), c(1.0, 0.5))]);
                let b = build(1, &[(q.to_string(), c(-2.0, 0.25))]);
                let d = build(1, &[(r.to_string(), c(0.75, -1.5))]);
                let mut left = a.multiply(&b).multiply(&d);
                let mut right = a.multiply(&b.multiply(&d));
                left.reduce();
                right.reduce();
                assert_same_support(&left, &right, &format!("({p}·{q})·{r}"));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 2. Randomized n-qubit coverage against the same oracle
// ---------------------------------------------------------------------------

/// `mat(A · B) == mat(A)·mat(B)` on random `n`-qubit sums, `n = 1..4`
/// (`2ⁿ ≤ 16`, so the dense oracle stays cheap).
#[test]
fn random_n_qubit_products_match_the_matrix_product() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 1..=4usize {
            for (na, nb) in [(1usize, 1usize), (3, 4), (6, 5), (8, 1)] {
                let a = build(n, &random_terms(&mut rng, n, na));
                let b = build(n, &random_terms(&mut rng, n, nb));
                assert_product_matches_matrices(&a, &b, &format!("seed {seed}, n {n}"));
            }
        }
    }
}

/// The `iᵏ` twist is **not** a global convention that a phase-blind test would
/// miss: `X·Y = +iZ` while `Y·X = −iZ`, so the product is non-commutative in
/// exactly the way the `phaseExp` 2-cocycle says.
#[test]
fn the_twist_is_antisymmetric_on_anticommuting_words() {
    let x = build(1, &[("X".to_string(), c(1.0, 0.0))]);
    let y = build(1, &[("Y".to_string(), c(1.0, 0.0))]);
    let xy = x.multiply(&y);
    let yx = y.multiply(&x);
    let z = NewWord::from("Z");
    assert_eq!(xy.get(&z), Some(c(0.0, 1.0)), "X·Y must be +iZ");
    assert_eq!(yx.get(&z), Some(c(0.0, -1.0)), "Y·X must be −iZ");
}

// ---------------------------------------------------------------------------
// 3. Biadditivity — the property old's `MulAssign<PauliSum>` loses
// ---------------------------------------------------------------------------

/// Biadditivity in the **left** argument: `(A + B)·C == A·C + B·C`.
///
/// `tests/pauli_sum_multiply_diff.rs` covers the right argument (the side old's
/// per-rhs-term `map_add` fold destroys); this is the other half of what makes
/// `twistedConv` bilinear — machine-checked as `twistedConv_add_left` /
/// `twistedConv_add_right` in `lean/PPVM/Algebra/Twisted.lean`, which is also
/// the step lifting the monomial `tmul_assoc` to `twistedConv_assoc`.
/// `A + B` is built by `from_terms`, whose
/// `accumulate_batch` genuinely *sums* colliding keys (unlike old's
/// `AddAssign<PauliSum>`, which routes through `Extend` and hence overwrites —
/// `ppvm-pauli-sum/src/sum/ops.rs`).
#[test]
fn product_is_biadditive_in_the_left_argument() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 1..=4usize {
            let ta = random_terms(&mut rng, n, 5);
            let tb = random_terms(&mut rng, n, 4);
            let tc = random_terms(&mut rng, n, 3);

            let a = build(n, &ta);
            let b = build(n, &tb);
            let cc = build(n, &tc);
            let a_plus_b = build(n, &[ta.clone(), tb.clone()].concat());

            let mut lhs = a_plus_b.multiply(&cc);
            let mut rhs = CSum::new(n);
            a.multiply_into(&cc, &mut rhs);
            b.multiply_into(&cc, &mut rhs);
            lhs.reduce();
            rhs.reduce();
            assert_same_support(&lhs, &rhs, &format!("(A+B)·C (seed {seed}, n {n})"));
        }
    }
}

/// Scalar homogeneity in each argument: `(λA)·B == λ(A·B) == A·(λB)`, over a
/// genuinely complex `λ` so a mishandled `i` would show.
#[test]
fn product_is_homogeneous_in_each_argument() {
    let lambda = c(0.75, -1.25);
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 1..=3usize {
            let a = build(n, &random_terms(&mut rng, n, 5));
            let b = build(n, &random_terms(&mut rng, n, 4));

            let mut scaled_a = a.clone();
            scaled_a.scale(&lambda);
            let mut scaled_b = b.clone();
            scaled_b.scale(&lambda);
            let mut scaled_product = a.multiply(&b);
            scaled_product.scale(&lambda);

            assert_same_support(&scaled_a.multiply(&b), &scaled_product, "(λA)·B");
            assert_same_support(&a.multiply(&scaled_b), &scaled_product, "A·(λB)");
        }
    }
}

// ---------------------------------------------------------------------------
// 4. Finite support, and the support bound
// ---------------------------------------------------------------------------

/// `twistedConv` of two finitely-supported maps is finitely supported, with
/// `|supp(A·B)| ≤ |supp A| · |supp B|` — and, because `p ↦ p·q` is a bijection for
/// fixed `q`, at least `max(|A|, |B|)` when the other operand is a single
/// monomial. The bound is what licenses the store's `aux.reserve(max(|A|,|B|))`.
#[test]
fn product_support_is_finite_and_bounded() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 2..=4usize {
            for (na, nb) in [(1usize, 1usize), (4, 6), (9, 3)] {
                let a = build(n, &random_terms(&mut rng, n, na));
                let b = build(n, &random_terms(&mut rng, n, nb));
                let p = a.multiply(&b);
                assert!(
                    p.len() <= a.len() * b.len(),
                    "|A·B| = {} exceeds |A|·|B| = {}",
                    p.len(),
                    a.len() * b.len()
                );
                assert!(p.len() <= 1usize << (2 * n), "|A·B| exceeds the key space");
                // Single-monomial rhs: a bijection, so the support is preserved
                // exactly (this is the `mul_word_assign` fast path's invariant).
                if b.len() == 1 {
                    assert_eq!(p.len(), a.len(), "monomial rhs must preserve |support|");
                }
            }
        }
    }
}

/// Every product key is the key product of *some* pair, and every pair's key
/// product appears — the support identity `supp(A·B) ⊆ {p·q}` from the
/// definition of `twistedConv`.
#[test]
fn every_product_key_is_a_key_product_of_the_operands() {
    let mut rng = seeded_rng(4242);
    for n in 2..=4usize {
        let ta = random_terms(&mut rng, n, 5);
        let tb = random_terms(&mut rng, n, 4);
        let a = build(n, &ta);
        let b = build(n, &tb);
        let product = a.multiply(&b);

        // Recompute the expected key set directly from the definition, using
        // XOR on the x/z bit-planes (independent of the engine's `key_mul`).
        let mut expected: BTreeMap<String, ()> = BTreeMap::new();
        for (pw, _) in a.iter() {
            for (qw, _) in b.iter() {
                expected.insert(xor_words(&pw.to_string(), &qw.to_string()), ());
            }
        }
        for (k, _) in product.iter() {
            assert!(
                expected.contains_key(&k.to_string()),
                "product key {k} is not any p·q"
            );
        }
        assert_eq!(
            product.len(),
            expected.len(),
            "every p·q must appear (no key silently dropped)"
        );
    }
}

/// The Pauli group product on *letters* (`p·q` up to phase), computed from the
/// `(x, z)` bit-plane XOR — deliberately not routed through the engine's
/// `key_mul`, so it is an independent witness of the key half of the product.
fn xor_words(p: &str, q: &str) -> String {
    fn xz(ch: char) -> (u8, u8) {
        match ch {
            'I' => (0, 0),
            'X' => (1, 0),
            'Y' => (1, 1),
            'Z' => (0, 1),
            other => panic!("not a Pauli letter: {other}"),
        }
    }
    fn letter(x: u8, z: u8) -> char {
        match (x, z) {
            (0, 0) => 'I',
            (1, 0) => 'X',
            (1, 1) => 'Y',
            _ => 'Z',
        }
    }
    p.chars()
        .zip(q.chars())
        .map(|(a, b)| {
            let ((ax, az), (bx, bz)) = (xz(a), xz(b));
            letter(ax ^ bx, az ^ bz)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// 5. `twistedConv_apply_id` — the L3↔L4 tie, against the matrix trace
// ---------------------------------------------------------------------------

/// `(A·B)[I] == ⟨A, B⟩ == Tr(mat(A)·mat(B)) / 2ⁿ`.
///
/// `pauli_sum_multiply_diff.rs` checks the first equality against the engine's own
/// `overlap`; this adds the second, against the honest matrix trace
/// (`PPVM.PauliMatrix.trace_toOperator_mul`), so the identity coefficient is
/// pinned to physics and not just to internal consistency.
#[test]
fn identity_coefficient_is_the_normalized_matrix_trace() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 1..=4usize {
            let a = build(n, &random_terms(&mut rng, n, 6));
            let b = build(n, &random_terms(&mut rng, n, 6));
            let product = a.multiply(&b);
            let id = NewWord::from("I".repeat(n).as_str());
            let got = product.get(&id).unwrap_or(zero());

            let m = sum_matrix(&a).matmul(&sum_matrix(&b));
            let dim = 1usize << n;
            let mut tr = zero();
            for i in 0..dim {
                tr += m.at(i, i);
            }
            let expected = tr / (dim as f64);
            assert!(
                (got - expected).norm() <= TOL,
                "(A·B)[I] {got} != Tr(AB)/2ⁿ {expected} (seed {seed}, n {n})"
            );
            // And the engine's own L3 pairing agrees with both.
            assert!((a.overlap(&b) - expected).norm() <= TOL, "overlap vs trace");
        }
    }
}

// ---------------------------------------------------------------------------
// 6. The in-place form is the same algebra
// ---------------------------------------------------------------------------

/// `A *= B` (which draws its accumulator from the store's persistent aux
/// double-buffer) computes exactly the allocating `A.multiply(&B)`, on random
/// operands — the store's buffer reuse must not change the algebra, and in
/// particular must not leak stale entries from a previous product.
#[test]
fn in_place_product_matches_the_allocating_form_repeatedly() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 1..=4usize {
            let a = build(n, &random_terms(&mut rng, n, 5));
            let b = build(n, &random_terms(&mut rng, n, 4));
            let d = build(n, &random_terms(&mut rng, n, 3));

            // Two products in a row through the SAME store: the second must not
            // see anything left in `aux` by the first.
            let mut got = a.clone();
            got *= &b;
            assert_same_support(&got, &a.multiply(&b), &format!("A*=B (seed {seed}, n {n})"));
            got *= &d;
            let expected = a.multiply(&b).multiply(&d);
            assert_same_support(&got, &expected, &format!("(A*=B)*=D (seed {seed}, n {n})"));
        }
    }
}

/// The in-place product against the matrix oracle too — buffer reuse included.
#[test]
fn in_place_product_matches_the_matrix_product() {
    let mut rng = seeded_rng(90210);
    for n in 1..=3usize {
        let a = build(n, &random_terms(&mut rng, n, 5));
        let b = build(n, &random_terms(&mut rng, n, 4));
        let expected = sum_matrix(&a).matmul(&sum_matrix(&b));
        let mut got = a.clone();
        got *= &b;
        let d = sum_matrix(&got).max_abs_diff(&expected);
        assert!(d <= TOL, "mat(A *= B) != mat(A)·mat(B), max |Δ| = {d}");
    }
}

// ---------------------------------------------------------------------------
// Store-buffer invariant: the product shares `aux` with the Clifford re-key and
// `scratch` with the rotations. Nothing above ever runs a product and a GATE
// through the same store, so a product that left `aux` non-empty (or a gate that
// left it dirty for the next product) would go unnoticed.
// ---------------------------------------------------------------------------

/// One step of the interleaved operation sequence below.
type Step = Box<dyn Fn(&mut CSum)>;

/// A product and the gate kernels draw on the **same** `HashMapStore` buffers:
/// `multiply_in_place`/`mul_word_assign` use `aux`, as does the Clifford re-key,
/// while the rotations stage branch terms in `scratch`. Each operation both
/// clears the buffer it uses and is documented to leave it empty afterwards
/// (`crates/ppvm-pauli-sum-2/src/store.rs`), so *any* interleaving must be
/// invisible.
///
/// Oracle: replay the identical operation sequence, but move the state into a
/// **freshly allocated** `Sum` after every step (a clone starts `aux`/`scratch`
/// empty by construction). If a buffer leaked state across operation *kinds*, the
/// interleaved run would diverge from the fresh-per-step run.
#[test]
fn products_and_gates_interleave_without_leaking_store_buffers() {
    use ppvm_traits_2::{Clifford, RotationOne};

    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 2..=4usize {
            let a = build(n, &random_terms(&mut rng, n, 6));
            let b = build(n, &random_terms(&mut rng, n, 4));
            let q = NewWord::from(random_pauli_string(&mut rng, n).as_str());

            // Every operation kind, interleaved so each buffer is used by two
            // different code paths back to back:
            //   rekey (aux) → product (aux) → rotation (scratch) → word-mul (aux)
            //   → rotation (scratch) → product (aux) → rekey (aux)
            let ops: Vec<Step> = vec![
                Box::new(|s: &mut CSum| s.cnot(0, 1)),
                Box::new({
                    let b = b.clone();
                    move |s: &mut CSum| *s *= &b
                }),
                Box::new(|s: &mut CSum| s.rz(0, 0.37)),
                Box::new({
                    let q = q.clone();
                    move |s: &mut CSum| s.mul_word_assign(&q)
                }),
                Box::new(|s: &mut CSum| s.rx(1, 0.21)),
                Box::new({
                    let b = b.clone();
                    move |s: &mut CSum| *s *= &b
                }),
                Box::new(|s: &mut CSum| s.h(0)),
            ];

            // Run 1: everything through ONE store, so buffers are reused.
            let mut shared = a.clone();
            for op in &ops {
                op(&mut shared);
            }

            // Run 2: after each step, continue in a fresh store.
            let mut fresh = a.clone();
            for op in &ops {
                op(&mut fresh);
                // `from_terms` allocates a new store with empty aux/scratch.
                fresh = CSum::from_terms(n, fresh.iter().collect::<Vec<_>>());
            }

            // Neither run canonicalizes: `from_terms` deliberately does NOT
            // `reduce` (it preserves exact-zero terms, as old does), so the two
            // supports must match *exactly* — zero-coefficient terms included.
            assert_same_support(
                &shared,
                &fresh,
                &format!("interleaved products+gates (seed {seed}, n {n})"),
            );
        }
    }
}
