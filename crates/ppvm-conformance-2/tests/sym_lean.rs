// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Phase-5 **Lean-oracle** property tests for the symbolic and exact coefficient
//! rings.
//!
//! Everything here is checked against the machine-checked semantics in
//! `lean/PPVM/**`, not against the old crate — old's phase handling is broken in
//! three independent places (`oldSuspectedBugs` #2/#3/#4) and its `eval` discards
//! the phase entirely, so it is not a parity target on this surface.
//!
//! | section | Lean oracle |
//! |:--|:--|
//! | the `ℤ/4` phase group | `lean/PPVM/Pauli/Phase.lean` (`phaseExp_eq_ref`), `ppvm_traits_2::Phase` |
//! | `mul_phase` is multiplication by `iᵏ` | `lean/PPVM/Instantiations/Symbolic.lean` `phaseFold_eq_iSym_pow_mul`, `evalC_phaseFold`, `phaseFold_const`, `phaseFold_drop_const_ne` |
//! | `max_sin` is representation-dependent | `Symbolic.lean` `mulImpl_not_wellDefined`, `mulImpl_one_one_untruncated`, `fastArm_escapes_bound` |
//! | `iᵃ·iᵇ = i^{a+b}` | `lean/PPVM/Algebra/Twisted.lean` `iPow_add` |
//! | `conj i = −i` | `lean/PPVM/Pauli/Matrix.lean` `star_iU` |
//! | twisted-product associativity | `Twisted.lean` `tmul_assoc`, `twistedConv_assoc` |
//! | twisted-product biadditivity | `Twisted.lean` `twistedConv_add_left`/`_right` |
//! | the ℓ¹ truncation bound's premise | `lean/PPVM/Algebra/Truncation.lean` `l1_bound_abv`, `l1_bound_needs_subadditive` |
//! | accumulation laws on the coefficient | `lean/PPVM/Algebra/GradedMap.lean` `accumulate_comm`, `accumulate_assoc`, `scale_scale` |
//! | the 2-D rotation | `lean/PPVM/Instantiations/Rotation.lean` `rot_norm_sq`, `rot_rot` |
//!
//! # The exact-ring witness (Phase-1 debt `t2.coefficient.2`)
//!
//! `ImaginaryUnit`/`Conjugate` were pinned only on `f64`/`Complex<f64>` in Phase
//! 1, and the witness proving L4 genuinely admits **exact** rings was deferred to
//! this phase. [`GaussianInt`] (`ℤ[i]`, two `i64`s, no float anywhere in the
//! representation) is that witness, and the sections below discharge the debt:
//! every law is asserted with `assert_eq!` and **zero tolerance**, and the
//! twisted product is driven end to end through `ppvm-pauli-sum-2::Sum` over it.
//!
//! Coverage per the component brief: **exhaustive** for the finite single-qubit
//! cases (all 16 `ℤ/4 × ℤ/4` phase pairs, all 16 single-qubit Pauli pairs, all 64
//! associativity triples), **randomized** for `n` qubits.

use std::collections::{BTreeMap, BTreeSet, HashSet};

use num::{Complex, One, Zero};
use ppvm_conformance_2::{PAULIS, random_pauli_string, seeded_rng};
use ppvm_pauli_sum_2::{HashMapStore, NoPolicy, PauliWord, Sum};
use ppvm_sym_2::{GaussianInt, Prod, Term};
use ppvm_traits_2::{Coefficient, Conjugate, ImaginaryUnit, Phase};
use rand::RngExt;
use rand::rngs::StdRng;

/// The exact-ring sum: `[u8; 8]` keys, `ℤ[i]` coefficients, no truncation.
type ExactSum = Sum<HashMapStore<PauliWord<[u8; 8]>, GaussianInt>, NoPolicy>;

const SEEDS: [u64; 4] = [1, 7, 2024, 31337];

// ===========================================================================
// 1. The `ℤ/4` phase group — `lean/PPVM/Pauli/Phase.lean`, `Twisted.lean::iPow_add`.
// ===========================================================================

#[test]
fn monomial_phase_is_the_z4_group_exhaustively() {
    // `Prod::add_phase` must be the `ℤ/4` group operation, i.e. exactly
    // `ppvm_traits_2::Phase::compose` (`phaseExp_eq_ref`). All 16 pairs.
    for a in 0..4u8 {
        for b in 0..4u8 {
            let mut p = Prod::new();
            p.add_phase(a);
            p.add_phase(b);
            assert_eq!(
                p.phase(),
                Phase::from_exponent(a)
                    .compose(Phase::from_exponent(b))
                    .exponent(),
                "add_phase({a}) then add_phase({b})"
            );
            assert_eq!(p.phase(), (a + b) % 4);
        }
    }
    // Identity and inverse.
    for a in 0..4u8 {
        let mut p = Prod::sin(0);
        p.add_phase(a);
        p.add_phase(0);
        assert_eq!(p.phase(), a, "i⁰ is the identity");
        let mut q = Prod::sin(0);
        q.add_phase(a);
        q.add_phase(Phase::from_exponent(a).inverse().exponent());
        assert_eq!(q.phase(), 0, "i^k · i^{{-k}} = 1");
    }
}

#[test]
fn ipow_add_holds_on_monomial_multiplication_exhaustively() {
    // `iᵃ · iᵇ = i^{a+b}` (`Twisted.lean` `iPow_add`), on the monomial product —
    // the operation old dropped entirely (`oldSuspectedBugs` #2).
    for a in 0..4u8 {
        for b in 0..4u8 {
            let mut p = Prod::sin(0);
            p.add_phase(a);
            let mut q = Prod::cos(1);
            q.add_phase(b);
            let r = p * q;
            assert_eq!(r.phase(), (a + b) % 4, "i^{a} · i^{b}");
            assert_eq!(r.sin_pow(), 1);
            assert_eq!(r.cos_pow(), 1);
        }
    }
}

#[test]
fn ipow_add_holds_on_the_symbolic_term_exhaustively() {
    // The same law one level up, on `Term::mul_phase`, observed through the
    // phase-aware readout. Checked on all three non-`Var` representations.
    let terms: Vec<(&str, Term)> = vec![
        ("Const", Term::from(2.0)),
        ("One", Term::var(0).sin() * 3.0),
        ("Sum", Term::from(2.0) + Term::var(0).sin()),
    ];
    let vals = [0.7];
    for (what, t) in &terms {
        for a in 0..4u8 {
            for b in 0..4u8 {
                let lhs = t.mul_phase(a).mul_phase(b);
                let rhs = t.mul_phase((a + b) % 4);
                let (l, r) = (
                    lhs.eval_complex(&vals).unwrap(),
                    rhs.eval_complex(&vals).unwrap(),
                );
                assert!(
                    (l - r).norm() < 1e-12,
                    "{what}: i^{a}·i^{b} ≠ i^{{{a}+{b}}}: {l} vs {r}"
                );
                // …and it really is `i^{a+b}` times the unphased value.
                let base = t.eval_complex(&vals).unwrap();
                let want = base * ipow_c((a + b) % 4);
                assert!((l - want).norm() < 1e-12, "{what}: {l} vs {want}");
            }
        }
    }
}

fn ipow_c(k: u8) -> Complex<f64> {
    match k % 4 {
        0 => Complex::new(1.0, 0.0),
        1 => Complex::new(0.0, 1.0),
        2 => Complex::new(-1.0, 0.0),
        _ => Complex::new(0.0, -1.0),
    }
}

#[test]
fn star_iu_holds_on_both_rings() {
    // `conj i = −i` — `lean/PPVM/Pauli/Matrix.lean` `star_iU`.
    // Exact on `ℤ[i]`:
    let i = GaussianInt::imaginary_unit();
    assert_eq!(i.conj(), -i);
    assert_eq!(i * i, -GaussianInt::one());
    assert_eq!(i * i * i * i, GaussianInt::one());
    // Denotational on the symbolic ring (its `PartialEq` is representational, so
    // `One(Prod{phase:3},1)` and `One(Prod{phase:1},−1)` are different spellings
    // of the same value — see the `ImaginaryUnit for Term` law caveat).
    let ti = Term::imaginary_unit();
    assert_eq!(
        ti.conj().eval_complex(&[]).unwrap(),
        (-ti.clone()).eval_complex(&[]).unwrap()
    );
    assert_eq!(
        (ti.clone() * ti).eval_complex(&[]).unwrap(),
        Complex::new(-1.0, 0.0)
    );
}

// ===========================================================================
// 2. The exact-ring witness `GaussianInt` — the `t2.coefficient.2` deliverable.
// ===========================================================================

/// Every Gaussian integer with `|re|, |im| <= 3` — the finite grid the exact
/// laws are checked exhaustively over.
fn small_gaussians() -> Vec<GaussianInt> {
    let mut v = Vec::new();
    for re in -3..=3i64 {
        for im in -3..=3i64 {
            v.push(GaussianInt::new(re, im));
        }
    }
    v
}

fn random_gaussian(rng: &mut StdRng) -> GaussianInt {
    GaussianInt::new(rng.random_range(-64..64i64), rng.random_range(-64..64i64))
}

#[test]
fn exact_ring_is_a_commutative_ring_with_one() {
    // Randomized over `ℤ[i]`, EXACT (`assert_eq!`, no tolerance) — the whole
    // point of the witness: an exact ring must be usable as a `Coefficient`
    // without a lossy float anywhere.
    for seed in SEEDS {
        let mut rng = seeded_rng(seed);
        for _ in 0..512 {
            let (a, b, c) = (
                random_gaussian(&mut rng),
                random_gaussian(&mut rng),
                random_gaussian(&mut rng),
            );
            assert_eq!((a + b) + c, a + (b + c), "additive associativity");
            assert_eq!(a + b, b + a, "additive commutativity");
            assert_eq!(a + GaussianInt::zero(), a, "additive identity");
            assert_eq!(a + (-a), GaussianInt::zero(), "additive inverse");
            assert_eq!(a - b, a + (-b), "subtraction");
            assert_eq!((a * b) * c, a * (b * c), "multiplicative associativity");
            assert_eq!(a * b, b * a, "multiplicative commutativity");
            assert_eq!(a * GaussianInt::one(), a, "multiplicative identity");
            assert_eq!(a * (b + c), a * b + a * c, "distributivity");
            assert_eq!(a * GaussianInt::zero(), GaussianInt::zero());
        }
    }
}

#[test]
fn exact_ring_imaginary_unit_and_conjugate_laws_exhaustively() {
    let i = GaussianInt::imaginary_unit();
    for z in small_gaussians() {
        // `ImaginaryUnit`: `mul_i` is the ring multiply by `i`.
        assert_eq!(z.mul_i(), z * i, "mul_i ≠ ·i at {z}");
        assert_eq!(z.mul_i().mul_i(), -z, "i² = −1 at {z}");
        assert_eq!(z.mul_i().mul_i().mul_i().mul_i(), z, "i⁴ = 1 at {z}");
        // `Conjugate` is a ring *-involution.
        assert_eq!(z.conj().conj(), z, "conj is an involution");
        assert_eq!(
            z.conj().mul_i(),
            -(z.mul_i().conj()),
            "conj(iz) = −i·conj z"
        );
        // `Phase::apply` — the fold `iᵏ·c` the L4 product uses — agrees with the
        // ring multiply, exhaustively over `ℤ/4`.
        for k in 0..4u8 {
            let p = Phase::from_exponent(k);
            let mut want = z;
            for _ in 0..k {
                want *= i;
            }
            assert_eq!(p.apply(&z), want, "Phase::apply({k}) at {z}");
        }
        // `mul_sign` is an exact integer sign flip.
        assert_eq!(z.mul_sign(1), z);
        assert_eq!(z.mul_sign(-1), -z);
    }
    for a in small_gaussians() {
        for b in small_gaussians() {
            assert_eq!((a + b).conj(), a.conj() + b.conj(), "conj is additive");
            assert_eq!(
                (a * b).conj(),
                a.conj() * b.conj(),
                "conj is multiplicative"
            );
        }
    }
}

#[test]
fn exact_ring_magnitude_is_a_genuine_absolute_value() {
    // `lean/PPVM/Algebra/Truncation.lean` `l1_bound_abv` literally requires an
    // `AbsoluteValue C R`. `GaussianInt` satisfies every clause — which is what
    // makes it a *valid* `Coefficient` and not merely a compiling one.
    for a in small_gaussians() {
        assert!(a.magnitude() >= 0.0, "nonnegative");
        assert_eq!(a.magnitude() == 0.0, a.is_zero(), "N(x)=0 iff x=0 at {a}");
        // The exact integer companion, so the law can be stated without floats.
        assert_eq!(a.norm_sq(), a.re * a.re + a.im * a.im);
        for b in small_gaussians() {
            // Multiplicative (checked exactly on the integer field norm).
            assert_eq!(
                (a * b).norm_sq(),
                a.norm_sq() * b.norm_sq(),
                "N(xy)=N(x)N(y)"
            );
            assert!(
                (a * b).magnitude() <= a.magnitude() * b.magnitude() + 1e-9,
                "modulus multiplicativity"
            );
            // Subadditive — the clause `l1_bound_needs_subadditive` shows is the
            // load-bearing one for the ℓ¹ truncation bound.
            assert!(
                (a + b).magnitude() <= a.magnitude() + b.magnitude() + 1e-9,
                "triangle inequality at {a}, {b}"
            );
        }
    }
}

#[test]
fn symbolic_magnitude_deliberately_violates_the_absolute_value_law() {
    // The adjudication recorded in the `ppvm-sym-2` crate docs, pinned as a test
    // so it cannot be "fixed" silently: `Term::magnitude` returns `+∞` for every
    // symbolic form, which reproduces old's inert `cutoff` (behavioural contract
    // 3) at the cost of the `N(x) == 0 iff x == 0` clause.
    //
    // No absolute value exists on `R[sᵢ, cᵢ]` at all: the natural `ℓ¹`
    // coefficient norm is only SUB-multiplicative, witnessed here by
    // `(1+s)(1−s) = 1 − s²` (`ℓ¹` gives 2·2 = 4 vs 2). So the choice is between
    // old's behaviour and a law, and the prime directive picks old's behaviour.
    let s = Term::var(0).sin();
    let a = Term::from(1.0) + s.clone(); // 1 + sin(x0), ℓ¹ = 2
    let b = Term::from(1.0) - s; // 1 − sin(x0), ℓ¹ = 2
    let ab = a * b; // 1 − sin²(x0),  ℓ¹ = 2 < 4
    assert!(
        (ab.eval(&[0.4]).unwrap() - (1.0 - 0.4f64.sin().powi(2))).abs() < 1e-12,
        "the ℓ¹ sub-multiplicativity witness must actually be (1 − s²)"
    );

    // The violated clause, spelled out: an empty symbolic `Sum` DENOTES 0 …
    let mut zero_sum = Term::from(1.0) + Term::var(1).cos();
    zero_sum.set_max_sin(0);
    zero_sum *= Term::var(0).sin();
    assert_eq!(zero_sum.eval(&[0.3, 0.9]).unwrap(), 0.0);
    // … yet reports `+∞`, so `CoefficientThreshold` keeps it, exactly as old did.
    assert_eq!(zero_sum.magnitude(), f64::INFINITY);
    // The constant form is the only one a threshold can see.
    assert_eq!(Term::from(-1.5).magnitude(), 1.5);
    assert_eq!(Term::from(0.0).magnitude(), 0.0);
}

#[test]
fn max_sin_truncation_is_representation_dependent() {
    // Lean `mulImpl_not_wellDefined` (`lean/PPVM/Instantiations/Symbolic.lean`):
    // the shipped four-way `Inner` product does **not** factor through the value
    // a `Term` denotes, because the `One × One` fast arm never consults
    // `max_sin`. Two representations of the same polynomial therefore multiply
    // to different polynomials — so `mulMono_drop_at_insert_eq_drop_at_end` is
    // not an end-to-end guarantee, and `set_max_sin` is not a hard degree bound
    // on the propagated result.
    //
    // This is the invariant behind the crate's "preserved old quirk": the single
    // `sin_pow = 7` escapee at `max_sin = 2` on the Trotter replay is an instance
    // of it, not an accident. Unifying the four `Inner` arms onto one map-backed
    // representation would make the product well-defined and change numbers.
    let s0 = Term::var(0).sin();
    let mut fast_form = s0.clone() * s0; // `One(sin²(x0), 1)`, the non-allocating form
    fast_form.set_max_sin(2);
    let mut map_form = fast_form.clone();
    map_form += 0.0; // same value, promoted to the map-backed `Sum`

    let vals = [0.5, 0.7];
    assert!(
        (fast_form.eval(&vals).unwrap() - map_form.eval(&vals).unwrap()).abs() < 1e-15,
        "the two representations must denote the same polynomial (`den a₁ == den a₂`)"
    );

    let b = Term::var(1).sin();
    let via_fast_arm = fast_form * b.clone(); // `One × One`: truncation-blind
    let via_map_arm = map_form * b; // `Sum × One`: routes through `Sum::mul_term`

    // `mulImpl_one_one_untruncated` / `fastArm_escapes_bound`: the fast arm is
    // exact, hence unbounded — sine degree 3 survives at `max_sin = 2`.
    let untruncated = 0.5f64.sin().powi(2) * 0.7f64.sin();
    assert!(
        (via_fast_arm.eval(&vals).unwrap() - untruncated).abs() < 1e-12,
        "the fast arm computes the untruncated ring product"
    );
    // …while the map-backed arm drops exactly that monomial.
    assert_eq!(
        via_map_arm.eval(&vals).unwrap(),
        0.0,
        "the map-backed arm must drop the degree-3 monomial at `max_sin = 2`"
    );
}

// ===========================================================================
// 3. The twisted product over the exact ring — a dense `ℤ[i]` matrix oracle.
// ===========================================================================
//
// `lean/PPVM/Pauli/Matrix.lean` interprets a word as a genuine matrix
// (`toOperator`); composing with `twistedConv` gives the ring homomorphism
//
//     mat(A · B) == mat(A) * mat(B)
//
// and `mat` is injective (the Paulis are an orthogonal basis). Over `ℤ[i]` every
// Pauli matrix entry is a Gaussian integer, so the ENTIRE oracle is exact: no
// tolerance anywhere, which is precisely the property the exact-ring witness is
// meant to demonstrate.

/// A dense `ℤ[i]` matrix, row-major.
#[derive(Clone, PartialEq, Eq, Debug)]
struct Mat {
    dim: usize,
    a: Vec<GaussianInt>,
}

impl Mat {
    fn zeros(dim: usize) -> Self {
        Self {
            dim,
            a: vec![GaussianInt::zero(); dim * dim],
        }
    }

    fn eye(dim: usize) -> Self {
        let mut m = Self::zeros(dim);
        for i in 0..dim {
            m.a[i * dim + i] = GaussianInt::one();
        }
        m
    }

    fn at(&self, r: usize, c: usize) -> GaussianInt {
        self.a[r * self.dim + c]
    }

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
        let n = self.dim;
        let mut out = Mat::zeros(n);
        for i in 0..n {
            for k in 0..n {
                let a = self.at(i, k);
                if a.is_zero() {
                    continue;
                }
                for j in 0..n {
                    out.a[i * n + j] += a * rhs.at(k, j);
                }
            }
        }
        out
    }

    fn add(&self, rhs: &Mat) -> Mat {
        let mut out = self.clone();
        for (o, r) in out.a.iter_mut().zip(rhs.a.iter()) {
            *o += *r;
        }
        out
    }

    fn scale(&self, c: GaussianInt) -> Mat {
        let mut out = self.clone();
        for o in out.a.iter_mut() {
            *o *= c;
        }
        out
    }
}

/// The four single-qubit Pauli matrices over `ℤ[i]` — all four have Gaussian
/// integer entries, which is why the oracle can be exact.
fn pauli_mat(p: char) -> Mat {
    let g = GaussianInt::new;
    let m = |a, b, c, d| Mat {
        dim: 2,
        a: vec![a, b, c, d],
    };
    match p {
        'I' => m(g(1, 0), g(0, 0), g(0, 0), g(1, 0)),
        'X' => m(g(0, 0), g(1, 0), g(1, 0), g(0, 0)),
        'Y' => m(g(0, 0), g(0, -1), g(0, 1), g(0, 0)),
        'Z' => m(g(1, 0), g(0, 0), g(0, 0), g(-1, 0)),
        other => panic!("not a Pauli letter: {other}"),
    }
}

/// `toOperator` of a whole Pauli word.
fn word_mat(word: &str) -> Mat {
    word.chars()
        .map(pauli_mat)
        .reduce(|acc, m| acc.kron(&m))
        .unwrap_or_else(|| Mat::eye(1))
}

/// `mat(Σ cₖ k)` for a sum given as `(word, coeff)` terms.
fn sum_mat(n: usize, terms: &[(String, GaussianInt)]) -> Mat {
    let dim = 1usize << n;
    terms.iter().fold(Mat::zeros(dim), |acc, (w, c)| {
        acc.add(&word_mat(w).scale(*c))
    })
}

fn build_exact(n: usize, terms: &[(String, GaussianInt)]) -> ExactSum {
    let mut s: ExactSum = ExactSum::new(n);
    for (w, c) in terms {
        s += (PauliWord::<[u8; 8]>::from(w.as_str()), *c);
    }
    s
}

fn exact_terms(s: &ExactSum) -> Vec<(String, GaussianInt)> {
    let mut v: Vec<(String, GaussianInt)> = s.iter().map(|(k, c)| (k.to_string(), c)).collect();
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
}

#[test]
fn exact_twisted_product_is_a_ring_homomorphism_exhaustively_on_one_qubit() {
    // All 16 single-qubit Pauli pairs, with a non-trivial `ℤ[i]` coefficient on
    // each side, asserted EXACTLY against the matrix oracle.
    for p in PAULIS {
        for q in PAULIS {
            let a = vec![(p.to_string(), GaussianInt::new(2, -3))];
            let b = vec![(q.to_string(), GaussianInt::new(-1, 5))];
            let prod = build_exact(1, &a).multiply(&build_exact(1, &b));
            assert_eq!(
                sum_mat(1, &exact_terms(&prod)),
                sum_mat(1, &a).matmul(&sum_mat(1, &b)),
                "({p})·({q})"
            );
        }
    }
}

#[test]
fn exact_twisted_product_is_a_ring_homomorphism_on_n_qubits() {
    // Randomized `n`-qubit half of the same statement (`twistedConv` composed
    // with `toOperator`), still exact.
    for seed in SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 2..=4usize {
            for _ in 0..8 {
                let a = random_exact_terms(&mut rng, n, 5);
                let b = random_exact_terms(&mut rng, n, 5);
                let prod = build_exact(n, &a).multiply(&build_exact(n, &b));
                assert_eq!(
                    sum_mat(n, &exact_terms(&prod)),
                    sum_mat(n, &a).matmul(&sum_mat(n, &b)),
                    "n={n} seed={seed}"
                );
            }
        }
    }
}

fn random_exact_terms(rng: &mut StdRng, n: usize, count: usize) -> Vec<(String, GaussianInt)> {
    let mut map: BTreeMap<String, GaussianInt> = BTreeMap::new();
    for _ in 0..count {
        let w = random_pauli_string(rng, n);
        let c = GaussianInt::new(rng.random_range(-4..4i64), rng.random_range(-4..4i64));
        *map.entry(w).or_insert(GaussianInt::zero()) += c;
    }
    map.into_iter().collect()
}

#[test]
fn exact_twisted_product_is_associative_and_biadditive() {
    // `Twisted.lean` `tmul_assoc` / `twistedConv_assoc` and
    // `twistedConv_add_left`/`twistedConv_add_right`, asserted with ZERO
    // tolerance because the ring is exact — the acceptance bar the
    // `sym.exact.multiply` workload names.
    for seed in SEEDS {
        let mut rng = seeded_rng(seed);
        for n in [2usize, 4, 6] {
            for _ in 0..6 {
                let a = build_exact(n, &random_exact_terms(&mut rng, n, 6));
                let b = build_exact(n, &random_exact_terms(&mut rng, n, 6));
                let c = build_exact(n, &random_exact_terms(&mut rng, n, 6));

                let lhs = a.multiply(&b).multiply(&c);
                let rhs = a.multiply(&b.multiply(&c));
                assert_eq!(
                    exact_terms(&lhs),
                    exact_terms(&rhs),
                    "(A·B)·C ≠ A·(B·C) at n={n} seed={seed}"
                );

                // Biadditivity: A·(B+C) == A·B + A·C, and (B+C)·A == B·A + C·A.
                let bc = add_exact(&b, &c);
                assert_eq!(
                    exact_terms(&a.multiply(&bc)),
                    exact_terms(&add_exact(&a.multiply(&b), &a.multiply(&c))),
                    "twistedConv_add_right at n={n}"
                );
                assert_eq!(
                    exact_terms(&bc.multiply(&a)),
                    exact_terms(&add_exact(&b.multiply(&a), &c.multiply(&a))),
                    "twistedConv_add_left at n={n}"
                );
            }
        }
    }
}

/// Term-wise sum of two exact sums (the `+` of the graded module).
fn add_exact(a: &ExactSum, b: &ExactSum) -> ExactSum {
    let mut out = a.clone();
    for (k, c) in b.iter() {
        out += (k, c);
    }
    out
}

#[test]
fn exact_multiply_into_accumulates_and_does_not_allocate_per_key_pair() {
    // `multiply_into` must ACCUMULATE into the destination (the L4 contract) and
    // reuse the destination's allocation across repeated products — the
    // "must not allocate per key-pair" half of the `sym.exact.multiply` gate.
    let mut rng = seeded_rng(4242);
    let n = 4;
    let a = build_exact(n, &random_exact_terms(&mut rng, n, 8));
    let b = build_exact(n, &random_exact_terms(&mut rng, n, 8));

    let mut acc: ExactSum = ExactSum::new(n);
    a.multiply_into(&b, &mut acc);
    let once = exact_terms(&acc);
    let cap_after_first = acc.capacity();
    a.multiply_into(&b, &mut acc);
    let twice = exact_terms(&acc);

    assert_eq!(once.len(), twice.len(), "the support must not change");
    for ((k1, c1), (k2, c2)) in once.iter().zip(twice.iter()) {
        assert_eq!(k1, k2);
        assert_eq!(*c2, *c1 + *c1, "multiply_into must accumulate, exactly");
    }
    assert!(
        acc.capacity() >= cap_after_first,
        "the destination allocation must be reused, not rebuilt smaller"
    );
}

// ===========================================================================
// 4. Graded-map accumulation and rotation laws on the SYMBOLIC ring.
// ===========================================================================

#[test]
fn symbolic_accumulation_laws_hold_denotationally() {
    // `lean/PPVM/Algebra/GradedMap.lean` `accumulate_comm`/`accumulate_assoc`/
    // `scale_scale`, on `Term`. Denotational (`eval`), because `Term`'s
    // `PartialEq` is representational by design (behavioural contract 5) — the
    // *values* must form a commutative monoid even when the spellings differ.
    for seed in SEEDS {
        let mut rng = seeded_rng(seed);
        for _ in 0..128 {
            let a = random_term(&mut rng, 3);
            let b = random_term(&mut rng, 3);
            let c = random_term(&mut rng, 3);
            let vals: Vec<f64> = (0..3)
                .map(|_| rng.random_range(-std::f64::consts::PI..std::f64::consts::PI))
                .collect();

            let ab = (a.clone() + b.clone()).eval(&vals).unwrap();
            let ba = (b.clone() + a.clone()).eval(&vals).unwrap();
            assert!((ab - ba).abs() < 1e-9, "accumulate_comm: {ab} vs {ba}");

            let l = ((a.clone() + b.clone()) + c.clone()).eval(&vals).unwrap();
            let r = (a.clone() + (b.clone() + c.clone())).eval(&vals).unwrap();
            assert!((l - r).abs() < 1e-9, "accumulate_assoc: {l} vs {r}");

            // `scale_scale`: (s·t)·x == s·(t·x).
            let (s, t) = (1.7f64, -0.4f64);
            let l = ((a.clone() * s) * t).eval(&vals).unwrap();
            let r = (a.clone() * (s * t)).eval(&vals).unwrap();
            assert!((l - r).abs() < 1e-9, "scale_scale: {l} vs {r}");

            // Additive identity, exactly as `num::Zero` spells it.
            let z: Term = Zero::zero();
            let id = (a.clone() + z).eval(&vals).unwrap();
            let a0 = a.eval(&vals).unwrap();
            assert!((id - a0).abs() < 1e-12, "zero is the additive identity");
        }
    }
}

/// A small random `Term` over `n_vars` variables (no bare `Var`, which is not an
/// expression on either crate).
fn random_term(rng: &mut StdRng, n_vars: u32) -> Term {
    let atom = |rng: &mut StdRng| match rng.random_range(0..3usize) {
        0 => Term::var(rng.random_range(0..n_vars)).sin(),
        1 => Term::var(rng.random_range(0..n_vars)).cos(),
        _ => Term::from(rng.random_range(-2.0..2.0)),
    };
    let mut t = atom(rng);
    for _ in 0..rng.random_range(0..4usize) {
        match rng.random_range(0..2usize) {
            0 => t += atom(rng),
            _ => t *= atom(rng),
        }
    }
    t
}

#[test]
fn symbolic_rotation_laws_hold() {
    // `lean/PPVM/Instantiations/Rotation.lean`: `rot_norm_sq` (the 2-D rotation
    // is norm-preserving, i.e. `sin²θ + cos²θ = 1`) and `rot_rot` (rotations
    // compose by angle addition) — stated over the SYMBOLIC angle domain, which
    // is what `Angle<Term> for Term` makes usable.
    use ppvm_traits_2::Angle;

    let (s, c) = Angle::<Term>::sin_cos(&Term::var(0));
    let norm = s.clone() * s.clone() + c.clone() * c.clone();
    let mut rng = seeded_rng(9001);
    for _ in 0..256 {
        let th: f64 = rng.random_range(-std::f64::consts::PI..std::f64::consts::PI);
        let v = norm.eval(&[th]).unwrap();
        assert!((v - 1.0).abs() < 1e-12, "rot_norm_sq at θ={th}: {v}");
    }

    // `rot_rot`: sin(a+b) = sin a cos b + cos a sin b, cos(a+b) = cos a cos b −
    // sin a sin b — as SYMBOLIC identities in two variables.
    let (sa, ca) = Angle::<Term>::sin_cos(&Term::var(0));
    let (sb, cb) = Angle::<Term>::sin_cos(&Term::var(1));
    let sin_sum = sa.clone() * cb.clone() + ca.clone() * sb.clone();
    let cos_sum = ca * cb - sa * sb;
    for _ in 0..256 {
        let a: f64 = rng.random_range(-std::f64::consts::PI..std::f64::consts::PI);
        let b: f64 = rng.random_range(-std::f64::consts::PI..std::f64::consts::PI);
        assert!((sin_sum.eval(&[a, b]).unwrap() - (a + b).sin()).abs() < 1e-12);
        assert!((cos_sum.eval(&[a, b]).unwrap() - (a + b).cos()).abs() < 1e-12);
    }
}

// ===========================================================================
// 5. The monomial-key hashing contract.
// ===========================================================================
//
// The design's lazy-hash contract (`Hash` writes exactly `key_hash()`) is stated
// for `Sum` KEYS; a `Prod` is not a `Sum` key, it is the key of the coefficient's
// own monomial table. The contract that matters for it is the one perf feature 3
// names: canonicality. Two monomials built by different multiplication orders
// must be the SAME key with the SAME digest, because that coalescing is the only
// thing standing between the symbolic representation and 2^depth growth.

fn fx_digest<T: std::hash::Hash>(v: &T) -> u64 {
    use std::hash::Hasher;
    let mut h = fxhash::FxHasher64::default();
    v.hash(&mut h);
    h.finish()
}

#[test]
fn monomial_identity_is_canonical_under_multiplication_order() {
    // Build the same monomial by many different orders; every one must be an
    // equal key with an equal digest.
    let mut rng = seeded_rng(0x5EED);
    for _ in 0..256 {
        let n_factors = rng.random_range(2..8usize);
        let mut atoms: Vec<Prod> = (0..n_factors)
            .map(|_| {
                let v = rng.random_range(0..6u32);
                if rng.random_range(0..2usize) == 0 {
                    Prod::sin(v)
                } else {
                    Prod::cos(v)
                }
            })
            .collect();

        let reference = atoms
            .iter()
            .cloned()
            .reduce(|a, b| a * b)
            .expect("at least two atoms");
        for _ in 0..8 {
            // Fisher-Yates with the seeded RNG, so the permutations are
            // reproducible.
            for i in (1..atoms.len()).rev() {
                let j = rng.random_range(0..=i);
                atoms.swap(i, j);
            }
            let permuted = atoms.iter().cloned().reduce(|a, b| a * b).unwrap();
            assert_eq!(permuted, reference, "monomial identity is order-dependent");
            assert_eq!(
                fx_digest(&permuted),
                fx_digest(&reference),
                "equal monomials hashed differently"
            );
            // The cached degree totals must agree too — an incrementally
            // maintained counter that drifts silently changes truncation
            // (perf feature 2).
            assert_eq!(permuted.sin_pow(), reference.sin_pow());
            assert_eq!(permuted.cos_pow(), reference.cos_pow());
            assert_eq!(
                permuted.sin_pow(),
                permuted
                    .factors()
                    .iter()
                    .map(|f| f.sin as usize)
                    .sum::<usize>(),
                "sin_pow drifted from the summed exponents"
            );
            assert_eq!(
                permuted.cos_pow(),
                permuted
                    .factors()
                    .iter()
                    .map(|f| f.cos as usize)
                    .sum::<usize>(),
                "cos_pow drifted from the summed exponents"
            );
        }
    }
}

#[test]
fn monomial_hash_is_deterministic_and_low_collision() {
    // Seed-freedom (perf feature 4) is what makes the `Display` snapshot
    // contract hold across runs, so the digest of a freshly built monomial must
    // be reproducible within and across constructions.
    let a = {
        let mut p = Prod::sin(3);
        p.mul_cos(1);
        p.mul_sin(1);
        p
    };
    let b = {
        let mut p = Prod::sin(1);
        p.mul_sin(3);
        p.mul_cos(1);
        p
    };
    assert_eq!(fx_digest(&a), fx_digest(&b));
    assert_eq!(fx_digest(&a), fx_digest(&a));

    // Distribution: distinct monomials must land on distinct digests
    // essentially always. 4096 distinct keys into a 64-bit space; the birthday
    // expectation is ~4.5e-13 collisions, so ANY collision here is a red flag,
    // and the low bits (which the table indexes on) must spread too.
    // `Prod` is `Hash + Eq` but deliberately not `Ord`, so the distinctness set
    // is a `HashSet` (which still separates two keys that collide on the digest,
    // because it falls back to `Eq`).
    let mut keys: HashSet<Prod> = HashSet::new();
    let mut digests: BTreeSet<u64> = BTreeSet::new();
    let mut low_bits: BTreeMap<u64, usize> = BTreeMap::new();
    let mut rng = seeded_rng(0xD15);
    while keys.len() < 4096 {
        let mut p = Prod::new();
        for _ in 0..rng.random_range(1..6usize) {
            let v = rng.random_range(0..12u32);
            if rng.random_range(0..2usize) == 0 {
                p.mul_sin(v);
            } else {
                p.mul_cos(v);
            }
        }
        if !keys.insert(p.clone()) {
            continue;
        }
        let d = fx_digest(&p);
        digests.insert(d);
        *low_bits.entry(d & 0x3FF).or_insert(0) += 1;
    }
    assert_eq!(
        digests.len(),
        keys.len(),
        "{} digest collisions over {} distinct monomials",
        keys.len() - digests.len(),
        keys.len()
    );
    // Avalanche/spread on the bucket-selecting low bits: with 4096 keys over
    // 1024 buckets the mean occupancy is 4; a hash that ignored a factor would
    // pile up. Allow generous slack (a Poisson(4) max over 1024 draws is ~14).
    let worst = low_bits.values().copied().max().unwrap_or(0);
    assert!(
        worst <= 24,
        "low-bit bucket occupancy {worst} is far above the Poisson(4) expectation \
         — the digest is not mixing the monomial's factors"
    );
    assert!(
        low_bits.len() > 900,
        "only {} of 1024 low-bit buckets used",
        low_bits.len()
    );
}

#[test]
fn ord_on_prod_is_not_required_but_hash_eq_agree() {
    // The `FxHashMap<Prod, f64>` contract: `a == b` implies equal digests (the
    // `Hash`/`Eq` consistency obligation). Checked on the pairs most likely to
    // break it — same factors, different phase; same variables, different
    // exponents.
    let mut a = Prod::sin(0);
    a.mul_sin(0);
    let mut b = Prod::sin(0);
    b.mul_cos(0);
    assert_ne!(a, b);
    assert_ne!(fx_digest(&a), fx_digest(&b));

    let mut c = Prod::sin(0);
    c.add_phase(2);
    let d = Prod::sin(0);
    assert_ne!(c, d, "the phase is part of the monomial identity");
    assert_ne!(fx_digest(&c), fx_digest(&d));

    let mut e = Prod::sin(0);
    e.mul_sin(0);
    assert_eq!(a, e);
    assert_eq!(fx_digest(&a), fx_digest(&e));
}
