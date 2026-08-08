// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`Term`] as a `ppvm-traits-2` coefficient ring.
//!
//! Ported from `ppvm-sym/src/coeff.rs`, re-cut along the `-2` trait split
//! (Design: `traits-2-configuration-and-hashing.md` §"Coefficient, angle, and
//! truncation"). Old bundled four responsibilities into one
//! `ppvm_traits::Coefficient`; here they land in four places:
//!
//! | old `Coefficient` method | `-2` home |
//! |---|---|
//! | `sin_cos` | [`Angle<Term> for Term`] — the *angle* domain |
//! | `cutoff` | [`Coefficient::magnitude`] + the `Policy` threshold |
//! | `half` | [`Halvable`](ppvm_traits_2::Halvable) — **not implemented**, see below |
//! | `mul_sign` | [`Coefficient::mul_sign`], unchanged |
//!
//! and old's separate `ComplexCoefficient::mul_phase` becomes
//! [`ImaginaryUnit`] (+ [`Conjugate`] for the sesquilinear pairing).

use ppvm_traits_2::{Angle, Coefficient, Conjugate, ImaginaryUnit};

use crate::term::{Inner, Prod, Sum, Term};

/// Old's `Coefficient for Term`, minus the two methods the `-2` split moved out.
///
/// Design: §"Coefficient, angle, and truncation". The module laws this value
/// domain feeds are machine-checked in `lean/PPVM/Algebra/GradedMap.lean`
/// (`accumulate_comm`, `accumulate_assoc`, `scale_scale`).
impl Coefficient for Term {
    /// A **scalar** sign flip on the stored `f64` coefficients — never a ring
    /// multiply by a `Term` denoting `±1`.
    ///
    /// This is old's body verbatim (`clone` then `*= sign as f64`), and the
    /// scalar form is load-bearing: it is called once per branch term per
    /// rotation gate (`sin.mul_sign(ε)` in `ppvm-pauli-sum-2::rotation`), where a
    /// generic `Term × Term` would allocate and walk a map for a sign flip
    /// (integration baseline, perf feature 9).
    ///
    /// Like old it is **partial**: a bare [`Term::var`] receiver panics, because
    /// the scaling routes through `MulAssign<f64>` (`oldSuspectedBugs` #7).
    #[inline]
    fn mul_sign(&self, sign: i8) -> Self {
        let mut ret = self.clone();
        ret *= sign as f64;
        ret
    }

    #[inline]
    fn mul_sign_assign(&mut self, sign: i8) {
        *self *= sign as f64;
    }

    #[inline]
    fn add_assign_ref(&mut self, rhs: &Self) {
        self.add_ref(rhs);
    }

    /// `|c|` for a constant, `f64::INFINITY` for every symbolic form.
    ///
    /// # Behaviour parity, and the law it breaks
    ///
    /// This reproduces old's `Coefficient::cutoff` exactly: old returned
    /// `c.abs() < threshold` **only** for `Inner::Const` and `false` otherwise,
    /// so `CoefficientThreshold`'s `retain(|_, v| !v.cutoff(t))` never dropped a
    /// symbolic coefficient however tiny its monomial coefficients were. The `-2`
    /// keep-rule is `magnitude() >= threshold`, so returning `+∞` for the
    /// symbolic forms reproduces "never dropped" and `|c|` reproduces the
    /// constant case — sum-level coefficient-threshold truncation stays **inert**
    /// on symbolic coefficients (behavioural contract 3).
    ///
    /// It **violates** the documented `magnitude` law `N(x) == 0 iff x == 0`
    /// (an empty symbolic `Sum` denotes `0` yet reports `+∞`), and it is not
    /// multiplicative. That is not a slip, it is a forced choice, and the
    /// alternative is worse:
    ///
    /// * No absolute value `N: R[sᵢ, cᵢ] → ℝ` exists at all — the natural `ℓ¹`
    ///   coefficient norm is only *sub*-multiplicative (`(1+x)(1−x) = 1−x²`:
    ///   `ℓ¹` gives `2·2 = 4` vs `2`), so the `N(xy) = N(x)N(y)` clause is
    ///   unsatisfiable on this ring.
    /// * Returning the `ℓ¹` norm instead would satisfy subadditivity — which
    ///   `lean/PPVM/Algebra/Truncation.lean` `l1_bound_needs_subadditive` shows is
    ///   the load-bearing clause — but would start **dropping terms old kept**,
    ///   a user-visible behaviour change, so the prime directive rules it out.
    ///
    /// The open adjudication is now settled in the Lean spec, and it settles it
    /// *against* the law, not against parity: `l1_bound_seminorm`
    /// (`lean/PPVM/Algebra/Truncation.lean`) shows the `ℓ¹` truncation bound
    /// survives weakening `AbsoluteValue` to a nonnegative, `0`-vanishing,
    /// subadditive, **sub**-multiplicative seminorm — so an `ℓ¹` `magnitude`
    /// *would* be fully guaranteed — while `l1_bound_seminorm_needs_zero` shows
    /// the one clause that must survive is `N 0 = 0`, which is precisely the one
    /// `+∞` breaks. Parity wins under the prime directive, so this impl keeps
    /// `+∞` and the law exemption is explicit: `CoefficientThreshold` is inert on
    /// symbolic coefficients (old's behaviour) and carries **no** `ℓ¹` error
    /// bound here. See the crate-level `# Deferrals`.
    fn magnitude(&self) -> f64 {
        if let Inner::Const(c) = self.inner {
            c.abs()
        } else {
            f64::INFINITY
        }
    }
}

/// A [`Term`] is its own rotation-angle domain: `θ ↦ (sin θ, cos θ)`, both
/// already in the coefficient domain.
///
/// This is old's `Coefficient::sin_cos` verbatim — two clones, and it **panics**
/// on a `One`/`Sum` angle, so only `Var`/`Const` angles are usable. It is what
/// makes `sum.rz(0, Term::var(0))` compile on a symbolic sum, i.e. what
/// `examples/symbolic.rs` does.
///
/// Design: §"Coefficient, angle, and truncation" — this is the "symbolic angle"
/// instantiation the split was introduced for.
///
/// # The rotation laws hold only *after* `eval`
///
/// The coefficient ring driven here is the **free** polynomial ring `ℝ[sᵢ, cᵢ]`:
/// `sin(x).square() + cos(x).square()` is a two-monomial [`Sum`], not `1`, and no
/// `Term`-level operation ever reduces it — machine-checked as
/// `pythagorean_ne_one` in `lean/PPVM/Instantiations/Symbolic.lean`. The 2-D
/// rotation guarantees of `lean/PPVM/Instantiations/Rotation.lean`
/// (`rot_norm_sq`, `rot_rot`) consume exactly that relation, so they transfer to
/// this domain **only under evaluation**, pointwise in `θ`: `evalHom_symRot` is
/// the commuting square `ev ∘ symRot = rot θ ∘ ev`, `symRot_norm_sq_after_eval`
/// the transferred norm preservation, and `symRot_norm_sq_ne_symbolically` the
/// witness that the unqualified claim is false. Citing `rot_norm_sq` for the
/// symbolic ring without the `eval` qualifier would be unsound.
impl Angle<Term> for Term {
    #[inline]
    fn sin_cos(&self) -> (Term, Term) {
        (self.clone().sin(), self.clone().cos())
    }
}

/// A **real** `f64` angle driving a symbolic sum.
///
/// # Behaviour parity
///
/// Old's rotations took `theta: impl Into<T::Coeff>` with `Coefficient:
/// From<f64>`, so `sum.rx(0, 0.1)` compiled on a `Term`-coefficient sum: the
/// `f64` widened to `Term::from(0.1)` and `sin_cos` then constant-folded it to
/// `(Const(sin 0.1), Const(cos 0.1))`. [`RotationOne`](ppvm_traits_2::RotationOne)
/// cannot take `impl Into<A>` (`A` is a free parameter, so it would be
/// uninferable at the call site), so that caller spelling is preserved here
/// instead — exactly as `ppvm-traits-2` does for `Angle<Complex<f64>> for f64`.
impl Angle<Term> for f64 {
    #[inline]
    fn sin_cos(&self) -> (Term, Term) {
        let (s, c) = f64::sin_cos(*self);
        (Term::from_f64(s), Term::from_f64(c))
    }
}

impl Term {
    /// Multiply by `i^phase` — old's `ComplexCoefficient::mul_phase`.
    ///
    /// Reproduced arm for arm from `ppvm-sym/src/coeff.rs:32-76`, with the
    /// `oldSuspectedBugs` #3 correction: old's `Sum` arm phased the constant part
    /// by feeding a phase-only monomial to `add_term`, which short-circuited on
    /// `pow() == 0` and folded the value into `c0`, **throwing the phase away**.
    /// [`Sum::add_term`] no longer takes that short-circuit for a phase-carrying
    /// monomial, so the constant summand is phased like every other.
    ///
    /// The oracle is `lean/PPVM/Instantiations/Symbolic.lean`, which models this
    /// method as `phaseFold k` — the key relabelling `(m, j) ↦ (m, j + k)` — on
    /// the ℤ/4-graded `PhasedSymRing`, and proves
    /// `phaseFold_eq_iSym_pow_mul`: the relabelling **equals** the ring product
    /// `iᵏ · x`, with read-out corollary `evalC_phaseFold`
    /// (`evalC θ (mul_phase k x) = iᵏ · evalC θ x`). `phaseFold_const` is the
    /// `(0, 0) ↦ (0, k)` arm this `Sum` branch relies on, and
    /// `phaseFold_drop_const_ne` proves old's "leave the constant on key `(0, 0)`"
    /// is a *different function*, not another representation of the same one — so
    /// the divergence is forced. (Additivity of the twisted product,
    /// `lean/PPVM/Algebra/Twisted.lean` `twistedConv_add_left` /
    /// `twistedConv_add_right` and `iPow_add`, is the Pauli-key analogue and
    /// supports but does not by itself state the symbolic fold.)
    #[inline]
    pub fn mul_phase(&self, phase: u8) -> Self {
        match self.inner {
            Inner::Sum(ref s) => {
                let mut ret = Sum::new();
                if let Some(maps) = &s.maps {
                    for (p, c) in &maps.terms {
                        let mut new_p = p.clone();
                        new_p.add_phase(phase);
                        ret.add_term(new_p, *c, self.max_sin, self.min_eps);
                    }
                }
                let mut c0 = Prod::new();
                c0.add_phase(phase);
                ret.add_term(c0, s.c0, self.max_sin, self.min_eps);
                Term {
                    inner: Inner::Sum(ret),
                    max_sin: self.max_sin,
                    min_eps: self.min_eps,
                }
            }
            Inner::Const(f) => {
                let mut ret = Prod::new();
                ret.add_phase(phase);
                Term {
                    inner: Inner::One(ret, f),
                    max_sin: self.max_sin,
                    min_eps: self.min_eps,
                }
            }
            Inner::One(ref p, c) => {
                let mut ret = p.clone();
                ret.add_phase(phase);
                Term {
                    inner: Inner::One(ret, c),
                    max_sin: self.max_sin,
                    min_eps: self.min_eps,
                }
            }
            Inner::Var(_) => {
                panic!(
                    "variable is not used in sin/cos expressions, bare variable is not allowed in expression"
                );
            }
        }
    }
}

/// The `i` capability L4 needs — old's `ComplexCoefficient for Term`, re-cut.
///
/// # Law caveat (representational vs denotational equality)
///
/// [`ImaginaryUnit`] asks for `i·i == −one()`. Here `i·i` is
/// `One(Prod{phase: 2}, 1.0)` while `−one()` is `Const(−1.0)`: the two *denote*
/// the same value (checkable with [`Term::eval_complex`]) but are **not**
/// `PartialEq`, because `Term`'s equality is representational — old's
/// `#[derive(PartialEq)]` already made `Const(1.0) != One(Prod::new(), 1.0) !=
/// Sum{c0: 1.0, terms: {}}` (behavioural contract 5), and the engine's exact-map
/// equality depends on that. Normalizing representations to make the law hold
/// syntactically would change `PartialEq` results, so it is *not* done; the law
/// is discharged literally by the exact ring
/// [`GaussianInt`](crate::GaussianInt) instead, which is the L4 witness Phase 5
/// asks for.
///
/// That argument is no longer prose. Because `Prod`'s `phase` byte is part of
/// its `Hash`/`Eq`, the ring this type implements is the **ℤ/4-graded**
/// `PhasedSymRing = AddMonoidAlgebra ℝ (Mono × ZMod 4)`
/// (`lean/PPVM/Instantiations/Symbolic.lean`), and there:
///
/// * `iSym_sq_ne_neg_one` — `i·i` (key `phase 2`) and `−one()` (`−1` on key
///   `phase 0`) are **different elements**, exactly as `PartialEq` reports;
/// * `evalC_iSym_sq_eq_neg_one` — yet they have the same complex value, so the
///   exemption is sound rather than a papered-over bug;
/// * `evalC_not_injective` — the two coexist because [`Term::eval_complex`] is a
///   surjective-but-**not-injective** algebra hom (`evalC_mul` /`evalC_add`), so
///   representational equality on this ring is strictly finer than denotational
///   equality in `ℂ`. `phaseTwo_cancel_ne_zero` spells out the operational
///   consequence: two summands that cancel in `ℂ` are distinct keys here, so
///   `min_eps` thresholds them independently and `i²·p` never cancels `−p`.
///
/// Design: §"The map is a graded algebra over `C[K]`" (`ImaginaryUnit`) and
/// §"The symbolic coefficient ring is *free*"; Lean
/// `lean/PPVM/Pauli/Matrix.lean` `iU_sq`.
impl ImaginaryUnit for Term {
    fn imaginary_unit() -> Self {
        Term::from_f64(1.0).mul_phase(1)
    }

    /// The direct `i^1` fold, not the generic `self * imaginary_unit()`: it
    /// avoids building a throwaway `Term` and, on the `Sum` arm, avoids the
    /// whole-table rebuild the generic product would run.
    fn mul_i(&self) -> Self {
        self.mul_phase(1)
    }

    /// `iᵏ · self` folded **into the monomial's phase byte** — old's
    /// `ComplexCoefficient::mul_phase(k)` verbatim, for every `k` including
    /// `k = 0`.
    ///
    /// This override is load-bearing for behaviour parity, not a micro-
    /// optimization. The default body would take the `self.clone()` arm at
    /// `k = 0` and the `-self` arm at `k = 2`, whereas old *always* built a
    /// phase-carrying monomial: `Const(c).mul_phase(0)` is `One(i⁰, c)`, not
    /// `Const(c)`. `Term`'s `PartialEq` is representational and its `Display` is
    /// a snapshot contract (behavioural contracts 5 and 8), so the two spellings
    /// are user-visible: `sum += ("ZIII", Term::from(2.0)); sum *= Z_III`
    /// renders `"2.000 * "` on old and would render `"2"` under the default
    /// fold. Every engine-side phase fold ([`Phase::apply`](ppvm_traits_2::Phase::apply),
    /// hence `Sum::mul_word_assign` and `multiply_into`) routes through here, so
    /// the symbolic engine lands in old's representation.
    fn mul_i_pow(&self, k: u8) -> Self {
        self.mul_phase(k & 3)
    }
}

/// Complex conjugation on the symbolic ring: negate every monomial's phase
/// exponent (`i^k ↦ i^{-k}`); the `f64` coefficients are real and unchanged.
///
/// Law `conj(i) == −i` holds denotationally (`i^3 = −i`); as with
/// [`ImaginaryUnit`] above, `−i` written as `One(Prod{phase: 1}, −1.0)` is a
/// *different representation* of the same value than `One(Prod{phase: 3}, 1.0)`
/// (`conjSym_iSym_ne_neg_iSym`).
///
/// This impl has no old counterpart, so what makes phase negation the *right*
/// map rather than an arbitrary one is proved on the graded ring in
/// `lean/PPVM/Instantiations/Symbolic.lean`: `conjSym` is a ring involution
/// (`conjSym_conjSym`, and an `AlgHom` by construction) satisfying
/// `evalC ∘ conj == star ∘ evalC` (`evalC_conjSym`), whence `conj i == −i`
/// (`evalC_conjSym_iSym`). `lean/PPVM/Pauli/Matrix.lean` `star_iU` is the
/// corresponding fact about the 2×2 matrix `iU`; `evalC_conjSym` is the
/// coefficient-ring statement this impl actually owes.
impl Conjugate for Term {
    fn conj(&self) -> Self {
        match self.inner {
            Inner::Const(_) => self.clone(),
            Inner::One(ref p, c) => {
                let mut q = p.clone();
                q.phase = (4 - q.phase) % 4;
                Term {
                    inner: Inner::One(q, c),
                    max_sin: self.max_sin,
                    min_eps: self.min_eps,
                }
            }
            Inner::Sum(ref s) => {
                let mut ret = Sum::new();
                ret.c0 = s.c0;
                if let Some(maps) = &s.maps {
                    for (p, c) in &maps.terms {
                        let mut q = p.clone();
                        q.phase = (4 - q.phase) % 4;
                        ret.add_term(q, *c, self.max_sin, self.min_eps);
                    }
                }
                Term {
                    inner: Inner::Sum(ret),
                    max_sin: self.max_sin,
                    min_eps: self.min_eps,
                }
            }
            Inner::Var(_) => {
                panic!(
                    "variable is not used in sin/cos expressions, bare variable is not allowed in expression"
                );
            }
        }
    }
}

/// Folds from `Term::from_f64(0.0)`, so a read-out (`trace`/`overlap`) **resets
/// truncation to the defaults**: the accumulator is a fresh `Const(0.0)` with
/// `max_sin = usize::MAX` and `min_eps = f64::EPSILON`, regardless of what the
/// propagated coefficients carried (behavioural contract 10). The first `Sum`
/// summand's table is then adopted wholesale by the `Const + Sum` arm and becomes
/// the accumulator.
impl std::iter::Sum for Term {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut total = Term::from_f64(0.0);
        for t in iter {
            total += t;
        }
        total
    }
}

impl std::ops::Neg for Term {
    type Output = Term;

    #[inline]
    fn neg(self) -> Self::Output {
        let mut ret = self;
        ret *= -1.0;
        ret
    }
}

impl From<f32> for Term {
    #[inline]
    fn from(value: f32) -> Self {
        Term::from_f64(value as f64)
    }
}

impl From<f64> for Term {
    #[inline]
    fn from(value: f64) -> Self {
        Term::from_f64(value)
    }
}

impl From<i32> for Term {
    #[inline]
    fn from(value: i32) -> Self {
        Term::from_f64(value as f64)
    }
}

impl From<i64> for Term {
    #[inline]
    fn from(value: i64) -> Self {
        Term::from_f64(value as f64)
    }
}

/// `zero()` is `Const(0.0)`; `is_zero()` is `true` **only** for a constant below
/// `min_eps` — every non-`Const` form reports `false`, including an empty `Sum`
/// (which denotes `0`) and a `Sum` all of whose coefficients are `0`
/// (behavioural contract 4).
impl num::Zero for Term {
    fn zero() -> Self {
        Term::from_f64(0.0)
    }

    fn is_zero(&self) -> bool {
        if let Inner::Const(c) = self.inner {
            c.abs() < self.min_eps
        } else {
            false
        }
    }
}

/// `one()` is `Const(1.0)` — the mirror of [`num::Zero::zero`]'s `Const(0.0)`.
/// Old had no `num::One for Term`; it is required by
/// [`ImaginaryUnit`](ppvm_traits_2::ImaginaryUnit) and by `ppvm-pauli-sum-2`'s
/// noise kernels, and `Const(1.0)` is the only choice consistent with
/// `From<f64>`.
impl num::One for Term {
    fn one() -> Self {
        Term::from_f64(1.0)
    }
}

impl Term {
    /// Convenience for the `-2` phase surface: multiply by
    /// [`Phase`](ppvm_traits_2::Phase).
    ///
    /// Exactly `phase.apply(self)`: [`Phase::apply`](ppvm_traits_2::Phase::apply)
    /// delegates to [`ImaginaryUnit::mul_i_pow`](ppvm_traits_2::ImaginaryUnit::mul_i_pow),
    /// which this ring overrides to keep the `iᵏ` *in the monomial* for every
    /// `k` — old's `mul_phase` representation, including the `k = 0` promotion
    /// of a `Const` to a `One`.
    pub fn mul_pauli_phase(&self, phase: ppvm_traits_2::Phase) -> Self {
        self.mul_phase(phase.exponent())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::{Complex, One, Zero};

    #[test]
    fn magnitude_is_inert_on_symbolic_forms() {
        // Behavioural contract 3: a symbolic coefficient is never dropped by a
        // coefficient threshold, however tiny its monomial coefficients.
        let tiny = Term::var(0).sin() * 1e-30;
        assert_eq!(tiny.magnitude(), f64::INFINITY);
        assert!(tiny.magnitude() >= 1e-6);
        // …but a literal constant is.
        assert_eq!(Term::from_f64(1e-30).magnitude(), 1e-30);
        assert!(Term::from_f64(1e-30).magnitude() < 1e-6);
    }

    #[test]
    fn is_zero_only_for_constants() {
        // Behavioural contract 4.
        assert!(Term::from_f64(0.0).is_zero());
        let mut t = Term::var(0).sin();
        t += Term::var(0).sin() * -1.0;
        assert!(!t.is_zero());
        assert!(t != Term::from_f64(0.0));
    }

    #[test]
    fn partial_eq_is_representational() {
        // Behavioural contract 5: all three denote 1, none are equal.
        let a = Term::from_f64(1.0);
        let b = Term {
            inner: Inner::One(Prod::new(), 1.0),
            max_sin: usize::MAX,
            min_eps: f64::EPSILON,
        };
        let mut s = Sum::new();
        s.c0 = 1.0;
        let c = Term {
            inner: Inner::Sum(s),
            max_sin: usize::MAX,
            min_eps: f64::EPSILON,
        };
        assert!(a != b);
        assert!(b != c);
        assert!(a != c);
    }

    #[test]
    fn phase_is_part_of_monomial_identity() {
        // Behavioural contract 5: `sin(x)` and `i·sin(x)` never coalesce.
        let p = Prod::sin(0);
        let mut q = Prod::sin(0);
        q.add_phase(1);
        assert!(p != q);
    }

    #[test]
    fn conversions_and_neg() {
        assert_eq!(Term::from(2i32), Term::from(2.0f64));
        assert_eq!(Term::from(2i64), Term::from(2.0f64));
        assert_eq!(Term::from(2.0f32), Term::from(2.0f64));
        assert_eq!(-Term::from(2.0), Term::from(-2.0));
    }

    #[test]
    fn i_squared_denotes_minus_one() {
        // The `ImaginaryUnit` law, checked denotationally (see the impl's law
        // caveat): `i·i` and `−one()` are equal as values, not as
        // representations.
        let i = Term::imaginary_unit();
        let sq = i.clone() * i;
        let minus_one: Term = -Term::one();
        assert_eq!(
            sq.eval_complex(&[]).unwrap(),
            minus_one.eval_complex(&[]).unwrap()
        );
        assert_eq!(sq.eval_complex(&[]).unwrap(), Complex::new(-1.0, 0.0));
    }

    #[test]
    fn conj_of_i_denotes_minus_i() {
        // Lean `star_iU`.
        let i = Term::imaginary_unit();
        assert_eq!(i.conj().eval_complex(&[]).unwrap(), Complex::new(0.0, -1.0));
    }

    #[test]
    fn mul_phase_keeps_the_constant_part_phased() {
        // Divergence from old (`oldSuspectedBugs` #3): old dropped the phase of
        // the constant summand.
        let t = Term::from_f64(2.0) + Term::var(0).sin();
        let phased = t.mul_phase(1);
        let v = phased.eval_complex(&[0.5]).unwrap();
        assert!((v.re).abs() < 1e-15, "real part should vanish: {v}");
        assert!((v.im - (2.0 + 0.5f64.sin())).abs() < 1e-12, "{v}");
    }

    #[test]
    fn sin_cos_is_the_angle_domain() {
        let (s, c) = Angle::<Term>::sin_cos(&Term::var(0));
        assert_eq!(s, Term::var(0).sin());
        assert_eq!(c, Term::var(0).cos());
        let (s, c) = Angle::<Term>::sin_cos(&0.3f64);
        assert_eq!(s, Term::from_f64(0.3f64.sin()));
        assert_eq!(c, Term::from_f64(0.3f64.cos()));
    }

    #[test]
    #[should_panic(expected = "only variable or constant can be input of sin")]
    fn sin_of_a_compound_panics() {
        let _ = Term::var(0).sin().sin();
    }

    #[test]
    #[should_panic(expected = "only variable or constant can be input of cos")]
    fn cos_of_a_compound_panics() {
        let _ = Term::var(0).cos().cos();
    }

    #[test]
    #[should_panic(expected = "bare variable is not allowed")]
    fn arithmetic_on_a_bare_variable_panics() {
        let _ = Term::var(0) * Term::from_f64(2.0);
    }

    #[test]
    #[should_panic(expected = "bare variable is not allowed")]
    fn mul_sign_on_a_bare_variable_panics() {
        // `oldSuspectedBugs` #7: `mul_sign` routes through `*= sign as f64`.
        let _ = Term::var(0).mul_sign(1);
    }

    #[test]
    #[should_panic(expected = "bare variable is not allowed")]
    fn mul_phase_on_a_bare_variable_panics() {
        let _ = Term::var(0).mul_phase(1);
    }

    #[test]
    #[should_panic(expected = "only variable or constant can be input of sin")]
    fn sin_cos_on_a_compound_angle_panics() {
        let t = Term::var(0).sin();
        let _ = Angle::<Term>::sin_cos(&t);
    }
}
