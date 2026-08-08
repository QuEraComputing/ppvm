// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Multiplication: the ring product on [`Prod`], [`Sum`] and [`Term`].
//!
//! Ported from `ppvm-sym/src/mul.rs`, arm for arm. Every **operand-adoption**
//! (move-not-clone) arm old had is kept — `Const × Sum` re-tags `self` and scales
//! in place, `One × Sum` *moves* the right-hand table and multiplies it by the
//! single monomial, and the `Sum × Sum` arm skips the two cross loops when the
//! respective `c₀` is below `min_eps` (integration baseline, perf feature 6).
//! The operators are **by value** for exactly that reason: a uniform
//! `&self × &rhs → new` signature would turn an `O(1)` adoption into an
//! `O(#monomials)` deep copy on every gate.

use crate::term::{Factor, Inner, Prod, Sum, Term, merge_factors};

// `pauli_sum *= Term` / `pauli_sum * Term` — old's
// `impl_op_mul_assign_coefficient!(Term)` (`ppvm-sym/src/mul.rs:13`), which is
// the *only* reason old `ppvm-sym` depended on `ppvm-pauli-sum` at all.
//
// The impl cannot be written here by hand: `Sum<S, P>` is foreign to this crate
// and its type parameters precede the local `Term`, so the orphan rule rejects
// it. Instantiating the engine's exported macro places the impl in
// `ppvm-pauli-sum-2`'s expansion, which is legal — the same escape hatch old
// used. Semantically it is `Sum::scale`, i.e. old's `*v *= value.clone()` per
// coefficient, so the capability was never lost; only this spelling was.
ppvm_pauli_sum_2::impl_scalar_mul!(Term);

impl Prod {
    /// Multiply this product by an additional `sin(x_u)`.
    ///
    /// The cached `sin_pow` total is updated incrementally here (perf feature 2)
    /// so [`Sum::add_term`]'s `max_sin` test stays `O(1)`.
    pub fn mul_sin(&mut self, u: u32) {
        self.merge_factor(Factor {
            var: u,
            sin: 1,
            cos: 0,
        });
        self.sin_pow += 1;
        self.debug_check();
    }

    /// Multiply this product by an additional `cos(x_u)`.
    pub fn mul_cos(&mut self, u: u32) {
        self.merge_factor(Factor {
            var: u,
            sin: 0,
            cos: 1,
        });
        self.cos_pow += 1;
        self.debug_check();
    }
}

impl Sum {
    /// Multiply this sum in place by `coeff · p`, respecting the same
    /// `max` / `min_eps` truncation bounds used elsewhere.
    ///
    /// # The two whole-sum-to-zero shortcuts
    ///
    /// If `p.sin_pow() > max` or `|coeff| < min_eps`, *every* product monomial
    /// would be dropped by [`Sum::add_term`] (each has
    /// `sin_pow ≥ p.sin_pow()`), so the whole table is cleared in one `clear()`
    /// instead of walking it (integration baseline, perf feature 8).
    ///
    /// **The two arms are not the same kind of statement.** The *degree* arm is
    /// exact, machine-checked as `mulMono_clear_sound` / `mulMono_retain_clear`
    /// in `lean/PPVM/Instantiations/Symbolic.lean`: the sine degree is additive
    /// (`sinDeg_add`), so an over-degree multiplier puts every product monomial
    /// in the truncation ideal and the truncated product really is the zero map.
    /// The `min_eps` arm is **not** exact — it discards monomials the
    /// per-monomial rule in [`Sum::add_term`] would keep, because a large stored
    /// coefficient can rescue a small multiplier
    /// (`epsClear_ne_retain_pointwise`, same file: one entry of magnitude `1e6`
    /// times `c = 1e-13` under `min_eps = 1e-12`). It is kept anyway because it
    /// is an `ℓ¹`-controlled *over*-truncation, not an unsound one: the mass it
    /// discards is exactly `|c|·ℓ¹` of the table (`epsClear_l1_eq`), hence
    /// strictly under `min_eps·ℓ¹` (`epsClear_l1_lt`), and the read-out error is
    /// bounded by the same quantity (`epsClear_error_lt`, via
    /// `PPVM.Truncation.l1_bound`). Do **not** "simplify" this branch into a
    /// per-monomial loop: that changes results, and the degree-arm theorems do
    /// not cover it.
    ///
    /// Old leaves
    /// an **empty `Sum`**, not a `Const(0.0)` — so the result is not
    /// `is_zero()`, is `!=` to `Term::from(0.0)`, and prints as `[]`
    /// (behavioural contract 9). That is reproduced exactly.
    ///
    /// # The rebuild buffer
    ///
    /// Old did `let mut old_terms = std::mem::take(&mut self.terms)`, leaving a
    /// zero-capacity map to re-grow from scratch on every symbolic multiply.
    /// Here the persistent `aux` buffer is swapped in instead, so both
    /// allocations survive the rebuild (integration baseline, perf feature 5 —
    /// the crate's aux-double-buffer gap).
    #[inline]
    pub fn mul_term(&mut self, p: Prod, coeff: f64, max: usize, min_eps: f64) {
        if p.sin_pow() > max || coeff.abs() < min_eps {
            if let Some(maps) = &mut self.maps {
                maps.terms.clear();
            }
            self.c0 = 0.0;
            return;
        }

        // A phase-carrying monomial is not a scalar even when `pow() == 0`
        // (`iPow_add`, `lean/PPVM/Algebra/Twisted.lean`); see `Sum::add_term`.
        if p.pow() == 0 && p.phase() == 0 {
            *self *= coeff;
            return;
        }

        // Ping-pong: `aux` is empty on entry, so after the swap `terms` holds
        // `aux`'s retained capacity and `aux` holds the table to rebuild from.
        let mut maps = self.maps.take().unwrap_or_default();
        debug_assert!(maps.aux.is_empty());
        std::mem::swap(&mut maps.terms, &mut maps.aux);
        let mut old_terms = std::mem::take(&mut maps.aux);
        maps.terms.reserve(old_terms.len());
        self.maps = Some(maps);

        let c0 = self.c0;
        self.add_term(p.clone(), c0 * coeff, max, min_eps);
        for (k, v) in old_terms.drain() {
            self.add_term(k * p.clone(), v * coeff, max, min_eps);
        }
        self.c0 = 0.0;
        // Hand the (now empty, still-allocated) buffer back for the next multiply.
        self.maps.as_mut().unwrap().aux = old_terms;
    }
}

impl std::ops::Mul<Prod> for Prod {
    type Output = Prod;

    #[inline]
    fn mul(self, rhs: Prod) -> Self::Output {
        let mut new = self;
        new *= rhs;
        new
    }
}

/// Monomial multiplication: exponents add per variable and the **phase
/// composes**.
///
/// # Divergence from old (`oldSuspectedBugs` #2)
///
/// Old (`ppvm-sym/src/mul.rs:63-73`) merged the sin/cos maps and the degree
/// totals but never combined `phase`: `self.phase` was left untouched and
/// `rhs.phase` discarded, so `(i·P)·(i·Q)` yielded `i` rather than `i² = −1` —
/// and, because `phase` is part of `Prod`'s `Hash`/`Eq`, mis-keyed the product
/// monomial. The `ℤ/4` phase group is homomorphic under monomial multiplication:
/// `i^a · i^b = i^{a+b}` (`lean/PPVM/Algebra/Twisted.lean` `iPow_add`, with
/// `tmul_assoc`/`gtmul_assoc` requiring it for associativity, and
/// `lean/PPVM/Pauli/Phase.lean` `phaseExp_eq_ref` for the packed encoding).
/// `ppvm-traits-2::Phase::compose` is the correct model and is what this
/// reproduces.
impl std::ops::MulAssign<Prod> for Prod {
    #[inline]
    fn mul_assign(&mut self, rhs: Prod) {
        self.phase = (self.phase + rhs.phase) % 4;
        self.sin_pow += rhs.sin_pow;
        self.cos_pow += rhs.cos_pow;

        if rhs.factors.is_empty() {
            self.debug_check();
            return;
        }
        if self.factors.is_empty() {
            self.factors = rhs.factors;
            self.debug_check();
            return;
        }
        if rhs.factors.len() == 1 {
            // The overwhelmingly common case on the propagation hot path: the
            // right operand is a single `sin(x_u)` or `cos(x_u)` atom.
            self.merge_factor(rhs.factors[0]);
            self.debug_check();
            return;
        }
        self.factors = merge_factors(&self.factors, &rhs.factors);
        self.debug_check();
    }
}

impl std::ops::MulAssign<f64> for Sum {
    #[inline]
    fn mul_assign(&mut self, rhs: f64) {
        self.c0 *= rhs;
        if let Some(maps) = &mut self.maps {
            for v in maps.terms.values_mut() {
                *v *= rhs;
            }
        }
    }
}

impl std::ops::MulAssign<f64> for Term {
    #[inline]
    fn mul_assign(&mut self, rhs: f64) {
        match self.inner {
            Inner::Sum(ref mut s) => {
                *s *= rhs;
            }
            Inner::One(_, ref mut c) => {
                *c *= rhs;
            }
            Inner::Var(_) => {
                panic!(
                    "variable is not used in sin/cos expressions,\
                    bare variable is not allowed in expression"
                );
            }
            Inner::Const(ref mut c) => {
                *c *= rhs;
            }
        }
    }
}

impl std::ops::MulAssign<Term> for Term {
    #[inline]
    fn mul_assign(&mut self, rhs: Term) {
        match self.inner {
            Inner::Sum(ref mut s1) => match rhs.inner {
                Inner::Sum(ref s2) => {
                    let mut new_sum = Sum::new();
                    new_sum.c0 = s1.c0 * s2.c0;
                    // Skip the `s1.terms × s2.c0` cross loop when `s2.c0` is
                    // below `min_eps` (perf feature 6).
                    if s2.c0.abs() > self.min_eps
                        && let Some(maps) = &s1.maps
                    {
                        for (p1, c1) in &maps.terms {
                            let p = p1.clone();
                            new_sum.add_term(p, c1 * s2.c0, self.max_sin, self.min_eps);
                        }
                    }

                    if s1.c0.abs() > self.min_eps
                        && let Some(maps) = &s2.maps
                    {
                        for (p2, c2) in &maps.terms {
                            let p = p2.clone();
                            new_sum.add_term(p, c2 * s1.c0, self.max_sin, self.min_eps);
                        }
                    }

                    if let (Some(maps1), Some(maps2)) = (&s1.maps, &s2.maps) {
                        for (p1, c1) in &maps1.terms {
                            for (p2, c2) in &maps2.terms {
                                new_sum.add_term(
                                    p1.clone() * p2.clone(),
                                    c1 * c2,
                                    self.max_sin,
                                    self.min_eps,
                                );
                            }
                        }
                    }
                    self.inner = Inner::Sum(new_sum);
                }
                Inner::Const(c) => {
                    *s1 *= c;
                }
                Inner::One(p2, c2) => {
                    s1.mul_term(p2, c2, self.max_sin, self.min_eps);
                }
                Inner::Var(_) => {
                    panic!(
                        "variable is not used in sin/cos expressions, bare variable is not allowed in expression"
                    );
                }
            },
            Inner::Const(c) => match rhs.inner {
                // Operand adoption: the large right-hand table is MOVED in.
                Inner::Sum(s) => {
                    self.inner = Inner::Sum(s);
                    *self *= c;
                }
                Inner::Const(c2) => {
                    self.inner = Inner::Const(c * c2);
                }
                Inner::One(p2, c2) => {
                    self.inner = Inner::One(p2, c * c2);
                }
                Inner::Var(_) => {
                    panic!(
                        "variable is not used in sin/cos expressions, bare variable is not allowed in expression"
                    );
                }
            },
            Inner::One(ref p, c) => match rhs.inner {
                // Operand adoption again: move `s`, then fold the single
                // monomial through it in place.
                Inner::Sum(s) => {
                    let mut new_sum = s;
                    new_sum.mul_term(p.clone(), c, self.max_sin, self.min_eps);
                    self.inner = Inner::Sum(new_sum);
                }
                Inner::Const(c2) => {
                    self.inner = Inner::One(p.clone(), c * c2);
                }
                Inner::One(p2, c2) => {
                    self.inner = Inner::One(p.clone() * p2, c * c2);
                }
                Inner::Var(_) => {
                    panic!(
                        "variable is not used in sin/cos expressions, bare variable is not allowed in expression"
                    );
                }
            },
            Inner::Var(_) => {
                panic!(
                    "variable is not used in sin/cos expressions, bare variable is not allowed in expression"
                );
            }
        }
    }
}

impl std::ops::Mul<f64> for Term {
    type Output = Term;

    #[inline]
    fn mul(self, rhs: f64) -> Self::Output {
        let mut ret = self;
        ret *= rhs;
        ret
    }
}

impl std::ops::Mul<Term> for f64 {
    type Output = Term;

    #[inline]
    fn mul(self, rhs: Term) -> Self::Output {
        let mut ret = rhs;
        ret *= self;
        ret
    }
}

impl std::ops::Mul<Term> for Term {
    type Output = Term;

    #[inline]
    fn mul(self, rhs: Term) -> Self::Output {
        let mut ret = self;
        ret *= rhs;
        ret
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn monomial_canonical_under_multiplication_order() {
        // Same monomial built by different multiplication orders must be the
        // SAME key (perf feature 3: canonicality is what stops the term count
        // exploding).
        let mut a = Prod::sin(3);
        a *= Prod::cos(1);
        a *= Prod::sin(1);

        let mut b = Prod::sin(1);
        b *= Prod::sin(3);
        b *= Prod::cos(1);

        assert_eq!(a, b);
        assert_eq!(hash_of(&a), hash_of(&b));
        assert_eq!(a.sin_pow(), 2);
        assert_eq!(a.cos_pow(), 1);
    }

    fn hash_of(p: &Prod) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut h = fxhash::FxHasher64::default();
        p.hash(&mut h);
        h.finish()
    }

    #[test]
    fn phase_composes_mod_four() {
        // Old dropped the phase here entirely (`oldSuspectedBugs` #2).
        let mut a = Prod::sin(0);
        a.add_phase(1);
        let mut b = Prod::cos(0);
        b.add_phase(3);
        let c = a * b;
        assert_eq!(c.phase(), 0);
    }

    #[test]
    fn mul_term_zero_shortcut_leaves_empty_sum() {
        // Behavioural contract 9: the shortcut leaves an EMPTY `Sum`.
        let mut s = Sum::new();
        s.add_term(Prod::cos(0), 1.0, usize::MAX, f64::EPSILON);
        let mut p = Prod::sin(0);
        p.mul_sin(1);
        s.mul_term(p, 1.0, 1, f64::EPSILON);
        assert!(s.is_empty());
        assert_eq!(s.c0(), 0.0);
    }

    #[test]
    fn mul_term_reuses_the_aux_allocation() {
        // Perf feature 5: after a rebuild both buffers are still allocated.
        let mut s = Sum::new();
        for i in 0..16u32 {
            s.add_term(Prod::cos(i), 1.0, usize::MAX, f64::EPSILON);
        }
        s.mul_term(Prod::sin(0), 1.0, usize::MAX, f64::EPSILON);
        let maps = s.maps.as_ref().expect("first monomial creates the map box");
        assert!(maps.aux.is_empty());
        assert!(maps.aux.capacity() > 0, "aux allocation was freed");
        assert_eq!(s.len(), 16);
    }

    #[test]
    fn one_times_one_stays_one() {
        // Perf feature 1: the single-monomial fast form never promotes.
        let t = Term::var(0).sin() * Term::var(1).cos();
        assert!(matches!(t.inner(), Inner::One(..)));
    }

    #[test]
    fn const_times_anything_stays_out_of_the_map() {
        let t = Term::from_f64(2.0) * Term::var(0).sin();
        assert!(matches!(t.inner(), Inner::One(_, c) if *c == 2.0));
    }
}
