// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Evaluation / substitution: bind every variable to a value and fold the
//! symbolic expression to a number.
//!
//! Ported from `ppvm-sym/src/eval.rs`. The **fold order is part of the
//! contract** and is reproduced exactly: within a monomial, the `sin` factors in
//! ascending variable order and then the `cos` factors in ascending variable
//! order, which fixes the floating-point association (behavioural contract 7).
//! The packed-vector layout preserves that order because the factors are stored
//! ascending by variable.

use anyhow::Result;
use num::Complex;

use crate::term::{Inner, Prod, Sum, Term};

/// How many variables [`AngleCache`] memoizes inline.
///
/// 32 covers every workload the integration baseline names (the deepest,
/// `sym.tfim.trotter` at 10 layers, uses 20 variables). Beyond it the cache
/// degrades gracefully to old's behaviour — recompute per access — rather than
/// allocating.
const INLINE_VARS: usize = 32;

/// A memo of `sin(vals[i])` / `cos(vals[i])`, shared across every monomial of
/// one [`Sum::eval`] call.
///
/// # Why
///
/// `sym.expectation.grid` is the VQE-shaped read-out: propagate once, then sweep
/// a 1000-point angle grid. Without this, [`Prod::eval`] recomputes
/// `vals[k].sin()` / `.cos()` *once per factor per monomial per grid point*, so
/// the transcendental dominates and the packed-monomial layout is invisible on
/// that workload. The trigonometric functions are pure, so memoizing them is
/// **bit-identical**: `res *= s.powi(e)` keeps exactly the same operands in
/// exactly the same fold order (behavioural contract 7).
///
/// One slot per variable, one validity bit, both trigonometric values filled on
/// a miss: a variable the sum never mentions is never evaluated, so the cost is
/// at most two transcendentals per *distinct variable of the whole sum*, against
/// old's two per factor per monomial. (The only shape that loses is a map-backed
/// sum holding a single monomial that uses each of its variables once — `2V` vs
/// `V`; a one-monomial coefficient is normally `Inner::One`, which goes through
/// the uncached [`Prod::eval`].) Three alternatives were measured on
/// `sym/expectation_eval` (new/old = 0.40 here) and are recorded so they are not
/// re-tried: independent `sin`/`cos` slots with two masks (+6.8% — twice the
/// working set to save a transcendental the sums ask for anyway); a 16-slot
/// window (+0.8%, and it stops covering a 10-layer Trotter circuit); and filling
/// the whole `vals[..=k]` prefix on a miss (−1.6% here, but it evaluates
/// variables the sum never mentions, so a term over a few *high* variable ids
/// would pay up to 32 useless transcendental pairs — rejected for the worst
/// case, not the average).
///
/// The bound check is deliberately kept *per access*, in traversal order, so the
/// `variable %{k} not found` error names exactly the variable old's `vals.get(k)`
/// would have tripped on first — the cache never touches a variable the fold
/// would not have reached.
///
/// No allocation and no interior mutability (perf feature 10): the cache is a
/// stack array owned by the one `eval` call that fills it.
struct AngleCache<'a> {
    vals: &'a [f64],
    /// Bit `i` set <=> `pairs[i]` is valid.
    valid: u32,
    pairs: [(f64, f64); INLINE_VARS],
}

impl<'a> AngleCache<'a> {
    #[inline]
    fn new(vals: &'a [f64]) -> Self {
        Self {
            vals,
            valid: 0,
            pairs: [(0.0, 0.0); INLINE_VARS],
        }
    }

    /// `(sin(vals[k]), cos(vals[k]))`, or old's error if `k` is out of range.
    #[inline]
    fn get(&mut self, k: u32) -> Result<(f64, f64)> {
        let i = k as usize;
        let v = *self
            .vals
            .get(i)
            .ok_or_else(|| anyhow::anyhow!("variable %{k} not found"))?;
        if i >= INLINE_VARS {
            return Ok((v.sin(), v.cos()));
        }
        if self.valid & (1 << i) == 0 {
            self.pairs[i] = (v.sin(), v.cos());
            self.valid |= 1 << i;
        }
        Ok(self.pairs[i])
    }
}

impl Prod {
    /// Evaluate the product at the variable assignment `vals`, where
    /// `vals[i]` is the value substituted for variable `i`.
    ///
    /// Returns `Ok(1.0)` for the empty product. An out-of-range variable is an
    /// `Err` carrying old's message `variable %{k} not found` — **not** a panic
    /// and **not** a silent zero.
    ///
    /// # Phase
    ///
    /// The `i^k` phase is deliberately **not** applied: the return type is `f64`,
    /// and old's `Prod::eval` ignored it too (`oldSuspectedBugs` #4). Keeping
    /// that here preserves every real-valued golden master. Use
    /// [`Prod::eval_complex`] to observe the phase.
    pub fn eval(&self, vals: &[f64]) -> Result<f64> {
        if self.pow() == 0 {
            return Ok(1.0);
        }

        let mut res = 1.0;
        if self.sin_pow() > 0 {
            for f in &self.factors {
                if f.sin == 0 {
                    continue;
                }
                let k = f.var;
                res *= vals
                    .get(k as usize)
                    .ok_or_else(|| anyhow::anyhow!("variable %{k} not found"))?
                    .sin()
                    .powi(f.sin as i32);
            }
        }

        if self.cos_pow() > 0 {
            for f in &self.factors {
                if f.cos == 0 {
                    continue;
                }
                let k = f.var;
                res *= vals
                    .get(k as usize)
                    .ok_or_else(|| anyhow::anyhow!("variable %{k} not found"))?
                    .cos()
                    .powi(f.cos as i32);
            }
        }
        Ok(res)
    }

    /// [`Prod::eval`] against a shared, already-warm [`AngleCache`].
    ///
    /// Bit-identical to `eval` by construction: same operands (the cache holds
    /// exactly `vals[k].sin()` / `vals[k].cos()`), same skip conditions, same
    /// left-to-right fold, same per-access bound check and therefore the same
    /// error. Only the transcendental is hoisted out of the monomial loop.
    /// `eval_cached_agrees_with_eval_bit_for_bit` pins that.
    #[inline]
    fn eval_cached(&self, cache: &mut AngleCache<'_>) -> Result<f64> {
        if self.pow() == 0 {
            return Ok(1.0);
        }

        let mut res = 1.0;
        if self.sin_pow() > 0 {
            for f in &self.factors {
                if f.sin == 0 {
                    continue;
                }
                res *= cache.get(f.var)?.0.powi(f.sin as i32);
            }
        }

        if self.cos_pow() > 0 {
            for f in &self.factors {
                if f.cos == 0 {
                    continue;
                }
                res *= cache.get(f.var)?.1.powi(f.cos as i32);
            }
        }
        Ok(res)
    }

    /// Evaluate the product **including** its `i^k` phase.
    ///
    /// # Divergence from old (`oldSuspectedBugs` #4)
    ///
    /// Old had no such method: its only `eval` returned `f64` and dropped the
    /// phase, which made `Term`'s whole `ComplexCoefficient` capability
    /// unobservable — an accumulated `i^k` evaluated as if `k == 0`. The phase is
    /// observable under the pairing (`lean/PPVM/Pauli/Matrix.lean` `star_iU`:
    /// `conj i = −i`) and carrying it is what makes the twisted product
    /// associative (`lean/PPVM/Algebra/Twisted.lean` `twistedConv_assoc`,
    /// `tmul_assoc`), so the phase-aware evaluation is added here rather than
    /// changing `eval`'s type.
    ///
    /// # Specification
    ///
    /// `lean/PPVM/Instantiations/Symbolic.lean` `evalC`: the `ℝ`-algebra
    /// homomorphism `PhasedSymRing → ℂ` sending `single (m, k) c` to
    /// `c · iᵏ · monoValue θ m`. `evalC_add` / `evalC_mul` are the hom laws —
    /// multiplicativity needs `Twisted.iPow_add` on the ℤ/4 grading and
    /// `monoValue_add` on the exponents, and it is what makes this function a
    /// valid oracle for `X·Y = iZ` over the symbolic ring. `evalC_not_injective`
    /// records the other half: distinct `Term`s can share a value, so
    /// `eval_complex` equality is strictly weaker than `PartialEq`.
    pub fn eval_complex(&self, vals: &[f64]) -> Result<Complex<f64>> {
        Ok(self.with_phase(self.eval(vals)?))
    }

    /// [`Prod::eval_complex`] against a shared [`AngleCache`].
    #[inline]
    fn eval_complex_cached(&self, cache: &mut AngleCache<'_>) -> Result<Complex<f64>> {
        Ok(self.with_phase(self.eval_cached(cache)?))
    }

    /// Multiply a real value by this monomial's `i^k`.
    #[inline]
    fn with_phase(&self, re: f64) -> Complex<f64> {
        match self.phase % 4 {
            0 => Complex::new(re, 0.0),
            1 => Complex::new(0.0, re),
            2 => Complex::new(-re, 0.0),
            _ => Complex::new(0.0, -re),
        }
    }
}

impl Sum {
    /// Evaluate the sum at `vals`, starting from `c₀` and folding the monomials
    /// in the backend's (deterministic, seed-free) iteration order.
    ///
    /// The per-variable `sin`/`cos` are computed **once per call** into an
    /// [`AngleCache`] rather than once per factor per monomial; see that type for
    /// why the result is bit-identical and the error shape unchanged.
    pub fn eval(&self, vals: &[f64]) -> Result<f64> {
        let mut cache = AngleCache::new(vals);
        let mut res = self.c0;
        for (p, c) in &self.terms {
            res += p.eval_cached(&mut cache)? * c;
        }
        Ok(res)
    }

    /// Phase-aware evaluation; see [`Prod::eval_complex`].
    pub fn eval_complex(&self, vals: &[f64]) -> Result<Complex<f64>> {
        let mut cache = AngleCache::new(vals);
        let mut res = Complex::new(self.c0, 0.0);
        for (p, c) in &self.terms {
            res += p.eval_complex_cached(&mut cache)? * *c;
        }
        Ok(res)
    }
}

impl Term {
    /// Evaluate this symbolic term at `vals`.
    ///
    /// A bare [`Term::var`] evaluates to `vals[u]` — one of the only two
    /// operations a bare variable supports (the other is `Display`).
    pub fn eval(&self, vals: &[f64]) -> Result<f64> {
        match self.inner {
            Inner::Const(c) => Ok(c),
            Inner::Var(u) => vals
                .get(u as usize)
                .copied()
                .ok_or_else(|| anyhow::anyhow!("variable %{u} not found")),
            Inner::One(ref p, c) => Ok(p.eval(vals)? * c),
            Inner::Sum(ref s) => s.eval(vals),
        }
    }

    /// Phase-aware evaluation; see [`Prod::eval_complex`].
    pub fn eval_complex(&self, vals: &[f64]) -> Result<Complex<f64>> {
        match self.inner {
            Inner::Const(c) => Ok(Complex::new(c, 0.0)),
            Inner::Var(u) => vals
                .get(u as usize)
                .map(|v| Complex::new(*v, 0.0))
                .ok_or_else(|| anyhow::anyhow!("variable %{u} not found")),
            Inner::One(ref p, c) => Ok(p.eval_complex(vals)? * c),
            Inner::Sum(ref s) => s.eval_complex(vals),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_product_is_one() {
        assert_eq!(Prod::new().eval(&[]).unwrap(), 1.0);
    }

    #[test]
    fn missing_variable_is_an_error_with_old_message() {
        let err = Term::var(5).sin().eval(&[0.1]).unwrap_err();
        assert!(
            err.to_string().contains("variable %5 not found"),
            "unexpected message: {err}"
        );
    }

    #[test]
    fn bare_variable_evaluates_to_its_binding() {
        assert_eq!(Term::var(3).eval(&[0.0, 1.0, 2.0, 9.0]).unwrap(), 9.0);
    }

    #[test]
    fn constant_folding_needs_no_variables() {
        assert_eq!(Term::from_f64(0.5).sin().eval(&[]).unwrap(), 0.5f64.sin());
    }

    #[test]
    fn fold_order_is_sin_ascending_then_cos_ascending() {
        // Behavioural contract 7: any reassociation changes the last ulp, so
        // this is a BIT-identical comparison against the hand-written fold in
        // old's order.
        let mut p = Prod::sin(0);
        p.mul_sin(0);
        p.mul_sin(0);
        p.mul_sin(1);
        p.mul_sin(1);
        p.mul_cos(0);
        let vals = [1.1, 2.1];
        let expected = 1.0f64 * 1.1f64.sin().powi(3) * 2.1f64.sin().powi(2) * 1.1f64.cos().powi(1);
        assert_eq!(p.eval(&vals).unwrap().to_bits(), expected.to_bits());
    }

    /// The cached fold must be **bit-for-bit** the uncached one: it is the fold
    /// order that behavioural contract 7 pins, and `sym_diff`'s
    /// `eval_fold_order_is_bit_identical_to_old` compares it against old.
    #[test]
    fn eval_cached_agrees_with_eval_bit_for_bit() {
        let vals = [0.3, 1.1, 2.1, -0.7, 4.2, 0.9, 3.3, 1.7];
        let mut cases = vec![Prod::new()];
        for (a, b) in [(0u32, 1u32), (1, 3), (2, 2), (0, 7), (5, 6)] {
            let mut p = Prod::sin(a);
            p.mul_sin(b);
            p.mul_sin(a);
            p.mul_cos(b);
            p.mul_cos(a);
            p.mul_cos(a);
            cases.push(p);
            // sin-only and cos-only, so both loop-skip guards are exercised.
            let mut s = Prod::sin(a);
            s.mul_sin(b);
            cases.push(s);
            let mut c = Prod::cos(a);
            c.mul_cos(b);
            cases.push(c);
        }
        for p in &cases {
            let direct = p.eval(&vals).unwrap();
            let cached = p.eval_cached(&mut AngleCache::new(&vals)).unwrap();
            assert_eq!(direct.to_bits(), cached.to_bits(), "{p} diverged");
        }
    }

    /// The cache must not pre-touch a variable the fold would not have reached:
    /// the reported variable is the first one the sin-then-cos walk trips on.
    #[test]
    fn cached_eval_reports_the_same_missing_variable_as_the_direct_fold() {
        let mut p = Prod::sin(7);
        p.mul_cos(3);
        let vals = [0.1, 0.2, 0.3, 0.4, 0.5];
        let direct = p.eval(&vals).unwrap_err().to_string();
        let cached = p
            .eval_cached(&mut AngleCache::new(&vals))
            .unwrap_err()
            .to_string();
        assert!(direct.contains("variable %7 not found"), "{direct}");
        assert_eq!(direct, cached);

        // The sum-level entry point carries the same message.
        let mut s = Sum::new();
        s.add_term(p, 1.0, usize::MAX, f64::EPSILON);
        assert!(
            s.eval(&vals)
                .unwrap_err()
                .to_string()
                .contains("variable %7 not found")
        );
    }

    /// Variables past the inline window fall back to recomputing per access —
    /// still correct, still the same value.
    #[test]
    fn variables_beyond_the_inline_window_still_evaluate() {
        let vals: Vec<f64> = (0..(INLINE_VARS + 8)).map(|i| i as f64 * 0.01).collect();
        let hi = (INLINE_VARS + 5) as u32;
        let mut p = Prod::sin(hi);
        p.mul_cos(1);
        let mut s = Sum::new();
        s.add_term(p.clone(), 2.0, usize::MAX, f64::EPSILON);
        assert_eq!(
            s.eval(&vals).unwrap().to_bits(),
            (0.0 + p.eval(&vals).unwrap() * 2.0).to_bits()
        );
    }

    #[test]
    fn eval_ignores_the_phase_but_eval_complex_does_not() {
        let mut p = Prod::sin(0);
        p.add_phase(1);
        let vals = [0.7];
        assert_eq!(p.eval(&vals).unwrap(), 0.7f64.sin());
        assert_eq!(
            p.eval_complex(&vals).unwrap(),
            Complex::new(0.0, 0.7f64.sin())
        );
    }
}
