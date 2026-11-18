// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Evaluation at concrete angles. Within a monomial the fold multiplies the
//! `sin` factors, then the `cos` factors, each in ascending variable order;
//! that order fixes the floating-point result and is part of the contract.

use std::mem::MaybeUninit;

use num::Complex;

use crate::term::{Inner, Prod, Sum, Term};

/// Evaluation needed a variable `vals` has no value for: `vals.len() <= k`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MissingVariable(pub u32);

impl std::fmt::Display for MissingVariable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "variable %{} not found", self.0)
    }
}

impl std::error::Error for MissingVariable {}

/// Variables [`AngleCache`] memoizes; higher ids are recomputed per access.
const INLINE_VARS: usize = 32;

/// Per-call memo of `sin`/`cos` of each variable, so a [`Sum`] with many
/// monomials computes each transcendental once. Bit-identical to the uncached
/// fold, and bounds-checked per access so the same [`MissingVariable`] is
/// reported first. Slots stay uninitialized until their validity bit is set.
struct AngleCache<'a> {
    vals: &'a [f64],
    /// Bit `i` set <=> the corresponding slot is initialized.
    sin_valid: u32,
    cos_valid: u32,
    sin: [MaybeUninit<f64>; INLINE_VARS],
    cos: [MaybeUninit<f64>; INLINE_VARS],
}

impl<'a> AngleCache<'a> {
    #[inline]
    fn new(vals: &'a [f64]) -> Self {
        Self {
            vals,
            sin_valid: 0,
            cos_valid: 0,
            sin: [MaybeUninit::uninit(); INLINE_VARS],
            cos: [MaybeUninit::uninit(); INLINE_VARS],
        }
    }

    /// `sin` or `cos` of `vals[k]`, memoized for `k < INLINE_VARS`.
    #[inline]
    fn get(&mut self, k: u32, trig: Trig) -> Result<f64, MissingVariable> {
        let i = k as usize;
        let v = *self.vals.get(i).ok_or(MissingVariable(k))?;
        let f = match trig {
            Trig::Sin => f64::sin,
            Trig::Cos => f64::cos,
        };
        if i >= INLINE_VARS {
            return Ok(f(v));
        }
        let (valid, slots) = match trig {
            Trig::Sin => (&mut self.sin_valid, &mut self.sin),
            Trig::Cos => (&mut self.cos_valid, &mut self.cos),
        };
        if *valid & (1 << i) == 0 {
            slots[i].write(f(v));
            *valid |= 1 << i;
        }
        // SAFETY: the validity bit is set only after this slot is initialized.
        Ok(unsafe { slots[i].assume_init() })
    }
}

/// Which trigonometric function a fold step needs.
#[derive(Clone, Copy)]
enum Trig {
    Sin,
    Cos,
}

impl Prod {
    /// The real part of [`Prod::eval_complex`] at `vals` (`vals[i]` binds
    /// variable `i`): `0.0` for an imaginary monomial, `1.0` for the empty one.
    #[cfg(test)]
    #[inline]
    pub(crate) fn eval(&self, vals: &[f64]) -> Result<f64, MissingVariable> {
        Ok(self.eval_complex(vals)?.re)
    }

    /// The fold in contract order; every evaluation path goes through it.
    #[inline]
    fn fold(
        &self,
        mut trig_of: impl FnMut(u32, Trig) -> Result<f64, MissingVariable>,
    ) -> Result<f64, MissingVariable> {
        let mut res = 1.0;
        for f in self.factors.iter().filter(|f| f.sin > 0) {
            res *= trig_of(f.var, Trig::Sin)?.powi(f.sin as i32);
        }
        for f in self.factors.iter().filter(|f| f.cos > 0) {
            res *= trig_of(f.var, Trig::Cos)?.powi(f.cos as i32);
        }
        Ok(res)
    }

    /// The value at `vals`, including the phase.
    pub(crate) fn eval_complex(&self, vals: &[f64]) -> Result<Complex<f64>, MissingVariable> {
        Ok(self.with_phase(self.fold_uncached(vals)?))
    }

    /// [`Prod::fold`] reading `vals` directly, for a one-off monomial.
    #[inline]
    fn fold_uncached(&self, vals: &[f64]) -> Result<f64, MissingVariable> {
        self.fold(|k, trig| {
            let v = vals.get(k as usize).ok_or(MissingVariable(k))?;
            Ok(match trig {
                Trig::Sin => v.sin(),
                Trig::Cos => v.cos(),
            })
        })
    }

    /// `v`, or `i·v` for an imaginary monomial.
    #[inline]
    fn with_phase(&self, v: f64) -> Complex<f64> {
        if self.imag {
            Complex::new(0.0, v)
        } else {
            Complex::new(v, 0.0)
        }
    }
}

impl Sum {
    /// Phase-aware evaluation; see [`Prod::eval_complex`].
    pub(crate) fn eval_complex(&self, vals: &[f64]) -> Result<Complex<f64>, MissingVariable> {
        let mut cache = AngleCache::new(vals);
        let mut res = Complex::new(self.c0, 0.0);
        for (p, c) in self.iter() {
            res += p.with_phase(p.fold(|k, trig| cache.get(k, trig))?) * c;
        }
        Ok(res)
    }
}

impl Term {
    /// The real part of [`Term::eval_complex`].
    #[inline]
    pub fn eval(&self, vals: &[f64]) -> Result<f64, MissingVariable> {
        Ok(self.eval_complex(vals)?.re)
    }

    /// The value at `vals` (`vals[i]` binds variable `i`), including the phase.
    pub fn eval_complex(&self, vals: &[f64]) -> Result<Complex<f64>, MissingVariable> {
        match self.inner {
            Inner::Const(c) => Ok(Complex::new(c, 0.0)),
            Inner::One(ref p, c) => Ok(p.eval_complex(vals)? * c),
            Inner::Sum(ref s) => s.eval_complex(vals),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::term::Var;
    use crate::term::{DEFAULT_MAX_SIN, DEFAULT_MIN_EPS};

    #[test]
    fn missing_variable_is_a_typed_error() {
        let err = Var(5).sin().eval(&[0.1]).unwrap_err();
        assert_eq!(err, MissingVariable(5));
        assert_eq!(err.to_string(), "variable %5 not found");

        // The cached `Sum` path reports the first missing variable the
        // sin-then-cos walk reaches (`%7`), not the lowest missing id (`%5`).
        let mut p = Prod::sin(7);
        p.mul_cos(5);
        let mut s = Sum::new();
        s.add_term(p, 1.0, DEFAULT_MAX_SIN, DEFAULT_MIN_EPS);
        assert_eq!(
            s.eval_complex(&[0.1, 0.2, 0.3, 0.4]).map(|z| z.re),
            Err(MissingVariable(7))
        );
    }

    #[test]
    fn fold_order_is_sin_ascending_then_cos_ascending() {
        // Bit-identical to a hand-written fold in contract order.
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

    /// table order, across the inline window and with repeated variables.
    #[test]
    fn cached_sum_eval_is_bit_identical_to_uncached() {
        let n = INLINE_VARS as u32 + 4;
        let vals: Vec<f64> = (0..n).map(|i| 0.37 * i as f64 - 3.0).collect();
        let mut s = Sum::new();
        s.c0 = 0.125;
        for i in 0..n {
            let mut p = Prod::sin(i);
            p.mul_cos((i * 7) % n);
            p.mul_sin((i * 3) % n);
            if i % 3 == 0 {
                p.mul_phase(1);
            }
            s.add_term(p, 1.0 + i as f64, DEFAULT_MAX_SIN, DEFAULT_MIN_EPS);
        }
        let mut re = s.c0;
        let mut z = Complex::new(s.c0, 0.0);
        for (p, c) in s.iter() {
            re += p.eval(&vals).unwrap() * c;
            z += p.eval_complex(&vals).unwrap() * c;
        }
        assert_eq!(s.eval_complex(&vals).unwrap().re.to_bits(), re.to_bits());
        let got = s.eval_complex(&vals).unwrap();
        assert_eq!(
            (got.re.to_bits(), got.im.to_bits()),
            (z.re.to_bits(), z.im.to_bits())
        );
    }

    #[test]
    fn small_forms_and_the_real_part() {
        assert_eq!(Prod::new().eval(&[]), Ok(1.0));
        assert_eq!(Sum::new().eval_complex(&[]), Ok(Complex::new(0.0, 0.0)));
        let mut p = Prod::cos(0);
        p.mul_phase(1);
        assert_eq!(p.eval(&[0.0]), Ok(0.0));
        assert_eq!(p.eval_complex(&[0.0]), Ok(Complex::new(0.0, 1.0)));
    }
}
