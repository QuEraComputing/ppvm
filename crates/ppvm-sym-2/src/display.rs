// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `Display` for the symbolic types.
//!
//! This output is **user-facing text and a snapshot contract**: the rendering,
//! the `{:.3}` precision, the separator placement and the term ordering are
//! ported byte-for-byte from `ppvm-sym/src/display.rs` (behavioural contract 8).
//! Two things keep it stable across runs:
//!
//! * the monomial table is keyed through a **seed-free** FxHash-class hasher (a
//!   `RandomState` would make the tie-break order run-dependent — integration
//!   baseline, perf feature 4); and
//! * the `Sum` sort key `(sin_pow, cos_pow)` is a **non-total** order, with ties
//!   broken by that deterministic iteration order. Imposing a total order here
//!   would be a visible divergence even though it is "nicer", so it is not done.
//!
//! The `phase` is **not** printed, matching old (`oldSuspectedBugs` #4).

use std::fmt::Display;

use crate::term::{Inner, Prod, Sum, Term};

impl Display for Prod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let sin_len = self.factors.iter().filter(|x| x.sin > 0).count();
        let cos_len = self.factors.iter().filter(|x| x.cos > 0).count();

        for (i, fac) in self.factors.iter().filter(|x| x.sin > 0).enumerate() {
            let (u, m) = (fac.var, fac.sin);
            write!(f, "sin^{m}(%{u})")?;
            if i + 1 < sin_len || cos_len > 0 {
                write!(f, " ")?;
            }
        }
        for (i, fac) in self.factors.iter().filter(|x| x.cos > 0).enumerate() {
            let (u, m) = (fac.var, fac.cos);
            write!(f, "cos^{m}(%{u})")?;
            if i + 1 < cos_len {
                write!(f, " ")?;
            }
        }
        Ok(())
    }
}

impl Display for Sum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[")?;
        if self.c0 != 0.0 {
            write!(f, "{:.3} ", self.c0)?;

            if self.is_empty() {
                return Ok(());
            } else {
                write!(f, "+ ")?;
            }
        }

        let Some(maps) = &self.maps else {
            return Ok(());
        };
        let mut sorted_keys = maps.terms.keys().collect::<Vec<_>>();
        sorted_keys.sort_by(|a, b| {
            a.sin_pow()
                .cmp(&b.sin_pow())
                .then(a.cos_pow().cmp(&b.cos_pow()))
        });

        for (i, p) in sorted_keys.iter().enumerate() {
            let c = maps.terms.get(*p).unwrap();
            write!(f, "{c:.3} * {p}")?;
            if i + 1 < sorted_keys.len() {
                write!(f, " + ")?;
            }
        }
        write!(f, "]")?;
        Ok(())
    }
}

impl Display for Term {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.inner {
            Inner::Sum(ref s) => s.fmt(f),
            Inner::One(ref p, c) => write!(f, "{c:.3} * {p}"),
            Inner::Var(u) => write!(f, "%{u}"),
            Inner::Const(c) => write!(f, "{c}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_display_prod() {
        let mut p = Prod::new();
        p.mul_sin(1);
        p.mul_sin(1);
        p.mul_cos(2);
        assert_eq!(p.to_string(), "sin^2(%1) cos^1(%2)");
    }

    #[test]
    fn test_display_term() {
        let mut s = Term::from_f64(3.0);
        s += Term::var(1).sin();
        s += Term::var(2).cos();
        assert_eq!(
            s.to_string(),
            "[3.000 + 1.000 * cos^1(%2) + 1.000 * sin^1(%1)]"
        );

        s.set_max_sin(2);
        s *= Term::var(2).sin();
        s *= Term::var(1).sin();
        assert_eq!(
            s.to_string(),
            "[3.000 * sin^1(%1) sin^1(%2) + 1.000 * sin^1(%1) sin^1(%2) cos^1(%2)]"
        );
    }

    #[test]
    fn display_of_the_small_forms() {
        assert_eq!(Term::var(7).to_string(), "%7");
        // `Const` prints with NO precision spec, unlike every other form.
        assert_eq!(Term::from_f64(1.5).to_string(), "1.5");
        assert_eq!(Term::var(0).sin().to_string(), "1.000 * sin^1(%0)");
    }

    #[test]
    fn empty_sum_prints_as_brackets() {
        // Behavioural contract 9: the `mul_term` zero shortcut leaves an empty
        // `Sum`, which prints as `[]` — not `0`.
        let mut t = Term::from_f64(1.0) + Term::var(1).cos();
        t.set_max_sin(1);
        t *= Term::var(0).sin() * Term::var(1).sin();
        assert_eq!(t.to_string(), "[]");
    }
}
