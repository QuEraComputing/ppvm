// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `Display` for the symbolic types. A `Sum` prints its monomials in a total
//! order (sine degree, cosine degree, real before imaginary, then the factors),
//! so the output never depends on the table. An imaginary monomial prints a
//! leading `i`.

use std::fmt::Display;

use crate::prod::Prod;
use crate::sum::Sum;
use crate::term::{Inner, Term};

impl Display for Prod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let sin_len = self.factors.iter().filter(|x| x.sin > 0).count();
        let cos_len = self.factors.iter().filter(|x| x.cos > 0).count();
        if self.imag {
            write!(f, "i")?;
            if sin_len + cos_len > 0 {
                write!(f, " ")?;
            }
        }

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
            write!(f, "{:.3}", self.c0)?;
            if !self.is_empty() {
                write!(f, " + ")?;
            }
        }

        let mut sorted = self.iter().collect::<Vec<_>>();
        sorted.sort_by_cached_key(|(p, _)| (p.sin_pow(), p.cos_pow(), p.imag, p.factors.clone()));

        for (i, (p, c)) in sorted.iter().enumerate() {
            write!(f, "{c:.3} * {p}")?;
            if i + 1 < sorted.len() {
                write!(f, " + ")?;
            }
        }
        write!(f, "]")
    }
}

impl Display for Term {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.inner {
            Inner::Sum(ref s) => s.fmt(f),
            Inner::One(ref p, c) => write!(f, "{c:.3} * {p}"),
            Inner::Const(c) => write!(f, "{c}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::term::Var;

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
        s += Var(1).sin();
        s += Var(2).cos();
        assert_eq!(
            s.to_string(),
            "[3.000 + 1.000 * cos^1(%2) + 1.000 * sin^1(%1)]"
        );

        s.set_max_sin(2);
        s *= Var(2).sin();
        s *= Var(1).sin();
        assert_eq!(
            s.to_string(),
            "[3.000 * sin^1(%1) sin^1(%2) + 1.000 * sin^1(%1) sin^1(%2) cos^1(%2)]"
        );
    }

    #[test]
    fn display_of_the_small_forms() {
        // `Const` prints with NO precision spec, unlike every other form.
        assert_eq!(Term::from_f64(1.5).to_string(), "1.5");
        assert_eq!(Var(0).sin().to_string(), "1.000 * sin^1(%0)");
    }

    #[test]
    fn empty_sum_prints_as_brackets() {
        // The `mul_term` zero shortcut leaves an empty
        // `Sum`, which prints as `[]` — not `0`.
        let mut t = Term::from_f64(1.0) + Var(1).cos();
        t.set_max_sin(1);
        t *= Var(0).sin() * Var(1).sin();
        assert_eq!(t.to_string(), "[]");
        // A sum that never allocated its table renders the same way.
        assert_eq!(Sum::new().to_string(), "[]");
        // A constant-only sum closes its bracket.
        let mut s = Sum::new();
        s.c0 = 3.0;
        assert_eq!(s.to_string(), "[3.000]");
    }

    #[test]
    fn phase_and_sign_rendering() {
        let x = Var(0).sin();
        assert_eq!(x.mul_phase(1).to_string(), "1.000 * i sin^1(%0)");
        assert_eq!(x.mul_phase(3).to_string(), "-1.000 * i sin^1(%0)");
        assert_eq!(Term::from(2.0).mul_phase(1).to_string(), "2.000 * i");
        assert_eq!((x.clone() * -2.5).to_string(), "-2.500 * sin^1(%0)");
        let mut p = Prod::cos(0);
        p.mul_sin(0);
        assert_eq!(p.to_string(), "sin^1(%0) cos^1(%0)");
        assert_eq!(Prod::new().to_string(), "");
    }
}
