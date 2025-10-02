use std::fmt::Display;

use crate::term::Prod;

impl Display for Prod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, (u, m)) in self.sin.iter().enumerate() {
            write!(f, "sin^{m}(%{u})")?;
            if i + 1 < self.sin.len() || !self.cos.is_empty() {
                write!(f, " ")?;
            }
        }
        for (i, (u, m)) in self.cos.iter().enumerate() {
            write!(f, "cos^{m}(%{u})")?;
            if i + 1 < self.cos.len() {
                write!(f, " ")?;
            }
        }
        Ok(())
    }
}

impl Display for crate::term::Sum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.c0 != 0.0 {
            write!(f, "{:.3} ", self.c0)?;

            if self.terms.is_empty() {
                return Ok(());
            } else {
                write!(f, "+ ")?;
            }
        }

        let mut sorted_keys = self.terms.keys().collect::<Vec<_>>();
        sorted_keys.sort_by(|a, b| {
            a.sin_pow()
                .cmp(&b.sin_pow())
                .then(a.cos_pow().cmp(&b.cos_pow()))
        });

        for (i, p) in sorted_keys.iter().enumerate() {
            let c = self.terms.get(p).unwrap();
            write!(f, "{:.3} * {}", c, p)?;
            if i + 1 < sorted_keys.len() {
                write!(f, " + ")?;
            }
        }
        Ok(())
    }
}

impl Display for crate::term::Term {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use crate::term::Inner::*;
        match self.inner {
            Sum(ref s) => s.fmt(f),
            One(ref p, c) => write!(f, "{:.3} * {}", c, p),
            Var(u) => write!(f, "%{u}"),
            Const(c) => write!(f, "{c}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::Term;

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
            "3.000 + 1.000 * cos^1(%2) + 1.000 * sin^1(%1)"
        );

        s.set_max_sin(2);
        s *= Term::var(2).sin();
        s *= Term::var(1).sin();
        assert_eq!(
            s.to_string(),
            "3.000 * sin^1(%1) sin^1(%2) + 1.000 * sin^1(%1) sin^1(%2) cos^1(%2)"
        );
    }
}
