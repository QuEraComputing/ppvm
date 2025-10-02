use crate::term::{Inner, Prod, Sum, Term};

impl std::ops::AddAssign<f64> for Sum {
    fn add_assign(&mut self, rhs: f64) {
        self.c0 += rhs;
    }
}

impl std::ops::AddAssign<Prod> for Sum {
    fn add_assign(&mut self, rhs: Prod) {
        *self.terms.entry(rhs).or_insert(0.0) += 1.0;
    }
}

impl std::ops::AddAssign<f64> for Term {
    fn add_assign(&mut self, rhs: f64) {
        match self.inner {
            Inner::Const(ref mut c) => {
                *c += rhs;
            }
            Inner::One(ref p, c) => {
                let mut sum = Sum::new();
                sum.c0 = rhs;
                sum.add_term(p.clone(), c, self.max_sin, self.min_eps);
            }
            Inner::Sum(ref mut s) => {
                *s += rhs;
            }
            Inner::Var(_) => {
                panic!(
                    "variable is not used in sin/cos expressions,\
                    bare variable is not allowed in expression"
                );
            }
        }
    }
}

impl std::ops::AddAssign<Prod> for Term {
    fn add_assign(&mut self, rhs: Prod) {
        match self.inner {
            Inner::Const(c) => {
                let mut sum = Sum::new();
                sum.c0 = c;
                sum.add_term(rhs, 1.0, self.max_sin, self.min_eps);
                self.inner = Inner::Sum(sum);
            }
            Inner::One(ref p, c) => {
                let mut sum = Sum::new();
                sum.add_term(p.clone(), c, self.max_sin, self.min_eps);
                sum.add_term(rhs, 1.0, self.max_sin, self.min_eps);
                self.inner = Inner::Sum(sum);
            }
            Inner::Sum(ref mut s) => {
                *s += rhs;
            }
            Inner::Var(_) => {
                panic!(
                    "variable is not used in sin/cos expressions,\
                    bare variable is not allowed in expression"
                );
            }
        }
    }
}

impl std::ops::AddAssign<Term> for Term {
    fn add_assign(&mut self, rhs: Term) {
        match self.inner {
            Inner::Sum(ref mut s1) => match rhs.inner {
                Inner::Sum(s2) => {
                    for (p, c) in s2.terms {
                        s1.add_term(p, c, self.max_sin, self.min_eps);
                    }
                    s1.c0 += s2.c0;
                }
                Inner::One(p, c) => {
                    s1.add_term(p, c, self.max_sin, self.min_eps);
                }
                Inner::Const(c) => {
                    s1.add_const(c, self.min_eps);
                }
                _ => {
                    panic!(
                        "variable is not used in sin/cos expressions,\
                            bare variable is not allowed in expression"
                    );
                }
            },
            Inner::One(ref p, c) => match rhs.inner {
                Inner::Sum(s) => {
                    let mut new_sum = s;
                    new_sum.add_term(p.clone(), c, self.max_sin, self.min_eps);
                    self.inner = Inner::Sum(new_sum);
                }
                Inner::Const(c2) => {
                    let mut sum = Sum::new();
                    sum.c0 = c2;
                    sum.add_term(p.clone(), c, self.max_sin, self.min_eps);
                    self.inner = Inner::Sum(sum);
                }
                Inner::One(p2, c2) => {
                    let mut sum = Sum::new();
                    sum.add_term(p.clone(), c, self.max_sin, self.min_eps);
                    sum.add_term(p2, c2, self.max_sin, self.min_eps);
                    self.inner = Inner::Sum(sum);
                }
                _ => {
                    panic!(
                        "variable is not used in sin/cos expressions,\
                            bare variable is not allowed in expression"
                    );
                }
            },
            Inner::Const(c) => match rhs.inner {
                Inner::Sum(s) => {
                    self.inner = Inner::Sum(s);
                    *self += c;
                }
                Inner::One(p2, c2) => {
                    let mut sum = Sum::new();
                    sum.c0 = c;
                    sum.add_term(p2, c2, self.max_sin, self.min_eps);
                    self.inner = Inner::Sum(sum);
                }
                Inner::Const(c2) => {
                    self.inner = Inner::Const(c + c2);
                }
                _ => {
                    panic!(
                        "variable is not used in sin/cos expressions,\
                            bare variable is not allowed in expression"
                    );
                }
            },
            Inner::Var(_) => {
                panic!(
                    "variable is not used in sin/cos expressions,\
                    bare variable is not allowed in expression"
                );
            }
        }
    }
}

impl std::ops::Add for Term {
    type Output = Term;

    fn add(mut self, rhs: Term) -> Self::Output {
        self += rhs;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_sum() {
        let mut sum = Sum::new();
        sum += 1.0;
        sum += Prod::sin(0);
        assert_eq!(sum.eval(&[1.1]).unwrap(), 1.0 + 1.1f64.sin());
    }

    #[test]
    fn test_add_term() {
        let mut t = Term::from_f64(2.0);
        t += Prod::sin(0);
        t += Prod::cos(1);
        assert_eq!(
            t.eval(&[1.1, 2.2]).unwrap(),
            2.0 + 1.1f64.sin() + 2.2f64.cos()
        );

        let mut t = Term::var(0).sin();
        t += Term::from_f64(3.0);
        t += Term::var(1).cos();
        assert_eq!(
            t.eval(&[1.1, 2.2]).unwrap(),
            3.0 + 1.1f64.sin() + 2.2f64.cos()
        );
    }
}
