use super::term::{Prod, Sum, Term};
use anyhow::Result;

impl Prod {
    pub fn eval(&self, vals: &[f64]) -> Result<f64> {
        if self.pow() == 0 {
            return Ok(1.0);
        }

        let mut res = 1.0;
        for (k, v) in &self.sin {
            res *= vals
                .get(*k as usize)
                .ok_or_else(|| anyhow::anyhow!("variable %{k} not found"))?
                .sin()
                .powi(*v as i32);
        }

        for (k, v) in &self.cos {
            res *= vals
                .get(*k as usize)
                .ok_or_else(|| anyhow::anyhow!("variable %{k} not found"))?
                .cos()
                .powi(*v as i32);
        }
        Ok(res)
    }
}

impl Sum {
    pub fn eval(&self, vals: &[f64]) -> Result<f64> {
        let mut res = self.c0;
        for (p, c) in &self.terms {
            res += p.eval(vals)? * c;
        }
        Ok(res)
    }
}

impl Term {
    pub fn eval(&self, vals: &[f64]) -> Result<f64> {
        use crate::term::Inner::*;
        match self.inner {
            Const(c) => Ok(c),
            Var(u) => vals
                .get(u as usize)
                .copied()
                .ok_or_else(|| anyhow::anyhow!("variable %{u} not found")),
            One(ref p, c) => Ok(p.eval(vals)? * c),
            Sum(ref s) => s.eval(vals),
        }
    }
}
