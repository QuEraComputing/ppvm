use fxhash::FxHashMap;

use crate::term::{Item, Prod, Sum};

impl Prod {
    pub fn mul_sin(&mut self, u: u32) {
        *self.sin.entry(u).or_insert(0) += 1;
    }

    pub fn mul_cos(&mut self, u: u32) {
        *self.cos.entry(u).or_insert(0) += 1;
    }
}

impl std::ops::MulAssign<f64> for Sum {
    fn mul_assign(&mut self, rhs: f64) {
        self.c0 *= rhs;
        for v in self.terms.values_mut() {
            *v *= rhs;
        }
    }
}

impl std::ops::MulAssign<Item> for Sum {
    fn mul_assign(&mut self, rhs: Item) {
        use Item::*;
        match rhs {
            Sin(u) => {
                let mut new_terms = FxHashMap::default();

                if self.c0.abs() > 1e-12 {
                    *new_terms.entry(Prod::sin(u)).or_insert(0.0) += self.c0;
                    self.c0 = 0.0;
                }

                for (mut p, coeff) in self.terms.drain() {
                    if self.max < p.sin_pow() + 1 {
                        continue;
                    }
                    p.mul_sin(u);
                    *new_terms.entry(p).or_insert(0.0) += coeff;
                }
                self.terms = new_terms;
            }
            Cos(u) => {
                let mut new_terms = FxHashMap::default();
                if self.c0.abs() > 1e-12 {
                    *new_terms.entry(Prod::cos(u)).or_insert(0.0) += self.c0;
                    self.c0 = 0.0;
                }
                for (mut p, coeff) in self.terms.drain() {
                    p.mul_cos(u);
                    *new_terms.entry(p).or_insert(0.0) += coeff;
                }
                self.terms = new_terms;
            }
        }
    }
}
