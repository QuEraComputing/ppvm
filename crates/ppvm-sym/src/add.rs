use crate::term::{Item, Prod, Sum};

impl std::ops::AddAssign<f64> for Sum {
    fn add_assign(&mut self, rhs: f64) {
        self.c0 += rhs;
    }
}

impl std::ops::AddAssign<Item> for Sum {
    fn add_assign(&mut self, rhs: Item) {
        match rhs {
            Item::Sin(u) => {
                let p = Prod::sin(u);
                *self.terms.entry(p).or_insert(0.0) += 1.0;
            },
            Item::Cos(u) => {
                let p = Prod::cos(u);
                *self.terms.entry(p).or_insert(0.0) += 1.0;
            },
        }
    }
}
