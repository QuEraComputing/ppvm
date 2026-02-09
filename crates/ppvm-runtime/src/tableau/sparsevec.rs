use crate::traits::Coefficient;

pub trait SparseVector<T>: Clone + IntoIterator<Item = (T, usize)> {
    fn new() -> Self;
    /// Inserts an element without checking whether the index already exists.
    fn unsafe_insert(&mut self, index: usize, value: T);
    fn add_or_insert(&mut self, index: usize, value: T);
    fn get(&self, index: usize) -> T;
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool;
    fn mul_element_by(&mut self, index: usize, factor: T);
}

impl<T: Coefficient> SparseVector<T> for Vec<(T, usize)> {
    fn new() -> Self {
        Vec::new()
    }

    fn unsafe_insert(&mut self, index: usize, value: T) {
        self.push((value, index));
    }

    fn add_or_insert(&mut self, index: usize, value: T) {
        for (v, i) in self.iter_mut() {
            if *i == index {
                *v += value;
                return;
            }
        }
        self.push((value, index));
    }

    fn get(&self, index: usize) -> T {
        for (v, i) in self.iter() {
            if *i == index {
                return v.clone();
            }
        }
        T::zero()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn mul_element_by(&mut self, index: usize, factor: T) {
        for (v, i) in self.iter_mut() {
            if *i == index {
                *v *= factor;
                return;
            }
        }
    }
}
