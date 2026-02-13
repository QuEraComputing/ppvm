use num::complex::{Complex, ComplexFloat};
use num::traits::{One, Zero};

pub trait SparseVector<T, I>: Clone + IntoIterator<Item = (T, I)> {
    fn new() -> Self;
    /// Inserts an element without checking whether the index already exists.
    fn unsafe_insert(&mut self, index: I, value: T);
    fn add_or_insert(&mut self, index: I, value: T);
    fn get(&self, index: &I) -> T;
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool;
    fn mul_element_by(&mut self, index: I, factor: T);
    fn trim(&mut self, cutoff: T);
    fn retain(&mut self, f: impl FnMut(&(T, I)) -> bool);
    fn normalize(&mut self);
}

impl<T, I> SparseVector<T, I> for Vec<(T, I)>
where
    T: std::ops::AddAssign + std::ops::MulAssign + One + ComplexFloat + Zero,
    f64: std::ops::Div<<T as ComplexFloat>::Real>,
    I: std::cmp::PartialEq + Clone,
{
    fn new() -> Self {
        Vec::new()
    }

    fn unsafe_insert(&mut self, index: I, value: T) {
        self.push((value, index));
    }

    fn add_or_insert(&mut self, index: I, value: T) {
        for (v, i) in self.iter_mut() {
            if *i == index {
                *v += value;
                return;
            }
        }
        self.push((value, index));
    }

    fn get(&self, index: &I) -> T {
        for (v, i) in self.iter() {
            if i == index {
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

    fn mul_element_by(&mut self, index: I, factor: T) {
        for (v, i) in self.iter_mut() {
            if *i == index {
                *v *= factor;
                return;
            }
        }
    }

    fn trim(&mut self, cutoff: T) {
        // TODO: make cutoff real
        self.retain(|(element, _)| element.abs() > cutoff.abs());
    }

    fn retain(&mut self, f: impl FnMut(&(T, I)) -> bool) {
        Vec::retain(self, f);
    }

    fn normalize(&mut self) {
        let norm: T::Real = self
            .iter()
            .fold(T::Real::zero(), |acc, (v, _)| acc + v.abs() * v.abs());

        if norm == T::Real::zero() {
            panic!("Zero norm encountered during normalization");
        }
        if norm == T::Real::one() {
            return;
        }
        let norm_sqrt = norm.sqrt();
        let inv_norm_sqrt = T::Real::one() / norm_sqrt;
        for (v, _) in self.iter_mut() {
            // Scale by multiplying by 1/norm_sqrt
            let scale = Complex::new(inv_norm_sqrt, T::Real::zero());
            *v = *v * T::from(scale).expect("Failed to convert scale factor");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::complex::Complex64;

    #[test]
    fn test_sparse_vector_new() {
        let vec: Vec<(Complex64, usize)> = SparseVector::new();
        assert!(vec.is_empty());
        assert_eq!(vec.len(), 0);
    }

    #[test]
    fn test_unsafe_insert() {
        let mut vec: Vec<(Complex64, usize)> = SparseVector::new();
        vec.unsafe_insert(0, Complex64::new(1.0, 0.0));
        vec.unsafe_insert(2, Complex64::new(2.0, 1.0));

        assert_eq!(vec.len(), 2);
        assert_eq!(vec.get(&0), Complex64::new(1.0, 0.0));
        assert_eq!(vec.get(&2), Complex64::new(2.0, 1.0));
        assert_eq!(vec.get(&1), Complex64::new(0.0, 0.0)); // Non-existent index returns zero
    }

    #[test]
    fn test_add_or_insert() {
        let mut vec: Vec<(Complex64, usize)> = SparseVector::new();

        // First insert
        vec.add_or_insert(0, Complex64::new(1.0, 0.0));
        assert_eq!(vec.get(&0), Complex64::new(1.0, 0.0));

        // Add to existing
        vec.add_or_insert(0, Complex64::new(2.0, 1.0));
        assert_eq!(vec.get(&0), Complex64::new(3.0, 1.0));

        // Insert new
        vec.add_or_insert(5, Complex64::new(0.5, -0.5));
        assert_eq!(vec.get(&5), Complex64::new(0.5, -0.5));
    }

    #[test]
    fn test_mul_element_by() {
        let mut vec: Vec<(Complex64, usize)> = SparseVector::new();
        vec.unsafe_insert(0, Complex64::new(2.0, 1.0));
        vec.unsafe_insert(1, Complex64::new(1.0, 0.0));

        vec.mul_element_by(0, Complex64::new(2.0, 0.0));
        assert_eq!(vec.get(&0), Complex64::new(4.0, 2.0));
        assert_eq!(vec.get(&1), Complex64::new(1.0, 0.0)); // Unchanged

        // Multiply non-existent element (should do nothing)
        vec.mul_element_by(99, Complex64::new(5.0, 0.0));
        assert_eq!(vec.get(&99), Complex64::new(0.0, 0.0));
    }

    #[test]
    fn test_trim() {
        let mut vec: Vec<(Complex64, usize)> = SparseVector::new();
        vec.unsafe_insert(0, Complex64::new(1.0, 0.0));
        vec.unsafe_insert(1, Complex64::new(0.001, 0.0));
        vec.unsafe_insert(2, Complex64::new(0.0001, 0.0));
        vec.unsafe_insert(3, Complex64::new(2.0, 1.0));

        assert_eq!(vec.len(), 4);

        vec.trim(Complex64::new(0.01, 0.0));
        assert_eq!(vec.len(), 2); // Only elements with abs > 0.01 remain
        assert_eq!(vec.get(&0), Complex64::new(1.0, 0.0));
        assert_eq!(vec.get(&3), Complex64::new(2.0, 1.0));
        assert_eq!(vec.get(&1), Complex64::new(0.0, 0.0)); // Trimmed
    }

    #[test]
    fn test_retain() {
        let mut vec: Vec<(Complex64, usize)> = SparseVector::new();
        vec.unsafe_insert(0, Complex64::new(1.0, 0.0));
        vec.unsafe_insert(1, Complex64::new(2.0, 0.0));
        vec.unsafe_insert(2, Complex64::new(3.0, 0.0));
        vec.unsafe_insert(3, Complex64::new(4.0, 0.0));

        // Keep only even indices
        vec.retain(|(_, idx)| idx % 2 == 0);

        assert_eq!(vec.len(), 2);
        assert_eq!(vec.get(&0), Complex64::new(1.0, 0.0));
        assert_eq!(vec.get(&2), Complex64::new(3.0, 0.0));
        assert_eq!(vec.get(&1), Complex64::new(0.0, 0.0)); // Removed
    }

    #[test]
    fn test_normalize() {
        let mut vec: Vec<(Complex64, usize)> = SparseVector::new();
        vec.unsafe_insert(0, Complex64::new(3.0, 0.0));
        vec.unsafe_insert(1, Complex64::new(4.0, 0.0));

        // Norm should be sqrt(3^2 + 4^2) = 5
        vec.normalize();

        // After normalization, norm should be 1
        let norm: f64 = vec
            .iter()
            .map(|(v, _)| v.abs() * v.abs())
            .sum::<f64>()
            .sqrt();

        assert!(
            (norm - 1.0).abs() < 1e-10,
            "Norm should be 1 after normalization, got {}",
            norm
        );

        // Check values
        assert!((vec.get(&0).re - 0.6).abs() < 1e-10); // 3/5 = 0.6
        assert!((vec.get(&1).re - 0.8).abs() < 1e-10); // 4/5 = 0.8
    }

    #[test]
    fn test_normalize_complex() {
        let mut vec: Vec<(Complex64, usize)> = SparseVector::new();
        vec.unsafe_insert(0, Complex64::new(1.0, 1.0));
        vec.unsafe_insert(1, Complex64::new(1.0, -1.0));

        // |1+i| + |1-i| = sqrt(2) + sqrt(2) = 2*sqrt(2)
        // After normalization: 1/2 * sqrt(2)
        vec.normalize();

        let norm: f64 = vec
            .iter()
            .map(|(v, _)| v.abs() * v.abs())
            .sum::<f64>()
            .sqrt();

        assert!(
            (norm - 1.0).abs() < 1e-10,
            "Norm should be 1 after normalization"
        );
    }

    #[test]
    fn test_normalize_already_normalized() {
        let mut vec: Vec<(Complex64, usize)> = SparseVector::new();
        vec.unsafe_insert(0, Complex64::new(1.0, 0.0));

        vec.normalize(); // Should return early since norm is already 1

        assert_eq!(vec.get(&0), Complex64::new(1.0, 0.0));
    }
}
