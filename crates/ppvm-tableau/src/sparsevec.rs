// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::complex::{Complex, ComplexFloat};
use num::traits::{One, Zero};

/// A sparse vector keyed by index type `I` with values of type `T`.
///
/// Used by [`GeneralizedTableau`](crate::data::GeneralizedTableau) to
/// store the coefficient over each branching bitstring. The default
/// implementation is `Vec<(T, I)>`; alternative backings (BTreeMap,
/// HashMap, etc.) can be provided by downstream code.
pub trait SparseVector<T, I>: Clone + IntoIterator<Item = (T, I)> {
    /// Construct an empty sparse vector.
    fn new() -> Self;
    /// Inserts an element without checking whether the index already exists.
    fn unsafe_insert(&mut self, index: I, value: T);
    /// Add `value` into the entry at `index`, creating it if absent.
    fn add_or_insert(&mut self, index: I, value: T);
    /// Retrieve the value at `index`, or zero if absent.
    fn get(&self, index: &I) -> T;
    /// Number of stored entries.
    fn len(&self) -> usize;
    /// `true` if no entries are stored.
    fn is_empty(&self) -> bool;
    /// Multiply every entry's value by `factor`.
    fn mul_by(&mut self, factor: T);
    /// Multiply the value at `index` by `factor`. No-op if absent.
    fn mul_element_by(&mut self, index: I, factor: T);
    /// Drop entries whose magnitude is at most `|cutoff|`.
    fn trim(&mut self, cutoff: T);
    /// Drop entries failing the predicate `f`.
    fn retain(&mut self, f: impl FnMut(&(T, I)) -> bool);
    /// L2-normalize the vector in place. Panics on zero norm.
    fn normalize(&mut self);
    /// Reserve capacity for at least `additional` more entries. Backings
    /// that don't support pre-allocation can leave this as a no-op.
    fn reserve(&mut self, _additional: usize) {}
    /// Borrow the stored entries as an iterator without consuming the vector.
    fn iter<'a>(&'a self) -> impl Iterator<Item = &'a (T, I)>
    where
        T: 'a,
        I: 'a;
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
                return *v;
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

    fn mul_by(&mut self, factor: T) {
        for (v, _) in self.iter_mut() {
            *v *= factor;
        }
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
        let c_re = cutoff.re();
        let c_im = cutoff.im();
        let cutoff_sq = c_re * c_re + c_im * c_im;
        self.retain(|(element, _)| {
            let e_re = element.re();
            let e_im = element.im();
            e_re * e_re + e_im * e_im > cutoff_sq
        });
    }

    fn retain(&mut self, f: impl FnMut(&(T, I)) -> bool) {
        Vec::retain(self, f);
    }

    fn reserve(&mut self, additional: usize) {
        Vec::reserve(self, additional);
    }

    fn iter<'a>(&'a self) -> impl Iterator<Item = &'a (T, I)>
    where
        T: 'a,
        I: 'a,
    {
        <[(T, I)]>::iter(self)
    }

    fn normalize(&mut self) {
        // `re*re + im*im` directly; `abs() * abs()` would compute
        // `hypot(re, im)` (a sqrt) only to immediately square it again.
        let norm: T::Real = self.iter().fold(T::Real::zero(), |acc, (v, _)| {
            let re = v.re();
            let im = v.im();
            acc + re * re + im * im
        });

        if norm == T::Real::zero() {
            panic!("Zero norm encountered during normalization");
        }
        if norm == T::Real::one() {
            return;
        }
        let norm_sqrt = norm.sqrt();
        let inv_norm_sqrt = T::Real::one() / norm_sqrt;
        let scale: T = T::from(Complex::new(inv_norm_sqrt, T::Real::zero()))
            .expect("Failed to convert scale factor");
        for (v, _) in self.iter_mut() {
            *v *= scale;
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

    use bnum::types::U256;
    #[test]
    fn test_bigint() {
        let mut vec: Vec<(Complex64, U256)> = SparseVector::new();

        vec.unsafe_insert(U256::from(0u8), Complex64::new(1.0, 0.0));

        vec.normalize(); // Should return early since norm is already 1

        assert_eq!(vec.get(&U256::from(0u8)), Complex64::new(1.0, 0.0));
    }
}
