// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::sum::PauliSum;
use num::Zero;
use ppvm_runtime::config::Config;
use ppvm_runtime::traits::{ACMapIter, Trace};

impl<'a, T: Config, Rhs> Trace<'a, Rhs> for PauliSum<T>
where
    <T as Config>::Coeff: Zero + Clone + std::ops::AddAssign + 'a,
    <T as Config>::Storage: 'a,
    <T as Config>::Map: Trace<'a, Rhs, Output = <T as Config>::Coeff>,
    <T as Config>::BuildHasher: 'a,
    <T as Config>::PauliWordType: 'a,
    Rhs: Trace<'a, T::PauliWordType, Output = bool> + 'a,
{
    type Output = T::Coeff;
    fn trace(&'a self, value: &'a Rhs) -> Self::Output {
        self.data().trace(value)
    }
}

impl<T: Config> PauliSum<T>
where
    for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    for<'a> T::Map: Trace<'a, T::PauliWordType, Output = T::Coeff>,
    T::Coeff: std::iter::Sum + Copy + std::ops::Mul<Output = T::Coeff>,
{
    /// Inner-product-style overlap with another sum: `Σ_k self[k] · other[k]`.
    pub fn overlap(&self, other: &Self) -> T::Coeff {
        other
            .data()
            .iter()
            .map(|(k, v)| *v * self.data().trace(k))
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_runtime::config::fxhash::ByteF64;

    fn sum(terms: &[(&str, f64)]) -> PauliSum<ByteF64<2>> {
        let mut s: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(4).build();
        for &(word, coeff) in terms {
            s += (word, coeff);
        }
        s
    }

    #[test]
    fn test_overlap_with_empty() {
        let a = sum(&[("IIII", 1.0), ("XIII", 2.0)]);
        let empty = sum(&[]);
        assert_eq!(a.overlap(&empty), 0.0);
        assert_eq!(empty.overlap(&a), 0.0);
    }

    #[test]
    fn test_overlap_single_term_with_itself() {
        let a = sum(&[("XIII", 3.0)]);
        assert_eq!(a.overlap(&a), 9.0);
    }

    #[test]
    fn test_overlap_orthogonal_paulis() {
        let a = sum(&[("XIII", 1.0)]);
        let b = sum(&[("YIII", 1.0)]);
        assert_eq!(a.overlap(&b), 0.0);
    }

    #[test]
    fn test_overlap_dot_product() {
        // overlap(a*I + b*X, c*I + d*X) = a*c + b*d
        let a = sum(&[("IIII", 1.0), ("XIII", 2.0)]);
        let b = sum(&[("IIII", 3.0), ("XIII", 4.0)]);
        assert_eq!(a.overlap(&b), 1.0 * 3.0 + 2.0 * 4.0);
    }

    #[test]
    fn test_overlap_partial_support() {
        // terms in `a` or `b` with no counterpart in the other contribute 0
        let a = sum(&[("IIII", 2.0), ("XIII", 3.0)]);
        let b = sum(&[("IIII", 5.0), ("YIII", 7.0)]);
        assert_eq!(a.overlap(&b), 2.0 * 5.0);
    }

    #[test]
    fn test_overlap_is_symmetric() {
        let a = sum(&[("IIII", 1.0), ("XIII", 2.0), ("YIII", 3.0)]);
        let b = sum(&[("IIII", 4.0), ("XIII", 5.0), ("ZIII", 6.0)]);
        assert_eq!(a.overlap(&b), b.overlap(&a));
    }
}
