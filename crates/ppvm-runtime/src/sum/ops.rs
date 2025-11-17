use std::ops::AddAssign;

use num::One;

use crate::phase::PhasedPauliWord;
use crate::traits::*;
use crate::word::PauliWord;
use crate::{config::Config, sum::PauliSum};

/// Implement multiplication by a scalar coefficient for PauliSum.
/// Use this macro to implement `MulAssign` for other scalar types
/// if needed. This macro will forward the multiplication to the underlying map's
/// `mul_assign` method.
#[macro_export]
macro_rules! impl_op_mul_assign_coefficient {
    ($ty:ty) => {
        impl<T: Config> std::ops::MulAssign<$ty> for PauliSum<T>
        where
            T::Map: ACMapMulAssign<$ty, T::BuildHasher>,
        {
            fn mul_assign(&mut self, rhs: $ty) {
                self.data_mut().mul_assign(rhs.clone());
            }
        }
    };
}

pub use impl_op_mul_assign_coefficient;

impl_op_mul_assign_coefficient!(f64);

impl<T: Config> std::ops::Mul<PauliSum<T>> for PauliSum<T>
where
    T::BuildHasher: Sync + Send,
    T::Coeff: ComplexCoefficient,
    T::Map: for<'a> ACMapIter<'a, Item = (&'a PauliWord<T::Storage, T::BuildHasher>, &'a T::Coeff)>,
{
    type Output = PauliSum<T>;

    fn mul(self, rhs: PauliSum<T>) -> Self::Output {
        let mut output = self.clone();
        output *= rhs;
        output
    }
}

impl<T: Config> std::ops::Mul<PauliWord<T::Storage, T::BuildHasher>> for PauliSum<T>
where
    T::BuildHasher: Sync + Send,
    T::Coeff: ComplexCoefficient,
{
    type Output = PauliSum<T>;

    fn mul(self, rhs: PauliWord<T::Storage, T::BuildHasher>) -> Self::Output {
        let mut output = self.clone();
        output *= rhs;
        output
    }
}

impl<T: Config> std::ops::MulAssign<PauliSum<T>> for PauliSum<T>
where
    T::BuildHasher: Sync + Send,
    T::Coeff: ComplexCoefficient,
    T::Map: for<'a> ACMapIter<'a, Item = (&'a PauliWord<T::Storage, T::BuildHasher>, &'a T::Coeff)>,
{
    fn mul_assign(&mut self, rhs: PauliSum<T>) {
        for (rhs_word, rhs_coeff) in rhs.iter() {
            let phased_rhs = PhasedPauliWord {
                word: rhs_word.clone(),
                phase: 0,
            };
            self.map_add(|word, coeff| {
                let mut phased_word = PhasedPauliWord {
                    word: word.clone(),
                    phase: 0,
                };
                phased_word *= phased_rhs.clone();
                let new_coeff = coeff.mul_phase(phased_word.phase);
                (phased_word.word, new_coeff * rhs_coeff.clone())
            });
        }
    }
}

impl<T: Config> std::ops::MulAssign<PauliWord<T::Storage, T::BuildHasher>> for PauliSum<T>
where
    T::BuildHasher: Sync + Send,
    T::Coeff: ComplexCoefficient,
{
    fn mul_assign(&mut self, rhs: PauliWord<T::Storage, T::BuildHasher>) {
        let phased_rhs = PhasedPauliWord {
            word: rhs,
            phase: 0,
        };
        self.map_add(|word, coeff| {
            let mut phased_word = PhasedPauliWord {
                word: word.clone(),
                phase: 0,
            };
            phased_word *= phased_rhs.clone();
            (phased_word.word, coeff.mul_phase(phased_word.phase))
        });
    }
}

impl<T: Config> AddAssign<PauliSum<T>> for PauliSum<T>
where
    T::Coeff: std::ops::AddAssign,
    T::Map: Extend<(PauliWord<T::Storage, T::BuildHasher>, T::Coeff)>
        + IntoIterator<Item = (PauliWord<T::Storage, T::BuildHasher>, T::Coeff)>,
{
    fn add_assign(&mut self, rhs: PauliSum<T>) {
        debug_assert_eq!(self.n_qubits(), rhs.n_qubits());
        self.extend(rhs.into_iter())
    }
}

impl<'a, T: Config> AddAssign<&'a PauliSum<T>> for PauliSum<T>
where
    T::Coeff: std::ops::AddAssign,
    T::Map: Extend<(PauliWord<T::Storage, T::BuildHasher>, T::Coeff)>
        + ACMapIter<'a, Item = (PauliWord<T::Storage, T::BuildHasher>, T::Coeff)>,
{
    fn add_assign(&mut self, rhs: &'a PauliSum<T>) {
        debug_assert_eq!(self.n_qubits(), rhs.n_qubits());
        self.extend(rhs.data().iter().map(|(k, v)| (k.clone(), v.clone())))
    }
}

impl<T: Config, P> AddAssign<(P, T::Coeff)> for PauliSum<T>
where
    P: Into<PauliWord<T::Storage, T::BuildHasher>>,
    T::Coeff: std::ops::AddAssign,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher>,
{
    fn add_assign(&mut self, rhs: (P, T::Coeff)) {
        let key = rhs.0.into();
        debug_assert_eq!(self.n_qubits(), key.n_qubits());
        ACMapAddAssign::add_assign(self.data_mut(), key, rhs.1);
    }
}

impl<T: Config, P> AddAssign<P> for PauliSum<T>
where
    P: Into<PauliWord<T::Storage, T::BuildHasher>>,
    T::Coeff: std::ops::AddAssign + One,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher>,
{
    fn add_assign(&mut self, rhs: P) {
        let key = rhs.into();
        debug_assert_eq!(self.n_qubits(), key.n_qubits());
        ACMapAddAssign::add_assign(self.data_mut(), key, T::Coeff::one());
    }
}
