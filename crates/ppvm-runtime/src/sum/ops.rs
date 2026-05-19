// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::ops::{AddAssign, MulAssign};

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
    T::Map: for<'a> ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>: for<'a> From<&'a T::PauliWordType>
        + MulAssign<PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>>
        + Send
        + Sync,
{
    type Output = PauliSum<T>;

    fn mul(self, rhs: PauliSum<T>) -> Self::Output {
        let mut output = self.clone();
        output *= rhs;
        output
    }
}

impl<T: Config, W: PauliWordTrait> std::ops::Mul<W> for PauliSum<T>
where
    T::BuildHasher: Sync + Send,
    T::Coeff: ComplexCoefficient,
    PauliSum<T>: std::ops::MulAssign<W>,
{
    type Output = PauliSum<T>;

    fn mul(self, rhs: W) -> Self::Output {
        let mut output = self.clone();
        output *= rhs;
        output
    }
}

impl<T: Config> std::ops::MulAssign<PauliSum<T>> for PauliSum<T>
where
    T::BuildHasher: Sync + Send,
    T::Coeff: ComplexCoefficient,
    T::Map: for<'a> ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>: for<'a> From<&'a T::PauliWordType>
        + MulAssign<PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>>
        + Send
        + Sync,
{
    fn mul_assign(&mut self, rhs: PauliSum<T>) {
        for (rhs_word, rhs_coeff) in rhs.iter() {
            let phased_rhs: PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType> =
                rhs_word.into();
            self.map_add(|word, coeff| {
                let mut phased_word: PhasedPauliWord<_, _, <T as Config>::PauliWordType> =
                    word.into();
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
    T::PauliWordType: From<PauliWord<T::Storage, T::BuildHasher>>,
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>:
        MulAssign + From<T::PauliWordType> + Send + Sync,
{
    fn mul_assign(&mut self, rhs: PauliWord<T::Storage, T::BuildHasher>) {
        let rhs_word: T::PauliWordType = rhs.into();
        let phased_rhs: PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType> =
            rhs_word.into();
        self.map_add(|word, coeff| {
            let mut phased_word: PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType> =
                word.clone().into();
            phased_word *= phased_rhs.clone();
            (phased_word.word, coeff.mul_phase(phased_word.phase))
        });
    }
}

impl<T: Config> AddAssign<PauliSum<T>> for PauliSum<T>
where
    T::Coeff: std::ops::AddAssign,
    T::Map:
        Extend<(T::PauliWordType, T::Coeff)> + IntoIterator<Item = (T::PauliWordType, T::Coeff)>,
{
    fn add_assign(&mut self, rhs: PauliSum<T>) {
        debug_assert_eq!(self.n_qubits(), rhs.n_qubits());
        self.extend(rhs)
    }
}

impl<'a, T: Config> AddAssign<&'a PauliSum<T>> for PauliSum<T>
where
    T::Coeff: std::ops::AddAssign,
    T::Map:
        Extend<(T::PauliWordType, T::Coeff)> + ACMapIter<'a, Item = (T::PauliWordType, T::Coeff)>,
{
    fn add_assign(&mut self, rhs: &'a PauliSum<T>) {
        debug_assert_eq!(self.n_qubits(), rhs.n_qubits());
        self.extend(rhs.data().iter().map(|(k, v)| (k.clone(), v.clone())))
    }
}

impl<T: Config, P> AddAssign<(P, T::Coeff)> for PauliSum<T>
where
    P: Into<T::PauliWordType>,
    T::Coeff: std::ops::AddAssign,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>,
{
    fn add_assign(&mut self, rhs: (P, T::Coeff)) {
        let key = rhs.0.into();
        debug_assert_eq!(self.n_qubits(), key.n_qubits());
        ACMapAddAssign::add_assign(self.data_mut(), key, rhs.1);
    }
}

// NOTE: to avoid conflicts with the (P, T::Coeff) impl above, we implement a custom trait
// which is not implemented for tuples

mod private_into_pauli_word {
    // NOTE: sealed trait pattern: so downstream crates can't implement the IntoPauliWord trait
    // we additionally require this private trait to be implemented; otherwise, we again get
    // conflicting implementations
    pub trait SealedIntoPauliWord {}
}

trait IntoPauliWord<T: Config>:
    Into<T::PauliWordType> + private_into_pauli_word::SealedIntoPauliWord
{
}
impl<S: PauliStorage, H: std::hash::BuildHasher + Default + Clone>
    private_into_pauli_word::SealedIntoPauliWord for PauliWord<S, H>
{
}
impl private_into_pauli_word::SealedIntoPauliWord for &str {}
impl private_into_pauli_word::SealedIntoPauliWord for String {}
impl<T: Config, P> IntoPauliWord<T> for P where
    P: Into<T::PauliWordType> + private_into_pauli_word::SealedIntoPauliWord
{
}

impl<T: Config, P> AddAssign<P> for PauliSum<T>
where
    P: IntoPauliWord<T>,
    T::Coeff: std::ops::AddAssign + One,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>,
{
    fn add_assign(&mut self, rhs: P) {
        let key = rhs.into();
        debug_assert_eq!(self.n_qubits(), key.n_qubits());
        ACMapAddAssign::add_assign(self.data_mut(), key, T::Coeff::one());
    }
}
