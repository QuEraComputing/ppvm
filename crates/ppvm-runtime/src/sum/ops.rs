use std::ops::{AddAssign, MulAssign};

use num::One;

use crate::traits::*;
use crate::word::PauliWord;
use crate::{config::Config, sum::PauliSum};

impl<T: Config> MulAssign<T::Coeff> for PauliSum<T>
where
    T::Map: ACMapMulAssign<T::Coeff>,
    T::Coeff: MulAssign + Clone,
{
    fn mul_assign(&mut self, rhs: T::Coeff) {
        self.data_mut().mul_assign(rhs.clone());
    }
}

impl<T: Config> AddAssign<PauliSum<T>> for PauliSum<T>
where
    T::Coeff: std::ops::AddAssign,
    T::Map: Extend<(PauliWord<T::Storage>, T::Coeff)>
        + IntoIterator<Item = (PauliWord<T::Storage>, T::Coeff)>,
{
    fn add_assign(&mut self, rhs: PauliSum<T>) {
        debug_assert_eq!(self.n_qubits(), rhs.n_qubits());
        debug_assert!(self.len() + rhs.len() <= self.capacity());
        self.extend(rhs.into_iter())
    }
}

impl<'a, T: Config> AddAssign<&'a PauliSum<T>> for PauliSum<T>
where
    T::Coeff: std::ops::AddAssign,
    T::Map: Extend<(PauliWord<T::Storage>, T::Coeff)>
        + ACMapIter<'a, Item = (PauliWord<T::Storage>, T::Coeff)>,
{
    fn add_assign(&mut self, rhs: &'a PauliSum<T>) {
        debug_assert_eq!(self.n_qubits(), rhs.n_qubits());
        debug_assert!(self.len() + rhs.len() <= self.capacity());
        self.extend(rhs.data().iter().map(|(k, v)| (k.clone(), v.clone())))
    }
}

impl<T: Config, P> AddAssign<(P, T::Coeff)> for PauliSum<T>
where
    P: Into<PauliWord<T::Storage>>,
    T::Coeff: std::ops::AddAssign,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff>,
{
    fn add_assign(&mut self, rhs: (P, T::Coeff)) {
        let key = rhs.0.into();
        debug_assert_eq!(self.n_qubits(), key.n_qubits());
        debug_assert!(self.len() + 1 <= self.capacity());
        ACMapAddAssign::add_assign(self.data_mut(), key, rhs.1);
    }
}

impl<T: Config, P> AddAssign<P> for PauliSum<T>
where
    P: Into<PauliWord<T::Storage>>,
    T::Coeff: std::ops::AddAssign + One,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff>,
{
    fn add_assign(&mut self, rhs: P) {
        let key = rhs.into();
        debug_assert_eq!(self.n_qubits(), key.n_qubits());
        debug_assert!(self.len() + 1 <= self.capacity());
        ACMapAddAssign::add_assign(self.data_mut(), key, T::Coeff::one());
    }
}
