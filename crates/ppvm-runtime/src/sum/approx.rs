use super::data::PauliSum;
use crate::{
    config::Config,
    traits::{ACMapContains, ACMapIter},
    word::PauliWord,
};

impl<T: Config> approx::AbsDiffEq for PauliSum<T>
where
    T::Coeff: approx::AbsDiffEq,
    T::Map: PartialEq
        + for<'a> ACMapIter<'a, Item = (&'a PauliWord<T::Storage, T::BuildHasher>, &'a T::Coeff)>,
    <T::Coeff as approx::AbsDiffEq>::Epsilon: Copy,
{
    type Epsilon = <T::Coeff as approx::AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        <T::Coeff as approx::AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        if self.n_qubits() != other.n_qubits() {
            return false;
        }
        if self.len() != other.len() {
            return false;
        }
        for (k, v) in self.iter() {
            if !other
                .data()
                .contains_with(k, |ov| v.abs_diff_eq(ov, epsilon))
            {
                return false;
            }
        }
        true
    }
}

impl<T: Config> approx::RelativeEq for PauliSum<T>
where
    T::Coeff: approx::RelativeEq,
    T::Map: PartialEq
        + for<'a> ACMapIter<'a, Item = (&'a PauliWord<T::Storage, T::BuildHasher>, &'a T::Coeff)>,
    <T::Coeff as approx::AbsDiffEq>::Epsilon: Copy,
{
    fn default_max_relative() -> Self::Epsilon {
        <T::Coeff as approx::RelativeEq>::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        if self.n_qubits() != other.n_qubits() {
            return false;
        }
        if self.len() != other.len() {
            return false;
        }
        for (k, v) in self.iter() {
            if !other
                .data()
                .contains_with(k, |ov| v.relative_eq(ov, epsilon, max_relative))
            {
                return false;
            }
        }
        true
    }
}
