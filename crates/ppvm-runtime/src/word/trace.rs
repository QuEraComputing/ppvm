use std::hash::BuildHasher;

use super::data::PauliWord;
use crate::traits::{PauliStorage, Trace};

impl<'a, A, H, T> Trace<'a, PauliWord<A, H>, T> for PauliWord<A, H>
where
    A: PauliStorage + 'a,
    H: Default + BuildHasher + Clone + 'a,
    T: From<bool>,
{
    fn trace(&'a self, value: &'a PauliWord<A, H>) -> T {
        debug_assert_eq!(
            self.n_qubits(),
            value.n_qubits(),
            "#qubits mismatch, got {} and {}",
            self.n_qubits(),
            value.n_qubits()
        );
        T::from(self == value)
    }
}
