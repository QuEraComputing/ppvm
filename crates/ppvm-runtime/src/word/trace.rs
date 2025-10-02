use std::hash::BuildHasher;

use super::data::PauliWord;
use crate::traits::{PauliStorage, Trace};

impl<'a, A, H> Trace<'a, PauliWord<A, H>> for PauliWord<A, H>
where
    A: PauliStorage + 'a,
    H: Default + BuildHasher + Clone + 'a,
{
    type Output = bool;
    fn trace(&'a self, value: &'a PauliWord<A, H>) -> Self::Output {
        debug_assert_eq!(
            self.n_qubits(),
            value.n_qubits(),
            "#qubits mismatch, got {} and {}",
            self.n_qubits(),
            value.n_qubits()
        );
        self == value
    }
}
