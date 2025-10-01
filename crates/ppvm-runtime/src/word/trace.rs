use super::data::PauliWord;
use crate::{
    pattern::{Contains, PauliPattern},
    traits::{PauliStorage, Trace},
};

impl<'a, A: PauliStorage + 'a> Trace<'a, PauliWord<A>> for PauliWord<A> {
    type Output = bool;
    fn trace(&'a self, value: &'a PauliWord<A>) -> Self::Output {
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

impl<'a, A: PauliStorage + 'a> Trace<'a, PauliPattern> for PauliWord<A> {
    type Output = bool;
    fn trace(&'a self, value: &'a PauliPattern) -> Self::Output {
        value.contains(&self)
    }
}
