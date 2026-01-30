use crate::char::Pauli;
use crate::traits::PauliStorage;
use std::hash::Hash;

pub trait PauliIter {
    fn iter(&self) -> impl Iterator<Item = Pauli>;
}

pub trait PauliWordTrait<S: PauliStorage, V = fxhash::FxBuildHasher>:
    Clone + Hash + Eq + PauliIter
{
    fn new(nqubits: usize) -> Self;

    fn n_qubits(&self) -> usize;

    fn weight(&self) -> usize;

    fn rehash(&mut self);

    fn get(&self, index: usize) -> Pauli;

    fn get_multiple<const Q: usize>(&self, indices: [usize; Q]) -> Self {
        let mut result = Self::new(Q);
        for (i, &idx) in indices.iter().enumerate() {
            result.set(i, self.get(idx));
        }
        result
    }

    fn set_multiple<const Q: usize, B: PauliStorage>(&mut self, indices: [usize; Q], values: &Self);

    fn get_slice(&self, slice: std::ops::Range<usize>) -> Self;

    fn is(&self, index: usize, pauli: Pauli) -> bool;

    fn set(&mut self, index: usize, pauli: Pauli) -> &mut Self;

    #[inline(always)]
    fn set_new(&self, index: usize, pauli: Pauli) -> Self {
        if index >= self.n_qubits() {
            panic!("Index out of bounds");
        }
        let mut new = self.clone();
        new.set(index, pauli);
        new
    }

    #[inline(always)]
    fn set_new_2(&self, index_0: usize, pauli_0: Pauli, index_1: usize, pauli_1: Pauli) -> Self {
        let mut new = self.clone();
        new.set(index_0, pauli_0);
        new.set(index_1, pauli_1);
        new
    }
}
