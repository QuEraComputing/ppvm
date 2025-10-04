use super::coefficient::Coefficient;
use super::map::ACMap;
use super::storage::PauliStorage;

pub trait Strategy: Default + Clone + Copy {
    /// Given the number of qubits, predict the initial capacity of the map.
    /// Ideally this is about guessing the maximum Pauli paths will be generated
    /// during the computation, the more precise the better.
    fn capacity(&self, n_qubits: usize) -> usize;
    fn truncate<S, V, H, M>(&self, map: &mut M)
    where
        S: PauliStorage,
        V: Coefficient,
        H: std::hash::BuildHasher + Clone + Default,
        M: ACMap<S, V, H>;
}

#[derive(Debug, Clone, Default, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NoStrategy;

impl Strategy for NoStrategy {
    fn capacity(&self, n_qubits: usize) -> usize {
        // in exact simulation, let's guess there will be 4^n / 2 = 2^(2n - 1) paths
        1 << (2 * n_qubits - 1)
    }

    fn truncate<S, V, H, M>(&self, _map: &mut M)
    where
        S: PauliStorage,
        V: Coefficient,
        H: std::hash::BuildHasher + Clone + Default,
        M: ACMap<S, V, H>,
    {
    }
}
