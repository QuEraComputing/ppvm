use crate::traits::*;

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct CombinedStrategy<S1: Strategy, S2: Strategy>(pub S1, pub S2);

impl<S1: Strategy, S2: Strategy> Strategy for CombinedStrategy<S1, S2> {
    fn capacity(&self, n_qubits: usize) -> usize {
        self.0.capacity(n_qubits).min(self.1.capacity(n_qubits))
    }

    fn truncate<S, V, H, M>(&self, map: &mut M)
    where
        S: crate::prelude::PauliStorage,
        V: Coefficient,
        H: std::hash::BuildHasher + Clone + Default,
        M: crate::prelude::ACMap<S, V, H>,
    {
        self.0.truncate(map);
        self.1.truncate(map);
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaxPauliWeight(pub usize);

impl MaxPauliWeight {
    pub fn max_weight(&self) -> usize {
        self.0
    }
}

impl Default for MaxPauliWeight {
    fn default() -> Self {
        Self(usize::MAX)
    }
}

impl Strategy for MaxPauliWeight {
    fn capacity(&self, n_qubits: usize) -> usize {
        // the number here should scale binomially, but that can get large
        // since the capacity has a direct impact on performance, let's be conservative
        n_qubits * 10
    }

    fn truncate<S, V, H, M>(&self, map: &mut M)
    where
        S: PauliStorage,
        V: Coefficient,
        H: std::hash::BuildHasher + Clone + Default,
        M: ACMap<S, V, H>,
    {
        map.retain(|k, _| k.weight() <= self.max_weight());
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoefficientThreshold(pub f64);

impl Default for CoefficientThreshold {
    fn default() -> Self {
        Self(1e-12)
    }
}

impl Strategy for CoefficientThreshold {
    fn capacity(&self, n_qubits: usize) -> usize {
        // clearing maps scales as O(capacity), so let's be conservative here
        n_qubits * 10
    }

    fn truncate<S, V, H, M>(&self, map: &mut M)
    where
        S: PauliStorage,
        V: Coefficient,
        H: std::hash::BuildHasher + Clone + Default,
        M: ACMap<S, V, H>,
    {
        map.retain(|_, v| !v.cutoff(self.0));
    }
}
