// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! CPU baseline: ppvm's `PauliSum` configured like the Python `ppvm.PauliSum`.

use ppvm_pauli_sum::prelude::*;
use ppvm_pauli_sum::strategy::{CoefficientThreshold, CombinedStrategy, MaxPauliWeight};

use crate::Propagator;

/// The config behind `ppvm.PauliSum` (`PauliSumIndexMapFxHash<k>` in `_core`).
pub type CpuConfig<const B: usize> =
    config::indexmap::ByteFxHashF64<B, CombinedStrategy<CoefficientThreshold, MaxPauliWeight>>;

/// Storage bytes `ppvm.PauliSum` picks for `n_qubits`: the smallest power of
/// two strictly greater than `ceil(n_qubits / 8)`.
pub fn storage_bytes(n_qubits: usize) -> usize {
    let need = n_qubits.div_ceil(8);
    (0..16)
        .map(|k| 1usize << k)
        .find(|&b| b > need)
        .expect("n_qubits too large for ppvm.PauliSum storage")
}

/// `PauliSum` with `B` storage bytes, truncating after every gate.
#[derive(Clone)]
pub struct CpuPauliSum<const B: usize>(pub PauliSum<CpuConfig<B>>);

impl<const B: usize> CpuPauliSum<B> {
    /// Mirrors `PauliSum.new(n_qubits, terms, min_abs_coeff=cutoff)`.
    pub fn new(n_qubits: usize, terms: &[String], cutoff: f64) -> Self {
        let strategy = CombinedStrategy(CoefficientThreshold(cutoff), MaxPauliWeight(usize::MAX));
        let mut ps: PauliSum<CpuConfig<B>> = PauliSum::builder()
            .n_qubits(n_qubits)
            .strategy(strategy)
            .capacity(n_qubits)
            .build();
        for term in terms {
            ps += (term.as_str(), 1.0);
        }
        Self(ps)
    }

    /// All `(dense Pauli string, coefficient)` terms, sorted by string.
    pub fn terms(&self) -> Vec<(String, f64)> {
        let mut terms: Vec<_> = self.0.iter().map(|(k, c)| (k.to_string(), *c)).collect();
        terms.sort_by(|a, b| a.0.cmp(&b.0));
        terms
    }
}

impl<const B: usize> Propagator for CpuPauliSum<B> {
    fn rx(&mut self, q: usize, theta: f64) {
        self.0.rx(q, theta);
        self.0.truncate();
    }

    fn rxx(&mut self, a: usize, b: usize, theta: f64) {
        self.0.rxx(a, b, theta);
        self.0.truncate();
    }

    fn ryy(&mut self, a: usize, b: usize, theta: f64) {
        self.0.ryy(a, b, theta);
        self.0.truncate();
    }

    fn rzz(&mut self, a: usize, b: usize, theta: f64) {
        self.0.rzz(a, b, theta);
        self.0.truncate();
    }

    fn pauli_error(&mut self, q: usize, p: [f64; 3]) {
        self.0.pauli_error(q, p);
        self.0.truncate();
    }

    fn overlap_with_zero(&self) -> f64 {
        let zero_state: PauliPattern = "Z?*".into();
        self.0.trace(&zero_state)
    }

    fn len(&self) -> usize {
        self.0.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn storage_bytes_matches_python() {
        // paulisum.py: N = ceil(n/8); first 2**k with 2**k > N.
        assert_eq!(storage_bytes(4), 2);
        assert_eq!(storage_bytes(8), 2);
        assert_eq!(storage_bytes(9), 4);
        assert_eq!(storage_bytes(36), 8);
        assert_eq!(storage_bytes(128), 32);
        assert_eq!(storage_bytes(512), 128);
    }

    #[test]
    fn string_index_is_qubit_index() {
        let mut s = CpuPauliSum::<2>::new(3, &["ZII".to_owned()], 0.0);
        s.rx(0, 0.3);
        assert!((s.overlap_with_zero() - 0.3f64.cos()).abs() < 1e-12);
        assert_eq!(s.len(), 2);
    }

    #[test]
    fn terms_are_dense_strings() {
        // ppvm rx: Z -> cos θ Z + sin θ Y.
        let mut s = CpuPauliSum::<16>::new(70, &[format!("{}Z", "I".repeat(69))], 0.0);
        s.rx(69, 0.3);
        let terms = s.terms();
        assert_eq!(terms[0].0, format!("{}Y", "I".repeat(69)));
        assert!((terms[0].1 - 0.3f64.sin()).abs() < 1e-15);
        assert_eq!(terms[1].0, format!("{}Z", "I".repeat(69)));
    }
}
