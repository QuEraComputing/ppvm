// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::traits::*;

/// Two strategies run in sequence: `S1` then `S2`.
///
/// The capacity hint is the smaller of the two; truncation applies both
/// policies in order, so use this when you want, e.g.,
/// "drop by coefficient threshold *and* cap maximum Pauli weight".
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct CombinedStrategy<S1: Strategy, S2: Strategy>(pub S1, pub S2);

impl<S1: Strategy, S2: Strategy> Strategy for CombinedStrategy<S1, S2> {
    fn capacity(&self, n_qubits: usize) -> usize {
        self.0.capacity(n_qubits).min(self.1.capacity(n_qubits))
    }

    fn truncate<S, V, H, M, W>(&self, map: &mut M)
    where
        S: crate::prelude::PauliStorage,
        V: Coefficient,
        H: std::hash::BuildHasher + Clone + Default,
        W: PauliWordTrait,
        M: crate::prelude::ACMap<S, V, H, W>,
    {
        self.0.truncate(map);
        self.1.truncate(map);
    }
}

/// Drop terms whose Pauli weight (number of non-identity slots) exceeds
/// the given bound.
///
/// # Examples
///
/// ```
/// use ppvm_runtime::prelude::*;
/// use ppvm_runtime::strategy::MaxPauliWeight;
///
/// // Keep only weight-1 terms.
/// let strat = MaxPauliWeight(1);
/// let mut state: PauliSum<config::indexmap::ByteFxHashF64<1, MaxPauliWeight>> =
///     PauliSum::builder().n_qubits(3).strategy(strat).build();
/// state += ("XII", 1.0);              // weight 1, kept
/// state += ("XYI", 1.0);              // weight 2, dropped
/// state.truncate();
/// assert_eq!(state.len(), 1);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaxPauliWeight(pub usize);

impl MaxPauliWeight {
    /// Maximum Pauli weight retained.
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

    fn truncate<S, V, H, M, W>(&self, map: &mut M)
    where
        S: PauliStorage,
        V: Coefficient,
        H: std::hash::BuildHasher + Clone + Default,
        W: PauliWordTrait,
        M: ACMap<S, V, H, W>,
    {
        // `usize::MAX` is the conventional "disabled" sentinel (matches
        // `Default::default()`). Skip the retain pass entirely — saves a
        // full bucket walk per `truncate()` call when callers (e.g. the
        // Python binding) pass `MaxPauliWeight(usize::MAX)` to opt out
        // of weight truncation without changing the strategy type.
        if self.0 == usize::MAX {
            return;
        }
        map.retain(|k, _| k.weight() <= self.max_weight());
    }
}

/// Drop terms whose coefficient magnitude falls below the given threshold.
/// Defaults to `1e-12`.
///
/// # Examples
///
/// ```
/// use ppvm_runtime::prelude::*;
/// use ppvm_runtime::strategy::CoefficientThreshold;
///
/// let strict = CoefficientThreshold(1e-6);
/// let mut state: PauliSum<config::indexmap::ByteFxHashF64<1, CoefficientThreshold>> =
///     PauliSum::builder().n_qubits(2).strategy(strict).build();
/// state += ("ZZ", 1.0);
/// state += ("XX", 1e-9);              // below threshold
/// state.truncate();
/// assert_eq!(state.len(), 1);         // the XX term was dropped
/// ```
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

    fn truncate<S, V, H, M, W>(&self, map: &mut M)
    where
        S: PauliStorage,
        V: Coefficient,
        H: std::hash::BuildHasher + Clone + Default,
        W: PauliWordTrait,
        M: ACMap<S, V, H, W>,
    {
        map.retain(|_, v| !v.cutoff(self.0));
    }
}

/// Drop terms whose loss weight (number of lost qubits) exceeds the
/// given bound. Only meaningful for [`LossyPauliWord`](crate::loss::LossyPauliWord).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaxLossWeight(pub usize);

impl Default for MaxLossWeight {
    fn default() -> Self {
        Self(10)
    }
}

impl Strategy for MaxLossWeight {
    fn capacity(&self, n_qubits: usize) -> usize {
        // the number here should scale binomially, but that can get large
        // since the capacity has a direct impact on performance, let's be conservative
        n_qubits * 10
    }

    fn truncate<S, V, H, M, W>(&self, map: &mut M)
    where
        S: PauliStorage,
        V: Coefficient,
        H: std::hash::BuildHasher + Clone + Default,
        M: ACMap<S, V, H, W>,
        W: PauliWordTrait,
    {
        // Skip the retain pass when callers pass `usize::MAX` to opt out
        // of loss-weight truncation without changing the strategy type.
        if self.0 == usize::MAX {
            return;
        }
        map.retain(|k, _| k.loss_weight() <= self.0);
    }
}
