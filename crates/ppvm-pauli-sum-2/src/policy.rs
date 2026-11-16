// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_traits_2::{Coefficient, Retain, Word, loss::LossInfo};

/// A truncation rule, applied only when the caller requests it.
pub trait Policy<K, C> {
    /// Initial capacity hint. Zero lets storage grow as terms are added.
    fn capacity(&self, _n_sites: usize) -> usize {
        0
    }

    /// Remove terms according to this policy.
    fn truncate<M: Retain<K, C>>(&self, storage: &mut M);
}

/// Keep all terms, including those with zero coefficients.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct NoPolicy;

impl<K, C> Policy<K, C> for NoPolicy {
    fn truncate<M: Retain<K, C>>(&self, _storage: &mut M) {}
}

/// Keep terms whose coefficient magnitude is at least `threshold`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoefficientThreshold {
    /// Inclusive lower bound on coefficient magnitude.
    pub threshold: f64,
}

impl Default for CoefficientThreshold {
    fn default() -> Self {
        Self { threshold: 1e-12 }
    }
}

impl<K, C: Coefficient> Policy<K, C> for CoefficientThreshold {
    fn capacity(&self, n_sites: usize) -> usize {
        n_sites.saturating_mul(10)
    }

    fn truncate<M: Retain<K, C>>(&self, storage: &mut M) {
        storage.retain(|_, coeff| coeff.magnitude() >= self.threshold);
    }
}

/// Keep words with at most this many non-identity sites, including lost sites.
/// The default, `usize::MAX`, disables this limit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaxPauliWeight(pub usize);

impl Default for MaxPauliWeight {
    fn default() -> Self {
        Self(usize::MAX)
    }
}

impl<K: Word, C> Policy<K, C> for MaxPauliWeight {
    fn capacity(&self, n_sites: usize) -> usize {
        n_sites.saturating_mul(10)
    }

    fn truncate<M: Retain<K, C>>(&self, storage: &mut M) {
        if self.0 != usize::MAX {
            storage.retain(|word, _| word.weight() <= self.0);
        }
    }
}

/// Keep words with at most this many lost sites; defaults to `10`.
/// Ordinary Pauli words have loss weight zero. `usize::MAX` disables the limit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaxLossWeight(pub usize);

impl Default for MaxLossWeight {
    fn default() -> Self {
        Self(10)
    }
}

impl<K: LossInfo, C> Policy<K, C> for MaxLossWeight {
    fn capacity(&self, n_sites: usize) -> usize {
        n_sites.saturating_mul(10)
    }

    fn truncate<M: Retain<K, C>>(&self, storage: &mut M) {
        if self.0 != usize::MAX {
            storage.retain(|word, _| word.loss_weight() <= self.0);
        }
    }
}

/// Apply two policies in order, retaining terms that survive both.
/// Uses the smaller initial capacity hint; policies can be nested.
///
/// ```
/// use ppvm_pauli_sum_2::{CoefficientThreshold, CombinedPolicy, MaxPauliWeight, PauliSum};
///
/// let policy = CombinedPolicy(CoefficientThreshold { threshold: 0.5 }, MaxPauliWeight(1));
/// let mut sum = PauliSum::from_terms_with_policy(
///     2,
///     policy,
///     [("II".into(), 0.0), ("XI".into(), 0.5), ("XY".into(), 2.0)],
/// ).preserving(["II".into()]);
///
/// sum.truncate();
/// assert_eq!(sum.get(&"II".into()), Some(0.0)); // Preserved despite the threshold.
/// assert_eq!(sum.get(&"XI".into()), Some(0.5)); // Meets both inclusive limits.
/// assert!(!sum.contains_key(&"XY".into())); // Exceeds the weight limit.
/// ```
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct CombinedPolicy<P1, P2>(pub P1, pub P2);

impl<K, C, P1, P2> Policy<K, C> for CombinedPolicy<P1, P2>
where
    P1: Policy<K, C>,
    P2: Policy<K, C>,
{
    fn capacity(&self, n_sites: usize) -> usize {
        self.0.capacity(n_sites).min(self.1.capacity(n_sites))
    }

    fn truncate<M: Retain<K, C>>(&self, storage: &mut M) {
        self.0.truncate(storage);
        self.1.truncate(storage);
    }
}
