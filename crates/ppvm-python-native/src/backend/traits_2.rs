// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_pauli_sum_2::{
    CoefficientThreshold, CombinedPolicy, IndexMapStore, LossyPauliWord, MaxLossWeight,
    MaxPauliWeight, PauliWord, Sum,
};

pub use ppvm_pauli_sum_2 as pauli_sum;
pub use ppvm_tableau_2 as tableau;

pub type OrdinaryPolicy = CombinedPolicy<CoefficientThreshold, MaxPauliWeight>;
pub type LossPolicy =
    CombinedPolicy<CombinedPolicy<CoefficientThreshold, MaxPauliWeight>, MaxLossWeight>;
pub type OrdinaryPauliSum<const N: usize> =
    Sum<IndexMapStore<PauliWord<[u8; N]>, f64>, OrdinaryPolicy>;
pub type LossyPauliSum<const N: usize> =
    Sum<IndexMapStore<LossyPauliWord<[u8; N]>, f64>, LossPolicy>;

pub type GeneralizedTableau<const N: usize, I> = ppvm_tableau_2::GeneralizedTableau<[usize; N], I>;
pub type GeneralizedTableauSum<const N: usize, I> =
    ppvm_tableau_2::GeneralizedTableauSum<[usize; N], I>;
pub type MixtureSampler<const N: usize, I> =
    ppvm_tableau_2::MixtureSampler<[usize; N], I, fxhash::FxBuildHasher>;

/// The generator a `#[pyclass]` wrapper owns (see `§ Where the randomness
/// lives` in `backend/mod.rs`). `None` seeds from OS entropy.
///
/// `SmallRng` matches what the legacy tableau embedded, so a given seed
/// reproduces the same stream across both backends.
pub fn make_rng(seed: Option<u64>) -> rand::rngs::SmallRng {
    use rand::SeedableRng;
    match seed {
        Some(s) => rand::rngs::SmallRng::seed_from_u64(s),
        None => rand::make_rng(),
    }
}
