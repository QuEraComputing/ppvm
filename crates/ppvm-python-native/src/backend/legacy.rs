// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::complex::Complex64;
use ppvm_pauli_sum_legacy::prelude::*;
use ppvm_pauli_sum_legacy::strategy::{
    CoefficientThreshold, CombinedStrategy, MaxLossWeight, MaxPauliWeight,
};
use ppvm_tableau_sum_legacy::storage::vec::VecStorage;

pub use ppvm_pauli_sum_legacy as pauli_sum;
pub use ppvm_tableau_legacy as tableau;
pub use ppvm_tableau_sum_legacy as tableau_sum;

pub type OrdinaryPolicy = CombinedStrategy<CoefficientThreshold, MaxPauliWeight>;
pub type LossPolicy =
    CombinedStrategy<CombinedStrategy<CoefficientThreshold, MaxPauliWeight>, MaxLossWeight>;
pub type OrdinaryConfig<const N: usize> = config::indexmap::ByteFxHashF64<N, OrdinaryPolicy>;
pub type LossConfig<const N: usize> =
    config::indexmap::ByteFxHashF64<N, LossPolicy, LossyPauliWord<[u8; N]>>;
pub type OrdinaryPauliSum<const N: usize> = PauliSum<OrdinaryConfig<N>>;
pub type LossyPauliSum<const N: usize> = PauliSum<LossConfig<N>>;

pub type TableauConfig<const N: usize> = config::fx64hash::Byte8F64<N>;
pub type GeneralizedTableau<const N: usize, I> =
    ppvm_tableau_legacy::data::GeneralizedTableau<TableauConfig<N>, I>;
pub type GeneralizedTableauSum<const N: usize, I> =
    ppvm_tableau_sum_legacy::data::GeneralizedTableauSum<
        TableauConfig<N>,
        I,
        Vec<(Complex64, I)>,
        VecStorage<TableauConfig<N>, I, Vec<(Complex64, I)>>,
    >;
pub type MixtureSampler<const N: usize, I> =
    ppvm_tableau_sum_legacy::sampler::Sampler<TableauConfig<N>, I>;
