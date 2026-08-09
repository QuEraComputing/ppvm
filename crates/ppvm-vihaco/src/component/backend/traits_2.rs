// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bnum::types::{U256, U512, U1024, U2048};
use ppvm_pauli_sum_2::{
    CoefficientThreshold, CombinedPolicy, HashMapStore, LossyPauliWord, MaxPauliWeight, PauliWord,
    Sum,
};
use ppvm_tableau_2::GeneralizedTableau;

use crate::device_info::PPVMDeviceInfo;

pub use ppvm_pauli_sum_2::PauliPattern;
pub use ppvm_tableau_2::prelude::*;

pub type Policy = CombinedPolicy<CoefficientThreshold, MaxPauliWeight>;

pub type Sum64 = Sum<
    HashMapStore<PauliWord<[u8; 8]>, f64>,
    CombinedPolicy<CoefficientThreshold, MaxPauliWeight>,
>;
pub type Sum128 = Sum<
    HashMapStore<PauliWord<[u8; 16]>, f64>,
    CombinedPolicy<CoefficientThreshold, MaxPauliWeight>,
>;
pub type Sum256 = Sum<
    HashMapStore<PauliWord<[u8; 32]>, f64>,
    CombinedPolicy<CoefficientThreshold, MaxPauliWeight>,
>;
pub type Sum512 = Sum<
    HashMapStore<PauliWord<[u8; 64]>, f64>,
    CombinedPolicy<CoefficientThreshold, MaxPauliWeight>,
>;
pub type Sum1024 = Sum<
    HashMapStore<PauliWord<[u8; 128]>, f64>,
    CombinedPolicy<CoefficientThreshold, MaxPauliWeight>,
>;
pub type Sum2048 = Sum<
    HashMapStore<PauliWord<[u8; 256]>, f64>,
    CombinedPolicy<CoefficientThreshold, MaxPauliWeight>,
>;

pub type Lossy64 = Sum<
    HashMapStore<LossyPauliWord<[u8; 8]>, f64>,
    CombinedPolicy<CoefficientThreshold, MaxPauliWeight>,
>;
pub type Lossy128 = Sum<
    HashMapStore<LossyPauliWord<[u8; 16]>, f64>,
    CombinedPolicy<CoefficientThreshold, MaxPauliWeight>,
>;
pub type Lossy256 = Sum<
    HashMapStore<LossyPauliWord<[u8; 32]>, f64>,
    CombinedPolicy<CoefficientThreshold, MaxPauliWeight>,
>;
pub type Lossy512 = Sum<
    HashMapStore<LossyPauliWord<[u8; 64]>, f64>,
    CombinedPolicy<CoefficientThreshold, MaxPauliWeight>,
>;
pub type Lossy1024 = Sum<
    HashMapStore<LossyPauliWord<[u8; 128]>, f64>,
    CombinedPolicy<CoefficientThreshold, MaxPauliWeight>,
>;
pub type Lossy2048 = Sum<
    HashMapStore<LossyPauliWord<[u8; 256]>, f64>,
    CombinedPolicy<CoefficientThreshold, MaxPauliWeight>,
>;

#[cfg(not(target_arch = "wasm32"))]
pub type Tab64 = GeneralizedTableau<[usize; 1], usize>;
#[cfg(target_arch = "wasm32")]
pub type Tab64 = GeneralizedTableau<[usize; 2], usize>;
#[cfg(not(target_arch = "wasm32"))]
pub type Tab128 = GeneralizedTableau<[usize; 2], u128>;
#[cfg(target_arch = "wasm32")]
pub type Tab128 = GeneralizedTableau<[usize; 4], u128>;
#[cfg(not(target_arch = "wasm32"))]
pub type Tab256 = GeneralizedTableau<[usize; 4], U256>;
#[cfg(target_arch = "wasm32")]
pub type Tab256 = GeneralizedTableau<[usize; 8], U256>;
#[cfg(not(target_arch = "wasm32"))]
pub type Tab512 = GeneralizedTableau<[usize; 8], U512>;
#[cfg(target_arch = "wasm32")]
pub type Tab512 = GeneralizedTableau<[usize; 16], U512>;
#[cfg(not(target_arch = "wasm32"))]
pub type Tab1024 = GeneralizedTableau<[usize; 16], U1024>;
#[cfg(target_arch = "wasm32")]
pub type Tab1024 = GeneralizedTableau<[usize; 32], U1024>;
#[cfg(not(target_arch = "wasm32"))]
pub type Tab2048 = GeneralizedTableau<[usize; 32], U2048>;
#[cfg(target_arch = "wasm32")]
pub type Tab2048 = GeneralizedTableau<[usize; 64], U2048>;

pub struct Backend;

impl Backend {
    pub fn policy(info: &PPVMDeviceInfo) -> Policy {
        CombinedPolicy(
            CoefficientThreshold {
                threshold: info.coefficient_threshold,
            },
            MaxPauliWeight(info.max_pauli_weight.unwrap_or(usize::MAX)),
        )
    }

    pub fn render(value: &impl std::fmt::Display) -> String {
        value.to_string()
    }
}

/// The generator an executor owns (see `§ Where the randomness lives` in
/// `component/backend/mod.rs`). `None` seeds from OS entropy.
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

macro_rules! tableau_constructor {
    ($new:ident, $ty:ty) => {
        /// The `-2` tableau is pure state; the seed goes to the executor's RNG.
        pub fn $new(n_qubits: usize, threshold: f64) -> $ty {
            <$ty>::new(n_qubits, threshold)
        }
    };
}

impl Backend {
    tableau_constructor!(tab64, Tab64);
    tableau_constructor!(tab128, Tab128);
    tableau_constructor!(tab256, Tab256);
    tableau_constructor!(tab512, Tab512);
    tableau_constructor!(tab1024, Tab1024);
    tableau_constructor!(tab2048, Tab2048);
}

macro_rules! sum_constructor {
    ($name:ident, $ty:ty, $word:ty) => {
        pub fn $name(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> $ty {
            let mut state: $ty = Sum::with_policy(info.n_qubits, Self::policy(info));
            for (word, coeff) in terms {
                state += (<$word>::from(word.as_str()), *coeff);
            }
            state
        }
    };
}

impl Backend {
    sum_constructor!(sum64, Sum64, PauliWord<[u8; 8]>);
    sum_constructor!(sum128, Sum128, PauliWord<[u8; 16]>);
    sum_constructor!(sum256, Sum256, PauliWord<[u8; 32]>);
    sum_constructor!(sum512, Sum512, PauliWord<[u8; 64]>);
    sum_constructor!(sum1024, Sum1024, PauliWord<[u8; 128]>);
    sum_constructor!(sum2048, Sum2048, PauliWord<[u8; 256]>);
    sum_constructor!(lossy64, Lossy64, LossyPauliWord<[u8; 8]>);
    sum_constructor!(lossy128, Lossy128, LossyPauliWord<[u8; 16]>);
    sum_constructor!(lossy256, Lossy256, LossyPauliWord<[u8; 32]>);
    sum_constructor!(lossy512, Lossy512, LossyPauliWord<[u8; 64]>);
    sum_constructor!(lossy1024, Lossy1024, LossyPauliWord<[u8; 128]>);
    sum_constructor!(lossy2048, Lossy2048, LossyPauliWord<[u8; 256]>);
}
