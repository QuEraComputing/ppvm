// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bnum::types::{U256, U512, U1024, U2048};
use num::complex::Complex64;
use ppvm_pauli_sum_legacy::config::fx64hash::Byte8F64;
use ppvm_pauli_sum_legacy::config::indexmap::ByteFxHashF64;
use ppvm_pauli_sum_legacy::strategy::{CoefficientThreshold, CombinedStrategy, MaxPauliWeight};
use ppvm_tableau_legacy::data::GeneralizedTableau;

use crate::device_info::PPVMDeviceInfo;

pub use ppvm_tableau_legacy::prelude::*;

pub type Policy = CombinedStrategy<CoefficientThreshold, MaxPauliWeight>;
type SumConfig<const N: usize> = ByteFxHashF64<N, Policy>;
type LossyConfig<const N: usize> = ByteFxHashF64<N, Policy, LossyPauliWord<[u8; N]>>;

pub type Sum64 = PauliSum<SumConfig<8>>;
pub type Sum128 = PauliSum<SumConfig<16>>;
pub type Sum256 = PauliSum<SumConfig<32>>;
pub type Sum512 = PauliSum<SumConfig<64>>;
pub type Sum1024 = PauliSum<SumConfig<128>>;
pub type Sum2048 = PauliSum<SumConfig<256>>;

pub type Lossy64 = PauliSum<LossyConfig<8>>;
pub type Lossy128 = PauliSum<LossyConfig<16>>;
pub type Lossy256 = PauliSum<LossyConfig<32>>;
pub type Lossy512 = PauliSum<LossyConfig<64>>;
pub type Lossy1024 = PauliSum<LossyConfig<128>>;
pub type Lossy2048 = PauliSum<LossyConfig<256>>;

pub type Tab64 = GeneralizedTableau<Byte8F64<1>, usize, Vec<(Complex64, usize)>>;
pub type Tab128 = GeneralizedTableau<Byte8F64<2>, u128, Vec<(Complex64, u128)>>;
pub type Tab256 = GeneralizedTableau<Byte8F64<4>, U256, Vec<(Complex64, U256)>>;
pub type Tab512 = GeneralizedTableau<Byte8F64<8>, U512, Vec<(Complex64, U512)>>;
pub type Tab1024 = GeneralizedTableau<Byte8F64<16>, U1024, Vec<(Complex64, U1024)>>;
pub type Tab2048 = GeneralizedTableau<Byte8F64<32>, U2048, Vec<(Complex64, U2048)>>;

pub struct Backend;

impl Backend {
    pub fn policy(info: &PPVMDeviceInfo) -> Policy {
        CombinedStrategy(
            CoefficientThreshold(info.coefficient_threshold),
            MaxPauliWeight(info.max_pauli_weight.unwrap_or(usize::MAX)),
        )
    }

    pub fn render(value: &impl std::fmt::Display) -> String {
        value.to_string()
    }
}

macro_rules! tableau_constructors {
    ($new:ident, $seeded:ident, $ty:ty) => {
        pub fn $new(n_qubits: usize, threshold: f64) -> $ty {
            <$ty>::new(n_qubits, threshold)
        }

        pub fn $seeded(n_qubits: usize, threshold: f64, seed: u64) -> $ty {
            <$ty>::new_with_seed(n_qubits, threshold, seed)
        }
    };
}

impl Backend {
    tableau_constructors!(tab64, tab64_seeded, Tab64);
    tableau_constructors!(tab128, tab128_seeded, Tab128);
    tableau_constructors!(tab256, tab256_seeded, Tab256);
    tableau_constructors!(tab512, tab512_seeded, Tab512);
    tableau_constructors!(tab1024, tab1024_seeded, Tab1024);
    tableau_constructors!(tab2048, tab2048_seeded, Tab2048);
}

macro_rules! sum_constructor {
    ($name:ident, $ty:ty) => {
        pub fn $name(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> $ty {
            let mut state = PauliSum::builder()
                .n_qubits(info.n_qubits)
                .strategy(Self::policy(info))
                .build();
            for (word, coeff) in terms {
                state += (word.as_str(), *coeff);
            }
            state
        }
    };
}

impl Backend {
    sum_constructor!(sum64, Sum64);
    sum_constructor!(sum128, Sum128);
    sum_constructor!(sum256, Sum256);
    sum_constructor!(sum512, Sum512);
    sum_constructor!(sum1024, Sum1024);
    sum_constructor!(sum2048, Sum2048);
    sum_constructor!(lossy64, Lossy64);
    sum_constructor!(lossy128, Lossy128);
    sum_constructor!(lossy256, Lossy256);
    sum_constructor!(lossy512, Lossy512);
    sum_constructor!(lossy1024, Lossy1024);
    sum_constructor!(lossy2048, Lossy2048);
}
