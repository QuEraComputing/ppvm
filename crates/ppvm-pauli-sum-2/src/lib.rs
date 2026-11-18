// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod policy;
mod store;
mod sum;

pub use policy::{
    CoefficientThreshold, CombinedPolicy, MaxLossWeight, MaxPauliWeight, NoPolicy, Policy,
};
pub use store::{HashMapStore, StoreAlloc};
pub use sum::{PreserveSet, Sum};

pub use ppvm_lossy_pauli_word_2::LossyPauliWord;
pub use ppvm_pauli_word_2::PauliWord;

/// A Pauli sum backed by a hash map, with real coefficients by default.
pub type PauliSum<C = f64, P = NoPolicy> = Sum<HashMapStore<PauliWord, C>, P>;

/// A Pauli sum whose words can represent lost qubits.
pub type LossyPauliSum<C = f64, P = NoPolicy> = Sum<HashMapStore<LossyPauliWord, C>, P>;
