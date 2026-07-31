// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `PauliSum` / `LossyPauliSum` Pauli-propagation engine, truncation
//! strategies, and pre-built `Config` bundles. Built on `ppvm-traits` and
//! `ppvm-pauli-word`.
pub mod config;
pub mod strategy;
pub mod sum;
pub mod symmetry;

/// Drop-in replacement for the old `ppvm_runtime::prelude`.
pub mod prelude {
    pub use crate::config;
    pub use crate::config::Config;
    pub use crate::sum::{PauliSum, impl_op_mul_assign_coefficient};
    pub use ppvm_pauli_word::loss::LossyPauliWord;
    pub use ppvm_pauli_word::pattern::PauliPattern;
    pub use ppvm_pauli_word::phase::PhasedPauliWord;
    pub use ppvm_pauli_word::word::PauliWord;
    pub use ppvm_traits::char::Pauli;
    pub use ppvm_traits::traits::*;
}
