// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `PauliSum` and its supporting machinery: feature-gated
//! [`Config`](config::Config) bundles, truncation strategies, and the
//! ACMap impls for non-std map backends.
//!
//! Built on top of `ppvm-runtime`, which carries the primitives
//! (`PauliWord`, `Config` trait, gate/noise traits, etc.) shared with
//! every other ppvm backend.

/// [`Config`](config::Config) bundles: re-exports from `ppvm-runtime`
/// plus the feature-gated ones defined here.
pub mod config;
/// Truncation strategies applied when
/// [`PauliSum::truncate`](sum::PauliSum::truncate) runs.
pub mod strategy;
/// The [`PauliSum`](sum::PauliSum) type and its gate/noise/projection impls.
pub mod sum;

/// Convenience re-exports of the most common types and traits.
pub mod prelude {
    pub use crate::config;
    pub use crate::config::Config;
    pub use crate::strategy::{
        CoefficientThreshold, CombinedStrategy, MaxLossWeight, MaxPauliWeight,
    };
    pub use crate::sum::{PauliSum, impl_op_mul_assign_coefficient};
    pub use ppvm_runtime::char::Pauli;
    pub use ppvm_runtime::loss::LossyPauliWord;
    pub use ppvm_runtime::pattern::PauliPattern;
    pub use ppvm_runtime::phase::PhasedPauliWord;
    pub use ppvm_runtime::traits::*;
    pub use ppvm_runtime::word::PauliWord;
}
