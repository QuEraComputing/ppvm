// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Core runtime for the ppvm quantum-circuit simulator.
//!
//! `ppvm-runtime` defines the low-level types every ppvm backend uses:
//! single-qubit [`char::Pauli`] symbols, packed [`word::PauliWord`] strings,
//! phased and lossy variants, the [`traits`] hierarchy of gates,
//! measurements, and noise channels, and the [`config::Config`] bundle
//! trait that fixes storage, coefficient type, hashing, and truncation
//! strategy at compile time. Higher-level crates (`ppvm-paulisum`,
//! `ppvm-tableau`, `ppvm-sym`) build on top of these types.

/// Single-qubit Pauli alphabet (`I`, `X`, `Y`, `Z`, and a loss marker `L`).
pub mod char;
/// Pre-bundled [`Config`](config::Config) implementations and the trait
/// itself, which fixes storage / coefficient / hashing / strategy at
/// compile time.
pub mod config;
/// Loss-aware Pauli words that carry a parallel bitmap of lost qubits.
pub mod loss;
mod map;
/// Pauli "patterns" — sets of compatible Pauli words used for masking
/// and selective gate application.
pub mod pattern;
/// Phased Pauli words: a Pauli word paired with one of the four
/// fourth-roots of unity `±1, ±i`.
pub mod phase;
/// Gate, measurement, and noise-channel traits shared across backends.
pub mod traits;
/// Packed [`PauliWord`](word::PauliWord) — a fixed-width Pauli string
/// stored as `u8`-packed `Pauli` digits.
pub mod word;

/// Convenience re-exports of the most common types and traits.
pub mod prelude {
    pub use crate::char::Pauli;
    pub use crate::config;
    pub use crate::config::Config;
    pub use crate::loss::LossyPauliWord;
    pub use crate::pattern::PauliPattern;
    pub use crate::phase::PhasedPauliWord;
    pub use crate::traits::*;
    pub use crate::word::PauliWord;
}
