// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Core runtime for the ppvm quantum-circuit simulator.
//!
//! `ppvm-runtime` defines the low-level types every ppvm backend uses:
//! single-qubit [`char::Pauli`] symbols, packed [`word::PauliWord`] strings,
//! phased and lossy variants, the [`sum::PauliSum`] dictionary of Pauli
//! strings to coefficients, and the [`traits`] hierarchy of gates,
//! measurements, and noise channels.
//!
//! The crate is generic over a [`config::Config`] bundle that fixes
//! storage, coefficient type, hashing, and truncation strategy at compile
//! time. Higher-level crates (`ppvm-tableau`, `ppvm-sym`) build on top of
//! these types.
//!
//! # Quick example
//!
//! Compute `<ZZ>` after the GHZ-preparation circuit `H(0); CNOT(0, 1)` by
//! propagating `ZZ` backwards in the Heisenberg picture:
//!
//! ```
//! use ppvm_runtime::prelude::*;
//!
//! let mut state: PauliSum<config::indexmap::ByteFxHashF64<1>> =
//!     PauliSum::builder().n_qubits(2).build();
//! state += ("ZZ", 1.0);
//!
//! // Circuit is H(0); CNOT(0, 1) — apply in reverse for Heisenberg propagation.
//! state.cnot(0, 1);
//! state.h(0);
//!
//! // ZZ has propagated to a single term.
//! assert_eq!(state.len(), 1);
//! ```

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
/// Truncation strategies that decide when low-weight or low-magnitude
/// terms are dropped from a [`PauliSum`](sum::PauliSum).
pub mod strategy;
/// Lattice translation-symmetry groups and the orbit-representative
/// merging primitive used by operator-space evolution methods.
pub mod symmetry;
/// The [`PauliSum`](sum::PauliSum) type and its builder.
pub mod sum;
/// Gate, measurement, and noise-channel traits shared across backends.
pub mod traits;
/// Packed [`PauliWord`](word::PauliWord) — a fixed-width Pauli string
/// stored as `u8`-packed `Pauli` digits.
pub mod word;

/// Convenience re-exports of the most common types and traits.
///
/// ```
/// use ppvm_runtime::prelude::*;
/// let state: PauliSum<config::indexmap::ByteFxHashF64<1>> =
///     PauliSum::builder().n_qubits(2).build();
/// assert_eq!(state.n_qubits(), 2);
/// ```
pub mod prelude {
    pub use crate::char::Pauli;
    pub use crate::config;
    pub use crate::config::Config;
    pub use crate::loss::LossyPauliWord;
    pub use crate::pattern::PauliPattern;
    pub use crate::phase::PhasedPauliWord;
    pub use crate::sum::{PauliSum, impl_op_mul_assign_coefficient};
    pub use crate::traits::*;
    pub use crate::word::PauliWord;
}
