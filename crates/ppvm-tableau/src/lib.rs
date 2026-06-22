// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Generalized stabilizer-tableau simulator built on top of the ppvm core
//! crates (`ppvm-traits` and `ppvm-pauli-word`).
//!
//! Provides a forward-evolving state representation using the stabilizer
//! formalism, extended to non-Clifford gates by tracking a sparse vector
//! of coefficients indexed over bitstrings. The two top-level types are
//! [`data::Tableau`] (pure Clifford) and [`data::GeneralizedTableau`]
//! (Clifford + non-Clifford with sparse coefficient tracking).
//!
//! # Quick example
//!
//! Prepare a Bell pair and verify the two measurements are perfectly
//! correlated:
//!
//! ```
//! use ppvm_tableau::prelude::*;
//! use ppvm_pauli_sum::config::fxhash::ByteF64;
//!
//! let mut tab: GeneralizedTableau<ByteF64<1>> =
//!     GeneralizedTableau::new_with_seed(2, 1e-12, 0);
//! tab.h(0);
//! tab.cnot(0, 1);
//!
//! let r0 = LossyMeasure::measure(&mut tab, 0);
//! let r1 = LossyMeasure::measure(&mut tab, 1);
//! assert_eq!(r0, r1);
//! ```

/// Core [`Tableau`](data::Tableau) and
/// [`GeneralizedTableau`](data::GeneralizedTableau) types.
pub mod data;
/// `Display` implementations for tableau types.
pub mod display;
/// Gate implementations (Clifford, T, rotations).
pub mod gates;
/// Z-basis measurement, including loss-aware variants.
pub mod measure;

pub mod measure_all;

/// Noise channels: depolarizing, Pauli error, loss.
pub mod noise;
/// [`SparseVector`](sparsevec::SparseVector) trait and implementations.
pub mod sparsevec;
/// [`TableauIndex`](tableau_index::TableauIndex) — abstraction over
/// the bitstring index type (`usize`, `u128`, `U256`, …).
pub mod tableau_index;
/// [`TableauLike`](tableau_like::TableauLike) — shared trait for
/// stabilizer-tableau backends.
pub mod tableau_like;

/// Convenience re-exports for downstream code.
pub mod prelude {
    pub use crate::data::{GeneralizedTableau, Tableau};
    pub use crate::sparsevec::SparseVector;
    pub use crate::tableau_index::TableauIndex;
    pub use crate::tableau_like::TableauLike;
    pub use ppvm_pauli_sum::prelude::*;
}
