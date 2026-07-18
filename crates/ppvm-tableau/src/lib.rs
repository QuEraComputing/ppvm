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
//!
//! ## Data-parallel / GPU-offloadable surface
//!
//! A handful of primitives in this crate operate on whole contiguous arrays
//! at once, with no per-element host closures — the shape a thread pool or a
//! CUDA kernel would want:
//!
//! - Every tableau row ([`PhasedPauliWord`](ppvm_pauli_word::phase::PhasedPauliWord))
//!   stores its Pauli word as fixed-size `bytemuck::Pod` integer bit-planes
//!   (`xbits`/`zbits`, backed by `[u8; N]` or `[u64; N]`), reachable as
//!   contiguous raw integer slices via `as_raw_slice`/`as_raw_mut_slice` —
//!   plain-old-data, not a bit-addressed abstraction.
//! - The Clifford single-/two-qubit gates
//!   (`crates/ppvm-tableau/src/gates/clifford.rs`) all run through one shared
//!   implementation per gate that loops over rows operating directly on
//!   those raw integer slices with a hoisted word-index/bit/mask, bypassing
//!   `bitvec`'s per-bit bounds checks inside the loop. Every caller —
//!   [`data::Tableau`], [`data::GeneralizedTableau`] (via the
//!   `impl_generalized_tableau_clifford*` macros), and the fused batch path
//!   below — funnels through that one implementation.
//! - [`CliffordBatch`](ppvm_traits::traits::CliffordBatch) methods (`x_many`,
//!   `cz_many`, …) and [`data::GeneralizedTableau::cz_block`] /
//!   `cz_block_pairs` / `cz_block_pairs_cross_word` apply a gate to a whole
//!   contiguous block of qubits as bulk masked slice operations — the entry
//!   points a thread pool or CUDA kernel would implement over the row
//!   planes.
//! - `GeneralizedTableau::branch_with_coefficients` (crate-private) transforms
//!   the whole coefficient array in one pass. Its `#[cfg(feature = "rayon")]`
//!   path (`branch_coefficients_parallel`, gated by `RAYON_COEFF_THRESHOLD`)
//!   computes the per-entry branch/non-branch coefficient math in parallel,
//!   with no shared mutable state — the accumulation of those results into
//!   the coefficient map afterwards remains sequential.

/// Core [`Tableau`](data::Tableau) and
/// [`GeneralizedTableau`](data::GeneralizedTableau) types.
pub mod data;
/// `Display` implementations for tableau types.
pub mod display;
/// Pauli-string expectation values and pattern traces.
pub mod expectation;
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
