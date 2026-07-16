// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Direct Heisenberg-picture Lindbladian evolution on an adaptive
//! Pauli-string basis.
//!
//! For a Hermitian Pauli Hamiltonian `H = Σ c_i P_i` and jump operators
//! `L_k = Σ_a λ_{k,a} P_{k,a}` (each a Hermitian-Pauli linear combination
//! with possibly complex coefficients) with rates `γ_k ≥ 0`, the adjoint
//! Lindbladian acts on a single Pauli string `p` as
//!
//! ```text
//! L*(p) = i [H, p] + Σ_k γ_k ( L_k† p L_k − 1/2 {L_k† L_k, p} ).
//! ```
//!
//! Two jump shapes are supported with separate code paths:
//!
//! - **Hermitian Pauli** (`L = P`, `λ ∈ ℝ`): the dissipator collapses to a
//!   diagonal `-2γ` on Pauli strings that anti-commute with `P`. Same fast
//!   path used by every dephasing-style model.
//!
//! - **General** (complex `λ_a`, e.g. `σ± = (X ± iY)/2`): the dissipator
//!   becomes a double sum `Σ_{a,b} λ_a* λ_b P_a p P_b` plus a Pauli-
//!   linear-combination anti-commutator with `L†L`, which is precomputed
//!   once at construction. Intermediate coefficients are complex; the
//!   result is real because `L*` preserves Hermiticity, so we cast back
//!   to `f64` at the boundary (with a debug-only check that `|Im|` is at
//!   FP noise).
//!
//! Pauli strings are stored as [`ppvm_pauli_word::word::PauliWord`] backed by
//! two 64-bit chunks (≤128 qubits; four 32-bit chunks on 32-bit targets)
//! with cached hashes for fast HashMap lookup. The hot-path commutator/
//! product loops bypass the higher-level word API and operate directly on
//! the raw chunks for speed.

mod algebra;
mod basis;
pub mod config;
pub mod error;
pub(crate) mod expm;
mod spec;
mod step;
mod word;

/// Matrix-free / quspin-expm-backed `exp(dt·L*)·b` engine. See module docs.
pub(crate) mod mf_expm;

/// Per-step orbit-rep evolution under translation symmetry, with a
/// phase-aware complex action. See module docs.
pub mod orbit_rep;

pub use basis::build_basis_index;
pub use config::PcStepConfig;
pub use error::Error;
pub use spec::{JumpInput, LindbladSpec};
pub use step::PcStepTimings;
pub use word::{MAX_QUBITS, Word, codes_from_word, parse_pauli_string, word_from_codes};

#[cfg(test)]
mod tests;

