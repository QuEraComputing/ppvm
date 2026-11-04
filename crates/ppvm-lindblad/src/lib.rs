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
//! Pauli strings are stored as [`Word`]: an X and a Z bit plane of `C`
//! chunks each (64-bit, or 32-bit on 32-bit targets) and nothing else, so a
//! 128-qubit word is 32 bytes. Hashes are computed from the planes when a table
//! needs them. [`GeneralWord`] is the corresponding `ppvm_pauli_word` type taken
//! by `ppvm_pauli_sum` symmetry routines. The crate is const-generic in `C`
//! ([`Word`], [`LindbladSpec`]); the default is the 128-qubit width,
//! and [`chunks_for`] picks the narrowest of the 128/256/512-qubit widths
//! for a register. The hot-path commutator/
//! product loops bypass the higher-level word API and operate directly on
//! the raw chunks for speed.
//!
//! Real predictor-corrector steps use [`ShortlistScoreProposer`]: shortlist
//! external words by individual Hamiltonian contributions, then rank their
//! complete signed leakage. [`StreamingProposer`] remains available through
//! [`LindbladSpec::pc_step_with_proposer`] and handles general dissipators.

mod algebra;
mod basis;
mod candidates;
mod compact_index;
pub mod config;
mod direct_expm;
pub mod error;
pub(crate) mod expm;
mod indexed_candidates;
mod kossakowski;
mod pauli_plan;
#[cfg(test)]
mod reference;
mod scalar;
pub mod sector;
mod shortlist;
mod spec;
mod step;
mod truncate;
mod unit_expm;
mod word;

/// Matrix-free / quspin-expm-backed `exp(dt·L*)·b` engine. See module docs.
pub(crate) mod mf_expm;

pub use basis::{PcStepWorkspace, build_basis_index};
pub use candidates::{CandidateBuffer, CandidateProposer, StreamingProposer};
pub use config::PcStepConfig;
pub use error::Error;
pub use sector::{Sector, canonicalize_basis_to_rep};
pub use shortlist::ShortlistScoreProposer;
pub use spec::{JumpInput, LindbladSpec};
pub use step::PcStepTimings;
pub use word::{
    CHUNK_BITS, GeneralWord, MAX_QUBITS, MAX_SUPPORTED_QUBITS, W_CHUNKS, WIDTHS, Word, chunks_for,
    codes_from_word, max_qubits, parse_pauli_string, word_from_codes,
};

#[cfg(test)]
mod tests;
