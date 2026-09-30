// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Canonicalization of the stabilizer generators of a [`ppvm_tableau_2::Tableau`].
//!
//! These operations use the default tableau type and always track phases.
//! They update the paired destabilizers to preserve the Clifford frame.
//! Inputs must be valid pure-state frames; arbitrary or redundant generator
//! lists and generalized-tableau amplitudes are outside this crate's scope.
//!
//! Algorithms and documentation are adapted from
//! [QuantumClifford.jl](https://github.com/QuantumSavory/QuantumClifford.jl),
//! `src/canonicalization.jl`, `src/entanglement.jl`, and
//! `docs/src/canonicalization.md` (MIT, copyright 2023 Stefan Krastanov).
//! The original license is included in `LICENSE-QuantumClifford`.

mod canonical;
mod clipped;
mod gottesman;
mod rows;
mod rref;

pub use canonical::canonicalize;
pub use clipped::canonicalize_clip;
pub use gottesman::{GottesmanForm, canonicalize_gott};
pub use rref::canonicalize_rref;
