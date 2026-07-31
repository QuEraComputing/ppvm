// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `ppvm-phased-pauli-word-2`: the phased Pauli word wrapper.
//!
//! `Phased<W>` wraps a base word `W` with an explicit `ℤ₄` phase — the standalone
//! phased Pauli operator (`i^φ · g(w)`). Unlike the bare `PauliWord` /
//! `LossyPauliWord`, whose Clifford drops every conjugation sign, `Phased<W>`
//! carries a **fused** [`Clifford`](ppvm_traits_2::Clifford) impl that *tracks*
//! the `ℤ₄` sign (`HYH = −Y`, `SXS† = Y`, the CNOT/CZ deltas). Rather than the
//! shared blanket `Clifford` (which would read each inner X/Z bit twice — once for
//! the sign, once for the bit op), it reads each bit once via `W: PauliBits`,
//! computes the sign, applies the bit update, and folds in the sign — so it
//! deliberately does **not** implement `BlanketClifford` (the marker gating the
//! blanket); see `clifford.rs`. It is also deliberately **not** `Indexable`: the
//! phase is part of its identity, so it is not a production map key
//! (`docs/design/traits-2-configuration-and-hashing.md`).
//!
//! `PhasedPauliWord = Phased<PauliWord>` is the domain alias. The phased product
//! reuses the base word's `KeyProduct` kernel (accumulating the emitted phase),
//! so the `phaseExp` boolean kernel is not duplicated; a lossy base has Clifford
//! conjugation but no product (loss breaks the twisted-product group).
//!
//! Trait contracts and their Lean-validated semantics (the `𝒫₁` group in
//! `Phase.lean`, the conjugation signs in `Conjugation.lean`, the ℤ[i] matrix
//! grounding in `Matrix.lean`) follow
//! [`docs/design/traits-2-configuration-and-hashing.md`] and `lean/PPVM/**`;
//! ported from `ppvm-pauli-word`'s `phase/` module.
//!
//! # Friction: concrete `PhasedPauliWord` alias
//!
//! `word-data-structures.md` sketches the alias as the generic renaming
//! `pub type PhasedPauliWord<W> = Phased<W>;`. This crate instead ships the
//! *concrete* `PhasedPauliWord = Phased<PauliWord>` the task assigns it (the
//! target of the signed-string parsers and the `KeyProduct`-backed product),
//! matching the domain name the old crate exported. The generic wrapper
//! [`Phased`] is still public, so the generic alias's intent — wrapping any
//! word, ordinary or lossy — is preserved; only the exported alias is pinned to
//! the ordinary [`PauliWord`] this component wraps.

mod clifford;
mod data;
mod product;

pub use data::Phased;

/// The phased ordinary Pauli word: `i^φ · g(w)` over a packed [`PauliWord`].
///
/// The concrete domain alias for [`Phased`]`<`[`PauliWord`]`>` — the phased type
/// the signed-string parsers (`"+XYZ"`, `"-XYZ"`, `"+iXYZ"`, `"-iXYZ"`) and the
/// [`KeyProduct`](ppvm_traits_2::KeyProduct)-backed product target. See the
/// crate-level friction note on why this is concrete rather than the design's
/// generic renaming.
pub type PhasedPauliWord = Phased<ppvm_pauli_word_2::PauliWord>;

pub use ppvm_pauli_word_2::PauliWord;
