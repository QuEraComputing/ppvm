// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `ppvm-pauli-word-2`: concrete packed Pauli word representations for the
//! second trait-system experiment.
//!
//! This crate hosts the concrete types the trait foundation in `ppvm-traits-2`
//! abstracts over. This component ships the ordinary [`PauliWord`] (packed X/Z
//! planes with a lazy `AtomicU64` structural hash); `LossyPauliWord` (adds a
//! loss plane) and the `Phased` wrapper are later components.
//!
//! `PauliWord` implements the trait surface the design assigns it:
//! [`Word`](ppvm_traits_2::Word)`<Site = Pauli>`,
//! [`PauliBits`](ppvm_traits_2::PauliBits),
//! [`SymplecticColumns`](ppvm_traits_2::SymplecticColumns) +
//! [`PhaseTrack`](ppvm_traits_2::PhaseTrack) (hence the blanket
//! [`Clifford`](ppvm_traits_2::Clifford)),
//! [`KeyProduct`](ppvm_traits_2::KeyProduct) (the twisted product
//! `v·w = iᵏ (v⊕w)`), [`Indexable`](ppvm_traits_2::Indexable), and
//! [`Columnar`](ppvm_traits_2::Columnar).
//!
//! The packed layout follows [`docs/design/word-data-structures.md`]; the phase
//! and Clifford kernels are ported from `ppvm-pauli-word` to keep hot paths at
//! parity; the trait contracts and their Lean-validated semantics follow
//! [`docs/design/traits-2-configuration-and-hashing.md`] and `lean/PPVM/**`.
//!
//! The backing fields (packed planes, hash cache, representation parameters `A`
//! and `H`) are **private**; behavior is exposed only through the implemented
//! traits and a small set of inherent constructors/accessors.

mod clifford;
mod column;
mod data;
mod hash;
mod product;
mod storage;

pub use column::PauliKeyColumn;
pub use data::PauliWord;
pub use storage::{DefaultStorage, HashFinalize, PauliStorage};
