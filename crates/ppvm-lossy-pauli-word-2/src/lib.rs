// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `ppvm-lossy-pauli-word-2`: the concrete **lossy** Pauli word representation.
//!
//! [`LossyPauliWord`] is another concrete implementation of the Pauli-word
//! concept (a sibling of `ppvm-pauli-word-2`'s `PauliWord`), kept in its own
//! crate. It adds a packed **loss plane** to the X/Z planes and implements
//! [`Word`](ppvm_traits_2::Word)`<Site = LossySite<Pauli>>`,
//! [`PauliBits`](ppvm_traits_2::PauliBits) (with `is_lost` overridden),
//! [`SymplecticColumns`](ppvm_traits_2::SymplecticColumns) +
//! [`PhaseTrack`](ppvm_traits_2::PhaseTrack) (phase-discarding ⇒ bit-only
//! [`Clifford`](ppvm_traits_2::Clifford)),
//! [`Indexable`](ppvm_traits_2::Indexable), and
//! [`Columnar`](ppvm_traits_2::Columnar). Loss *writes* (`set_lost`,
//! `clear_loss`) and `loss_weight()` are **inherent** methods per the design,
//! not part of any trait.
//!
//! Concrete layout follows [`docs/design/word-data-structures.md`] §"Lossy Pauli
//! word"; trait contracts and Lean-validated semantics follow
//! [`docs/design/traits-2-configuration-and-hashing.md`] and `lean/PPVM/**`. The
//! packed-storage/hash utilities [`PauliStorage`](ppvm_pauli_word_2::PauliStorage)
//! and [`HashFinalize`](ppvm_pauli_word_2::HashFinalize) are reused from
//! `ppvm-pauli-word-2`.
//!
//! The backing fields (packed X/Z/loss planes, one lazy finalized-digest
//! `AtomicU64`, and the representation parameters `A`/`H`) are **private**; behavior
//! is exposed only through the implemented traits and a small set of inherent
//! constructors/accessors.

mod clifford;
mod column;
mod data;
mod hash;

pub use column::LossyPauliKeyColumn;
pub use data::LossyPauliWord;
