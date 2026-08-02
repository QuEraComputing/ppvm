// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `ppvm-pauli-sum-2`: the graded sparse-sum engine.
//!
//! `Sum<S, P>` is a sparse linear combination over the graded algebra `C[K]` —
//! the free `C`-module on the keys `K` — with a storage backend `S` (the graded
//! traits `Support`/`Accumulate`/`Scale`/`Pair`/`Multiply` are impl'd directly on
//! `Vec<(K,C)>` and `HashMap<K, C, IdentityBuildHasher>`) and a truncation
//! `Policy` `P`.
//!
//! A gate is a `TermProducer` fed into [`Sum::apply`] → `reset` →
//! `accumulate_batch`; the four hot gate families (Clifford re-key, rotation,
//! per-key scale, per-key sign flip) instead drive the storage's in-place fast
//! paths, which produce and merge in one pass over the support. **No gate path
//! runs `reduce` and no gate path truncates**: both are caller-driven
//! ([`Sum::reduce`], [`Sum::truncate`]), which is the old crate's contract — a
//! zero-coefficient term survives every operation, and two sub-threshold
//! contributions to one key may merge before any truncation sees them.
//! `PauliSum = Sum<HashMapStore<PauliWord, C>, P>` is the domain alias.
//!
//! Trait contracts and their Lean-validated semantics (the `C[K]` module/algebra
//! laws in `GradedMap.lean`, truncation bounds in `Truncation.lean`) follow
//! [`docs/design/traits-2-configuration-and-hashing.md`] and `lean/PPVM/**`; the
//! algorithm is ported from `ppvm-pauli-sum` to hold the hot paths at parity.
//!
//! Phase 3: the `Sum` core + graded traits + Clifford propagation, the
//! non-Clifford [`RotationOne`](ppvm_traits_2::RotationOne) branch, the diagonal
//! [`PauliError`](ppvm_traits_2::PauliError) channel, and the L4
//! [`Multiply`](ppvm_traits_2::Multiply) operator product (the twisted
//! convolution — see [`multiply`]).
//!
//! # Deferrals (this component)
//!
//! - **Columnar `ColumnStore` (SoA) backend** — Phase 6. `accumulate_batch`
//!   takes an array-of-structs `TermBatch` here.
//! - The graded-algebra **container impls** (`Support`/`Accumulate`/`Scale`/
//!   `Pair`/`Multiply`/`Retain` on `Vec`/`HashMap`) live in `ppvm-traits-2`
//!   (orphan rule); see `ppvm_traits_2::containers`.

mod clifford;
pub mod multiply;
mod noise;
mod ops;
mod policy;
mod producer;
mod rotation;
mod store;
mod sum;
mod trace;

pub use policy::{CoefficientThreshold, CombinedPolicy, MaxPauliWeight, NoPolicy, Policy};
pub use producer::RekeyProducer;
pub use store::{AddTerm, ApplyProducer, HashMapStore, MultiplyInPlace, StoreAlloc};
pub use sum::{PreserveSet, Sum};
pub use trace::{PauliPattern, SiteSet};

// Re-exports so downstream can name the storage contract without depending on
// `ppvm-traits-2` directly.
pub use ppvm_pauli_word_2::PauliWord;
pub use ppvm_traits_2::IdentityBuildHasher;

/// Support module for [`impl_scalar_mul!`] — the `ppvm-traits-2` bounds its
/// expansion mentions, re-exported so the macro's `$crate::`-absolute paths
/// resolve in a downstream crate that does not itself depend on `ppvm-traits-2`.
///
/// Not a stability surface: name these traits through `ppvm_traits_2` instead.
#[doc(hidden)]
pub mod reexport {
    pub use ppvm_traits_2::{Accumulate, Indexable, Scale, Word};
}

/// The Pauli-propagation sum: `C[PauliWord]` over the pass-through hash-map
/// storage. `C` defaults to `f64` and the policy to [`NoPolicy`].
///
/// Design: §"The convenience bundle" (`PauliSum<C, P> = Sum<HashMapStore<PauliWord,
/// C>, P>`).
pub type PauliSum<C = f64, P = NoPolicy> = Sum<HashMapStore<PauliWord, C>, P>;
