// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `ppvm-pauli-sum-2`: the graded sparse-sum engine.
//!
//! `Sum<S, P>` is a sparse linear combination over the graded algebra `C[K]` —
//! the free `C`-module on the keys `K` — with a storage backend `S` (the graded
//! traits `Support`/`Accumulate`/`Scale`/`Pair`/`Multiply` are impl'd directly on
//! `Vec<(K,C)>` and `HashMap<K, C, IdentityBuildHasher>`) and a truncation
//! `Policy` `P`. Gates are `TermProducer`s fed into `apply` →
//! `accumulate_batch` → `reduce` → `policy.truncate`. `PauliSum = Sum<HashMapStore
//! <PauliWord, C>, P>` is the domain alias.
//!
//! Trait contracts and their Lean-validated semantics (the `C[K]` module/algebra
//! laws in `GradedMap.lean`, truncation bounds in `Truncation.lean`) follow
//! [`docs/design/traits-2-configuration-and-hashing.md`] and `lean/PPVM/**`; the
//! algorithm is ported from `ppvm-pauli-sum` to hold the hot paths at parity.
//!
//! Phase 3: the `Sum` core + graded traits + Clifford propagation, plus the
//! non-Clifford [`RotationOne`](ppvm_traits_2::RotationOne) branch and the
//! diagonal [`PauliError`](ppvm_traits_2::PauliError) channel; the L4 `Multiply`
//! product follows.
//!
//! # Deferrals (this component)
//!
//! - **L4 `Multiply`** (the twisted operator product) — later component.
//! - **Columnar `ColumnStore` (SoA) backend** — Phase 6. `accumulate_batch`
//!   takes an array-of-structs `TermBatch` here.
//! - The graded-algebra **container impls** (`Support`/`Accumulate`/`Scale`/
//!   `Pair`/`Retain` on `Vec`/`HashMap`) live in `ppvm-traits-2` (orphan rule);
//!   see `ppvm_traits_2::containers`.

use std::collections::HashMap;

mod clifford;
mod noise;
mod policy;
mod producer;
mod rotation;
mod store;
mod sum;

pub use policy::{CoefficientThreshold, CombinedPolicy, MaxPauliWeight, NoPolicy, Policy};
pub use producer::RekeyProducer;
pub use store::StoreAlloc;
pub use sum::Sum;

// Re-exports so downstream can name the storage contract without depending on
// `ppvm-traits-2` directly.
pub use ppvm_pauli_word_2::PauliWord;
pub use ppvm_traits_2::IdentityBuildHasher;

/// The provided hash-map storage: a `HashMap` whose hasher is the pass-through
/// [`IdentityBuildHasher`], so a key's finalized `key_hash()` digest reaches
/// hashbrown untouched (Design: §"The pass-through storage contract").
pub type HashMapStore<K, C> = HashMap<K, C, IdentityBuildHasher>;

/// The Pauli-propagation sum: `C[PauliWord]` over the pass-through hash-map
/// storage. `C` defaults to `f64` and the policy to [`NoPolicy`].
///
/// Design: §"The convenience bundle" (`PauliSum<C, P> = Sum<HashMapStore<PauliWord,
/// C>, P>`).
pub type PauliSum<C = f64, P = NoPolicy> = Sum<HashMapStore<PauliWord, C>, P>;
