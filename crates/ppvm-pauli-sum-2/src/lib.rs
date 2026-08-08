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
//! Phase 6: the columnar [`ColumnStore`] (structure-of-arrays) backend — the
//! same engine, the same behaviour, over parallel key-plane / coefficient /
//! digest columns instead of a `HashMap`. [`ColumnPauliSum`] is its domain alias.
//!
//! # Deferrals (this component)
//!
//! - `accumulate_batch` still takes an array-of-structs-ish `TermBatch` (a
//!   scalar `Vec<W>` key column plus a coefficient column) rather than the
//!   design's `W::Column`; see the friction note in `ppvm_traits_2::batch`. The
//!   SoA key column *is* used by [`ColumnStore`]'s own support.
//! - **`ColumnStore` plane granularity/alignment** (design open question 2,
//!   which the implementation plan defers "to Phase 6 benchmarks"). The key
//!   column keeps each key's X/Z blob as one contiguous `A` slot per plane —
//!   `ppvm-pauli-word-2`'s Phase-2 `PauliKeyColumn`, consumed as-is. A
//!   bit-sliced layout (qubit-major planes across keys) would make the weight
//!   predicate and the plane compare vectorizable too; it is not needed for the
//!   wins this phase measured (which are on the *coefficient* column) and it
//!   would change `KeyColumn`'s layout contract, so it stays open.
//! - `Pair::probe_batch` consumes the batch's precomputed digest column but
//!   probes **scalar**; the design's "coalesced gathers" (group prefetch over a
//!   batch of buckets) are not implemented. This is the same side on which the
//!   backend currently *loses* (see the [`ColumnStore`] measurement table), so
//!   it is the obvious next lever.
//! - The graded-algebra **container impls** (`Support`/`Accumulate`/`Scale`/
//!   `Pair`/`Multiply`/`Retain` on `Vec`/`HashMap`) live in `ppvm-traits-2`
//!   (orphan rule); see `ppvm_traits_2::containers`.

mod clifford;
mod column_store;
mod display;
mod indexmap_store;
mod loss;
pub mod multiply;
mod noise;
mod ops;
mod pattern;
mod policy;
mod producer;
mod proj;
mod rotation;
mod store;
mod sum;
mod trace;

pub use column_store::ColumnStore;
pub use indexmap_store::IndexMapStore;
pub use pattern::{EnumMatchesPauliPattern, PatternParseError, PatternSite, PauliPattern, SiteSet};
pub use policy::{
    CoefficientThreshold, CombinedPolicy, MaxLossWeight, MaxPauliWeight, NoPolicy, Policy,
};
pub use producer::RekeyProducer;
pub use store::{
    AddTerm, ApplyProducer, BranchInPlace, HashMapStore, InsertTerm, MultiplyInPlace, StoreAlloc,
};
pub use sum::{PreserveSet, Sum};

// Re-exports so downstream can name the storage contract without depending on
// `ppvm-traits-2` directly.
pub use ppvm_lossy_pauli_word_2::LossyPauliWord;
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

/// An insertion-ordered Pauli sum with explicit fixed-width byte storage.
///
/// This is the storage shape used by the Python ABI's width-specialized sums;
/// [`PauliSum`] remains the default hash-map alias.
pub type IndexPauliSum<const N: usize, C = f64, P = NoPolicy> =
    Sum<IndexMapStore<PauliWord<[u8; N]>, C>, P>;

/// The Pauli-propagation sum over the **columnar** (structure-of-arrays)
/// storage: key planes, coefficients and digests in parallel columns behind an
/// open-addressed index.
///
/// A drop-in swap for [`PauliSum`] — observationally identical (Design:
/// §"Backends are containers; columnar is expressible from day one"; the
/// implementation plan's Phase 6 gate is "a backend swap must be
/// observationally identical"). See [`ColumnStore`] for the layout and the
/// measured op-class trade-off.
pub type ColumnPauliSum<C = f64, P = NoPolicy> = Sum<ColumnStore<PauliWord, C>, P>;

/// The **lossy** Pauli-propagation sum: `C[LossyPauliWord]` over the same
/// hash-join storage — the neutral-atom configuration, where a site may be
/// `Lost` in addition to carrying a Pauli.
///
/// Old's `LossyPauliSum` is spelled as a `Config` whose `PauliWordType` is
/// `LossyPauliWord<[u8; N], FxBuildHasher>` (`ppvm-pauli-sum/tests/loss.rs`);
/// with `-2`'s storage-carries-the-key shape it is a plain type alias. Only this
/// key admits [`ResetLossChannel`](ppvm_traits_2::ResetLossChannel) and only here
/// does the [`MaxLossWeight`] policy bite — for [`PauliSum`] the loss branches are
/// dead code and `loss_weight()` is a const `0` (architecture feature 11).
pub type LossyPauliSum<C = f64, P = NoPolicy> = Sum<HashMapStore<LossyPauliWord, C>, P>;
