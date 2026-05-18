// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

use std::hash::BuildHasher;

use crate::traits::{ACMap, Coefficient, PauliStorage, PauliWordTrait, Strategy};

/// Compile-time configuration bundle for a [`PauliSum`](crate::sum::PauliSum).
///
/// `Config` is a zero-cost generic dispatch mechanism: it pins down the
/// storage backing (`Storage`), coefficient type (`Coeff`), truncation
/// strategy (`Strategy`), hasher (`BuildHasher`), Pauli-word representation
/// (`PauliWordType`), and concrete map (`Map`) used by a single
/// `PauliSum` instantiation. Pre-built bundles live in the submodules
/// ([`fxhash`], [`indexmap`], …); user code may define its own.
pub trait Config: Clone {
    /// Backing storage for an individual Pauli word — typically a `[u8; N]`
    /// where `N` is the number of bytes packed.
    type Storage: PauliStorage;
    /// Numeric coefficient type (e.g. `f64`, `Complex<f64>`).
    type Coeff: Coefficient;
    /// Truncation strategy applied when
    /// [`PauliSum::truncate`](crate::sum::PauliSum::truncate) runs.
    type Strategy: Strategy;
    /// Hasher used by the underlying map.
    type BuildHasher: BuildHasher + Clone + Default;
    /// Concrete [`PauliWordTrait`] implementation used as map keys.
    type PauliWordType: PauliWordTrait;
    /// Concrete map type satisfying [`ACMap`] over the chosen
    /// storage / coefficient / hasher / word combination.
    type Map: ACMap<Self::Storage, Self::Coeff, Self::BuildHasher, Self::PauliWordType>;
}

/// Pre-built configs using rustc-hash's `FxHasher` and `f64` coefficients.
pub mod fx64hash;
/// Pre-built configs using `FxHasher`.
pub mod fxhash;

/// Pre-built configs backed by `dashmap::DashMap` for concurrent access.
#[cfg(feature = "dashmap")]
pub mod dashmap;

/// Pre-built configs backed by `indexmap::IndexMap`, preserving
/// insertion order — useful for deterministic iteration and snapshot
/// testing.
#[cfg(feature = "indexmap")]
pub mod indexmap;

/// Pre-built configs using `gxhash` — fast on platforms with AES
/// hardware acceleration. Requires the `gxhash` feature.
#[cfg(all(
    feature = "gxhash",
    any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")
))]
pub mod gxhash;
