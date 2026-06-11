// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::BuildHasher;

use crate::traits::{ACMap, Coefficient, PauliStorage, PauliWordTrait, Strategy};

/// Compile-time configuration bundle for a `PauliSum`.
///
/// `Config` is a zero-cost generic dispatch mechanism: it pins down the
/// storage backing (`Storage`), coefficient type (`Coeff`), truncation
/// strategy (`Strategy`), hasher (`BuildHasher`), Pauli-word representation
/// (`PauliWordType`), and concrete map (`Map`) used by a single
/// `PauliSum` instantiation. Pre-built bundles live in the submodules
/// ([`fxhash`], [`fx64hash`], plus the optional-feature bundles in
/// `ppvm_paulisum::config`); user code may define its own.
pub trait Config: Clone {
    /// Backing storage for an individual Pauli word — typically a `[u8; N]`
    /// where `N` is the number of bytes packed.
    type Storage: PauliStorage;
    /// Numeric coefficient type (e.g. `f64`, `Complex<f64>`).
    type Coeff: Coefficient;
    /// Truncation strategy applied when `PauliSum::truncate` runs.
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
