// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

pub use ppvm_traits::config::Config;

/// Pre-built configs using rustc-hash's `FxHasher` and `f64` coefficients.
pub mod fx64hash;
/// Pre-built configs using `FxHasher`.
pub mod fxhash;

/// Pre-built configs backed by `dashmap::DashMap` for concurrent access.
/// Native-only: `dashmap` pulls `rayon`, which needs OS threads, so this
/// module is unavailable on `wasm32`.
#[cfg(all(feature = "dashmap", not(target_arch = "wasm32")))]
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
