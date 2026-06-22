// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

pub use ppvm_traits::config::Config;

/// Pre-built configs using rustc-hash's `FxHasher` and `f64` coefficients.
///
/// Uses `[u64; N]` Pauli storage, which `bitvec` only supports on 64-bit
/// targets; not available on 32-bit targets such as wasm32. Use the `[u8; N]`
/// [`fxhash`] configs there instead.
#[cfg(target_pointer_width = "64")]
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
