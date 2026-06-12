// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`Config`] bundles for `PauliSum`.
//!
//! Re-exports the [`Config`] trait and the always-on bundles
//! ([`fxhash`], [`fx64hash`]) from `ppvm-runtime`, plus the
//! feature-gated bundles that live here ([`indexmap`], [`gxhash`],
//! [`dashmap`]).

pub use ppvm_runtime::config::{Config, fx64hash, fxhash};

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
