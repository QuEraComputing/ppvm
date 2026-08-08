// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Specialized classical mixtures of complete generalized-tableau states.
//!
//! Approximate amplitude equality is deliberately not exposed as `Eq`/`Hash`:
//! it is non-transitive. A structural fingerprint selects candidates and a full
//! collision check decides whether probabilities may be coalesced.

mod data;
mod equality;
mod fingerprint;
mod gates;
mod measure;
mod noise;
mod sampler;
#[cfg(test)]
mod tests;

pub use data::{GeneralizedTableauMixture, GeneralizedTableauSum};
pub use sampler::MixtureSampler;

pub(crate) use data::{Branch, LazyBranch};
