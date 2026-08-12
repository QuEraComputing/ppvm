// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! True API blocker: neither old nor new `Sum` or `Term` implements
//! `std::hash::Hash`.
//!
//! Only `Prod` is a hash key and exposes hashing on both crates. Inventing a
//! benchmark that hashes `Display` output would measure string formatting and
//! bytes rather than the requested public symbolic operation, so Sum/Term
//! hashing is deliberately reported as unavailable instead of given a false
//! paired benchmark.
