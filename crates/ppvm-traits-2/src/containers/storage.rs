// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// Implemented by storage backends, independently of coefficient arithmetic.
/// Both traversals must preserve the same keys and coefficient values.
pub trait RekeyStrategy {
    /// Prefer draining entries and moving coefficients when true.
    /// When false, prefer borrowing entries and cloning coefficients.
    const PREFER_MOVED_REKEY: bool = false;
}
