// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::config::Config;

/// The general single-qubit `U3(θ, φ, λ)` gate.
pub trait U3Gate<T: Config> {
    /// Apply `U3(θ, φ, λ)` to qubit `addr`.
    fn u3(&mut self, addr: usize, theta: T::Coeff, phi: T::Coeff, lambda: T::Coeff);
}
