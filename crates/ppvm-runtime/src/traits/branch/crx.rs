// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

use crate::config::Config;

/// Controlled `RX` rotation.
pub trait CRx<T: Config> {
    /// Apply `CRX(θ)` with the given control and target.
    fn crx(&mut self, control: usize, target: usize, theta: T::Coeff);
}
