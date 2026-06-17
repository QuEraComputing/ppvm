// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::config::Config;

/// The non-Clifford `T` gate and its adjoint.
///
/// `T = diag(1, e^{iπ/4})`. Implemented by the simulator backends; see
/// the example in [`ppvm_tableau`](https://docs.rs/ppvm-tableau).
pub trait TGate<T: Config> {
    /// Apply `T` (`diag(1, e^{iπ/4})`) to each target.
    fn t(&mut self, targets: impl crate::traits::Targets);
    /// Apply `T†` to each target.
    fn t_dag(&mut self, targets: impl crate::traits::Targets);
}
