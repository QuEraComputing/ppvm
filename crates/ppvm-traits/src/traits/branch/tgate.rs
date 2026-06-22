// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::config::Config;

/// The non-Clifford `T` gate and its adjoint.
///
/// `T = diag(1, e^{iπ/4})`. Implemented by the simulator backends; see
/// the example in [`ppvm_tableau`](https://docs.rs/ppvm-tableau).
pub trait TGate<T: Config> {
    /// Apply `T` (`diag(1, e^{iπ/4})`) to one qubit.
    fn t(&mut self, addr0: usize);
    /// Apply `T†` to one qubit.
    fn t_dag(&mut self, addr0: usize);

    /// Explicit batched `T`.
    fn t_batch(&mut self, targets: &[usize]) {
        for &q in targets {
            self.t(q);
        }
    }

    /// Explicit batched `T†`.
    fn t_dag_batch(&mut self, targets: &[usize]) {
        for &q in targets {
            self.t_dag(q);
        }
    }
}
