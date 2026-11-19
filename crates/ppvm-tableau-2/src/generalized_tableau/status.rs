// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// A qubit relative to the computational subspace. Lost and leaked qubits are
/// inactive: gates skip them. Measuring a lost qubit gives `None`, a leaked
/// one its pinned bit.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) enum QubitStatus {
    #[default]
    Live,
    Lost,
    Leaked,
}

impl QubitStatus {
    pub(super) fn is_inactive(self) -> bool {
        self != QubitStatus::Live
    }
}
