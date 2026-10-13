// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::complex::Complex;

use crate::data::GeneralizedTableau;
use crate::sparsevec::SparseVector;
use ppvm_traits::config::Config;

/// Per-qubit status relative to the computational subspace.
///
/// Loss and leakage both take a qubit out of the computational subspace, so
/// gates skip any status other than [`QubitStatus::Live`]. They differ only
/// at measurement: a lost qubit reports `None`, a leaked qubit reports the
/// pinned computational bit. Loss overwrites leakage — a leaked qubit that is
/// later lost is [`QubitStatus::Lost`], indistinguishable from a qubit lost
/// from the computational subspace.
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum QubitStatus {
    /// Qubit is in the computational subspace; gates apply normally.
    #[default]
    Live = 0,
    /// Qubit is pinned to a computational basis state. Gates skip it;
    /// measurement returns the pinned `0`/`1`.
    Leaked = 1,
    /// Qubit has been lost. Gates skip it; measurement returns `None`.
    Lost = 2,
}

impl QubitStatus {
    /// Outside the computational subspace — lost or leaked.
    #[inline]
    pub const fn is_inactive(self) -> bool {
        (self as u8) != 0
    }

    #[inline]
    pub const fn is_lost(self) -> bool {
        matches!(self, QubitStatus::Lost)
    }

    #[inline]
    pub const fn is_leaked(self) -> bool {
        matches!(self, QubitStatus::Leaked)
    }
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> GeneralizedTableau<T, I, C> {
    /// Whether qubit `addr0` is outside the computational subspace — either lost
    /// or leaked. Gates skip such qubits.
    #[inline]
    pub fn is_inactive(&self, addr0: usize) -> bool {
        self.qubit_status[addr0].is_inactive()
    }

    /// Whether qubit `addr0` is lost. Measurement of a lost qubit returns `None`.
    #[inline]
    pub fn is_lost(&self, addr0: usize) -> bool {
        self.qubit_status[addr0].is_lost()
    }

    /// Whether qubit `addr0` is leaked. Measurement of a leaked qubit returns
    /// the pinned computational bit.
    #[inline]
    pub fn is_leaked(&self, addr0: usize) -> bool {
        self.qubit_status[addr0].is_leaked()
    }

    /// Whether any qubit in `indices` is lost or leaked.
    #[inline]
    pub(crate) fn any_inactive(&self, indices: &[usize]) -> bool {
        indices.iter().any(|&i| self.is_inactive(i))
    }

    /// Whether any pair has a lost or leaked control or target.
    #[inline]
    pub(crate) fn any_inactive_pair(&self, pairs: &[(usize, usize)]) -> bool {
        pairs
            .iter()
            .any(|&(c, t)| self.is_inactive(c) || self.is_inactive(t))
    }
}
