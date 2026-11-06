// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// Reset the state of a live qubit to `|0⟩`, leaving lost and leaked qubits
/// untouched.
pub trait ResetState {
    /// Reset the qubit at `addr0` to `|0⟩` if it is live; no-op otherwise.
    fn reset_state(&mut self, addr0: usize);
}

/// Full reset to `|0⟩` of a single Pauli word, in the Heisenberg picture.
///
/// Pauli-propagation backends build [`Reset`] from this: each term maps to
/// the words returned here, each carrying the term's coefficient.
pub trait WordReset: crate::traits::PauliWordTrait {
    /// Words this word maps to under a full reset of qubit `addr0`.
    ///
    /// `I` and `Z` on `addr0` map to the identity on the full qubit space
    /// (for lossy words, `I + L`); `X`, `Y` and `L` drop out.
    fn word_reset(&self, addr0: usize) -> impl Iterator<Item = Self> + use<Self>;
}

/// Reset one qubit to a computational/Pauli basis state.
pub trait Reset: ResetState + crate::traits::Clifford + crate::traits::CliffordExtensions {
    /// Reset one qubit to `|0⟩` (stim `R`/`RZ`).
    ///
    /// This is a full reset: a lost or leaked qubit is first brought back
    /// to the computational subspace, then its state is reset via
    /// [`ResetState::reset_state`].
    fn reset(&mut self, addr0: usize);
    /// stim `RZ` alias — reset to `|0⟩`.
    fn reset_z(&mut self, addr0: usize) {
        self.reset(addr0)
    }
    /// stim `RX` — reset to `|+⟩`.
    fn reset_x(&mut self, addr0: usize) {
        self.reset(addr0);
        self.h(addr0);
    }
    /// stim `RY` — reset to `|i⟩`.
    fn reset_y(&mut self, addr0: usize) {
        self.reset(addr0);
        self.h(addr0);
        self.s(addr0);
    }

    /// Explicit batched reset to `|0⟩`.
    fn reset_many(&mut self, targets: &[usize]) {
        for &q in targets {
            self.reset(q);
        }
    }

    /// Explicit batched `RZ` alias.
    fn reset_z_many(&mut self, targets: &[usize]) {
        self.reset_many(targets)
    }

    /// Explicit batched `RX`.
    fn reset_x_many(&mut self, targets: &[usize]) {
        for &q in targets {
            self.reset_x(q);
        }
    }

    /// Explicit batched `RY`.
    fn reset_y_many(&mut self, targets: &[usize]) {
        for &q in targets {
            self.reset_y(q);
        }
    }
}
