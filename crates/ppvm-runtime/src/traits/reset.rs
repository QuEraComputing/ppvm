// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// Reset a qubit to the `|0⟩` computational-basis state.
pub trait Reset {
    /// Reset qubit `addr0` to `|0⟩`.
    fn reset(&mut self, addr0: usize);
}
