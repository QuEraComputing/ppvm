// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Hybrid stabilizer simulation: `|ψ⟩ = C |χ⟩`, an inverse tableau for the
//! Clifford `C` plus the state `|χ⟩` of the virtual qubits in an exchangeable
//! representation. Design: `docs/design/hybrid-tableau.md`.
//!
//! - `inverse_tableau`: [`InverseTableau`] and its [`VirtualView`].
//! - `virtual_state`: the [`VirtualState`] trait and [`ActiveDormant`], which
//!   keeps amplitudes only for the active virtual qubits.
//! - `generalized_tableau`: [`GeneralizedTableau`], the simulator.

mod generalized_tableau;
mod inverse_tableau;
mod virtual_state;

pub use generalized_tableau::GeneralizedTableau;
pub use inverse_tableau::{InverseTableau, OwnedRow, RowRef, VirtualView};
pub use virtual_state::{
    ActiveDormant, ActiveState, AxisPauli, DenseState, Rotation, VirtualState,
};
