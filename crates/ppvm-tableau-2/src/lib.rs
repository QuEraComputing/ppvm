// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Hybrid stabilizer simulation: an inverse tableau plus a dense vector over
//! the active qubits. Design: `docs/design/hybrid-tableau.md`.

mod inverse_tableau;

pub use inverse_tableau::{InverseTableau, RowRef};
