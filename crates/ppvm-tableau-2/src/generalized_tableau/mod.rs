// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`GeneralizedTableau`]: the simulator, `|ψ⟩ = C |χ⟩`.
//!
//! - `data`: the struct, constructors and shared helpers.
//! - `gates`: Clifford, T, rotation and `U3` gates, and reset.
//! - `measure`: Z measurement.
//! - `loss`: loss and leakage channels.
//! - `noise`: Pauli and depolarizing channels.

mod data;
mod gates;
mod loss;
mod measure;
mod noise;
mod status;

pub use data::GeneralizedTableau;
