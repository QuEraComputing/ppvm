// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Foundational trait system, `Config` trait, single-qubit `Pauli` alphabet,
//! and the `ACMap`/`Trace` map implementations shared across ppvm crates.
pub mod char;
pub mod config;
mod map;
pub mod traits;

pub use char::Pauli;
pub use config::Config;
pub use traits::*;

pub mod prelude {
    pub use crate::char::Pauli;
    pub use crate::config::Config;
    pub use crate::traits::*;
}
