// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Packed Pauli words, explicit phases, and hashing utilities.

mod data;
mod hash;
mod pauli_parser;
mod phase;
mod storage;
mod product;

pub use data::PauliWord;
pub use hash::HashFinalize;
pub use phase::Phased;
pub use storage::{DefaultStorage, PauliStorage};
