// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bnum::types::{U256, U512, U1024, U2048};
use paste::paste;
#[cfg(feature = "traits-2")]
use ppvm_pauli_sum_2::PauliPattern;
#[cfg(feature = "legacy")]
use ppvm_tableau::prelude::*;
#[cfg(feature = "traits-2")]
use ppvm_tableau_2::prelude::*;
use pyo3::prelude::*;
use pyo3::types::{PyComplex, PyDict};

#[macro_use]
mod gates;
#[macro_use]
mod macros;
#[macro_use]
mod noise;
#[macro_use]
mod state;
#[macro_use]
mod stim;

pub(crate) fn measurement_to_u8(m: Option<bool>) -> u8 {
    match m {
        Some(false) => 0,
        Some(true) => 1,
        None => 2,
    }
}

// up to 64 qubits
create_interface_range!(IndexMapFxHash, usize, 1);

// 64 - 128 qubits
create_interface_range!(IndexMapFxHash, u128, 2);

// 128 - 256 qubits
create_interface_range!(IndexMapFxHash, U256, 3, 4);

create_interface_range!(IndexMapFxHash, U512, 5, 6, 7, 8);

create_interface_range!(IndexMapFxHash, U1024, 9, 10, 11, 12, 13, 14, 15, 16);

create_interface_range!(
    IndexMapFxHash,
    U2048,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    32
);
