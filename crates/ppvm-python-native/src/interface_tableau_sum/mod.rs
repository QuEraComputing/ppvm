// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bnum::types::{U256, U512, U1024, U2048};
use paste::paste;
#[cfg(feature = "legacy")]
use ppvm_pauli_sum::prelude::*;
#[cfg(feature = "traits-2")]
use ppvm_tableau_2::prelude::*;
use pyo3::prelude::*;

use crate::interface_tableau::measurement_to_u8;

#[macro_use]
mod gates;
#[macro_use]
mod macros;
#[macro_use]
mod noise;
#[macro_use]
mod sampler;
#[macro_use]
mod state;

// up to 64 qubits
create_sum_interface_range!(usize, 1);

// 64 - 128 qubits
create_sum_interface_range!(u128, 2);

// 128 - 256 qubits
create_sum_interface_range!(U256, 3, 4);

create_sum_interface_range!(U512, 5, 6, 7, 8);

create_sum_interface_range!(U1024, 9, 10, 11, 12, 13, 14, 15, 16);

create_sum_interface_range!(
    U2048, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32
);
