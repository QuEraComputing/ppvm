// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use paste::paste;
#[cfg(feature = "legacy")]
use ppvm_pauli_sum::prelude::*;
#[cfg(feature = "legacy")]
use ppvm_pauli_sum::strategy::{
    CoefficientThreshold, CombinedStrategy, MaxLossWeight, MaxPauliWeight,
};
#[cfg(feature = "traits-2")]
use ppvm_pauli_sum_2::*;
#[cfg(feature = "traits-2")]
use ppvm_traits_2::prelude::*;
use pyo3::prelude::*;

#[macro_use]
mod backend;
#[macro_use]
mod gates;
#[macro_use]
mod loss;
#[macro_use]
mod macros;
#[macro_use]
mod noise;
#[macro_use]
mod python;
#[macro_use]
mod rotations;
#[macro_use]
mod state;

create_interface_range!(
    IndexMapFxHash,
    false,
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15
);

create_interface_range!(
    IndexMapFxHash,
    true,
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15
);
