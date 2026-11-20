// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Symbolic coefficients for `ppvm-traits-2`: polynomials in `sin`/`cos` of
//! circuit parameters, evaluated once angle values are known.

#[macro_use]
mod macros;

mod add;
mod coeff;
mod display;
mod eval;
mod mul;
mod prod;
mod sum;
mod term;

pub use eval::MissingVariable;
pub use term::{DEFAULT_MAX_SIN, DEFAULT_MIN_EPS, Term, Var};
