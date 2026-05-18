// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

//! Extended Stim dialect - interprets PPVM tag-based extensions into a
//! typed AST.

pub mod ast;
mod interpret;
pub mod parser;

pub use ast::{Axis, ExtendedInstruction, ExtendedProgram, RawPassthrough};
pub use parser::{ExtendedParseError, parse_extended};
