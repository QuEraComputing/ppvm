// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

//! Public entry point for the extended-dialect parse.

use crate::ast::ParseError;
use crate::extended::ast::ExtendedProgram;
use crate::extended::interpret::interpret;

#[derive(Debug, thiserror::Error)]
pub enum ExtendedParseError {
    #[error(transparent)]
    Parse(#[from] ParseError),
    #[error("invalid tag '{tag}' on '{instruction}' at line {line}: {message}")]
    InvalidTag {
        tag: String,
        instruction: String,
        line: usize,
        message: String,
    },
    #[error("'MPAD' at line {line} bit at position {index} (zero-based) = {value}, must be 0 or 1")]
    InvalidMPadBit {
        line: usize,
        index: usize,
        value: usize,
    },
}

pub fn parse_extended(src: &str) -> Result<ExtendedProgram, ExtendedParseError> {
    // Both `parse_impl` and `interpret` recurse through REPEAT bodies,
    // so run the whole pipeline on the oversized parser stack.
    crate::parser::run_on_parser_stack(|| {
        let prog = crate::parser::parse_impl(src)?;
        interpret(prog)
    })
}
