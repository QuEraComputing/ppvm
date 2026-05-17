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
    #[error("'MPAD' at line {line} target #{index} = {value}, must be 0 or 1")]
    InvalidMPadBit {
        line: usize,
        index: usize,
        value: usize,
    },
}

pub fn parse_extended(src: &str) -> Result<ExtendedProgram, ExtendedParseError> {
    let prog = crate::parser::parse(src)?;
    interpret(prog)
}
