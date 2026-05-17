//! Extended Stim dialect - interprets PPVM tag-based extensions into a
//! typed AST.

pub mod ast;
mod interpret;
pub mod parser;

pub use ast::{Axis, ExtendedInstruction, ExtendedProgram};
pub use parser::{ExtendedParseError, parse_extended};
