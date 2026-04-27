pub mod ast;

use ast::{ParseError, Program};

/// Parse Stim source into a [`Program`]. Performs name, arg-count, and
/// target-arity validation; preserves tags verbatim.
pub fn parse(_src: &str) -> Result<Program, ParseError> {
    todo!("implemented in Task 4 onward")
}
