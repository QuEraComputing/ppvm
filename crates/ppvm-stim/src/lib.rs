//! Stim circuit format parser and executor for the `ppvm` simulator.
//!
//! Three-stage pipeline (filled in across Tasks 2–14):
//!
//! 1. `parse` — `&str` → `Program` (pure-Stim AST, tags preserved).
//! 2. `normalize::to_tableau` — `Program` → `TableauProgram` (dialect resolved,
//!    phase-1-unsupported instructions rejected).
//! 3. `execute` — apply a `TableauProgram` to a `GeneralizedTableau`.
//!
//! Multi-shot usage should call `parse` and `normalize::to_tableau` once and
//! call `sample` for the shot loop. The `run_string` / `run_file`
//! convenience helpers re-parse on every call and are intended for
//! single-shot demos only.

pub mod parser;
pub mod tableau_program;
pub mod normalize;
pub mod executor;

pub use parser::ast::{
    AnnotationKind, GateName, MeasureName, NoiseName, ParseError, Program,
    RawInstruction, Tag, TagParam,
};
pub use parser::parse;

pub use tableau_program::{
    GateKind, Instruction, MeasureKind, NoiseKind, TableauProgram,
};
