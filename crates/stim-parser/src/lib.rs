// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

pub mod ast;
mod display;
pub mod extended;
mod grammar;
mod line_map;
mod parser;
mod table;

use line_map::LineMap;

pub mod prelude {
    pub use crate::ast::{
        AnnotationKind, GateName, MeasureName, NoiseName, ParseError, Program, RawInstruction, Tag,
        TagParam, Target,
    };
    pub use crate::extended::{
        Axis, ExtendedInstruction, ExtendedParseError, ExtendedProgram, RawPassthrough,
        parse_extended,
    };
    pub use crate::parser::parse;
}
