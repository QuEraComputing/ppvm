// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use smallvec::SmallVec;
use vihaco::Message;

vihaco::component! {
    #[derive(Debug, Default)]
    pub component Circuit {}

    instruction {
        TwoQubitPauliError,
        Truncate,
        Trace,
        X,
        Y,
        Z,
        H,
        SqrtXAdj,
        SqrtX,
        SqrtYAdj,
        SqrtY,
        SAdj,
        S,
        CNOT,
        CZ,
        TAdj,
        T,
        RXX,
        RYY,
        RZZ,
        RX,
        RY,
        RZ,
        U3,
        Measure,
        Reset,
        R,
        Loss,
        CorrelatedLoss,
        PauliError,
        Depolarize2,
        Depolarize,
    }
}

pub use circuit::{runtime, syntax};
pub use runtime::Instruction as CircuitInstruction;

pub type CircuitSurfaceInstruction = syntax::Instruction;

impl std::fmt::Display for runtime::Instruction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use runtime::Instruction::*;
        match self {
            TwoQubitPauliError => write!(f, "TwoQubitPauliError"),
            Truncate => write!(f, "Truncate"),
            Trace => write!(f, "Trace"),

            X => write!(f, "X"),
            Y => write!(f, "Y"),
            Z => write!(f, "Z"),
            H => write!(f, "H"),

            SqrtXAdj => write!(f, "SqrtXAdj"),

            SqrtX => write!(f, "SqrtX"),

            SqrtYAdj => write!(f, "SqrtYAdj"),

            SqrtY => write!(f, "SqrtY"),

            SAdj => write!(f, "SAdj"),
            S => write!(f, "S"),

            CNOT => write!(f, "CNOT"),
            CZ => write!(f, "CZ"),

            TAdj => write!(f, "TAdj"),
            T => write!(f, "T"),

            RXX => write!(f, "RXX"),
            RYY => write!(f, "RYY"),
            RZZ => write!(f, "RZZ"),

            RX => write!(f, "RX"),
            RY => write!(f, "RY"),
            RZ => write!(f, "RZ"),

            U3 => write!(f, "U3"),

            Measure => write!(f, "Measure"),
            Reset => write!(f, "Reset"),

            R => write!(f, "R"),

            Loss => write!(f, "Loss"),
            CorrelatedLoss => write!(f, "CorrelatedLoss"),

            PauliError => write!(f, "PauliError"),
            Depolarize2 => write!(f, "Depolarize2"),
            Depolarize => write!(f, "Depolarize"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Message)]
pub enum CircuitMessage {
    None,                                           // Truncate (no operand)
    Qubit(usize),                                   // X, Y, Z, ...
    QubitAndFloat(usize, f64),                      // RX, depolarize, ...
    QubitAndTwoFloats(usize, f64, f64),             // R
    TwoQubit(usize, usize),                         // CX, CZ
    TwoQubitAndFloat(usize, usize, f64),            // RXX, ...
    QubitU3(usize, f64, f64, f64),                  // U3
    QubitAndFloatArr3(usize, [f64; 3]),             // PauliError
    TwoQubitAndFloatArr3(usize, usize, [f64; 3]),   // Correlated loss
    TwoQubitAndFloatArr15(usize, usize, [f64; 15]), // TwoQubitPauliError

    PauliPatternStr(String), // Trace (resolved Pauli-pattern source)

    // batched instructions
    QubitBatch(SmallVec<[usize; 8]>),              // X, Y, Z, ...
    QubitBatchAndFloat(SmallVec<[usize; 8]>, f64), // RX, depolarize, ...
    TwoQubitBatch(SmallVec<[(usize, usize); 8]>),  // CX, CZ
    TwoQubitBatchAndFloat(SmallVec<[(usize, usize); 8]>, f64), // RXX, ...
    QubitBatchU3(SmallVec<[usize; 8]>, f64, f64, f64), // U3
    QubitBatchAndFloatArr3(SmallVec<[usize; 8]>, [f64; 3]), // PauliError
    TwoQubitBatchAndFloatArr3(SmallVec<[(usize, usize); 8]>, [f64; 3]), // Correlated loss
    TwoQubitBatchAndFloatArr15(SmallVec<[(usize, usize); 8]>, [f64; 15]), // TwoQubitPauliError
}

#[derive(Debug, Clone)]
pub struct CircuitEffect {
    pub inst: CircuitInstruction,
    pub msg: CircuitMessage,
}

#[cfg(test)]
mod tests {
    use super::syntax::Instruction::*;
    use super::*;

    use chumsky::Parser as _;
    use vihaco::Parse;
    fn parse(src: &str) -> CircuitSurfaceInstruction {
        let src = format!("circuit.{src}");
        CircuitSurfaceInstruction::parser()
            .parse(src.as_str())
            .into_result()
            .unwrap_or_else(|e| panic!("parse of `{src}` failed: {e:?}"))
    }

    // ─── Parse: tokens are the lowercased variant name ────────────────────

    #[test]
    fn parses_simple_lowercase_tokens() {
        assert_eq!(parse("x"), X);
        assert_eq!(parse("y"), Y);
        assert_eq!(parse("z"), Z);
        assert_eq!(parse("h"), H);
        assert_eq!(parse("cnot"), CNOT);
        assert_eq!(parse("cz"), CZ);
        assert_eq!(parse("u3"), U3);
        assert_eq!(parse("measure"), Measure);
        assert_eq!(parse("reset"), Reset);
        assert_eq!(parse("r"), R);
        assert_eq!(parse("rx"), RX);
        assert_eq!(parse("ry"), RY);
        assert_eq!(parse("rz"), RZ);
        assert_eq!(parse("rxx"), RXX);
        assert_eq!(parse("depolarize2"), Depolarize2);
        assert_eq!(parse("depolarize"), Depolarize);
    }

    // ─── Parse: prefix-sensitive disambiguation ───────────────────────────
    //
    // These pairs share a prefix, so the declaration order in the enum is
    // load-bearing: the longer token must win. These tests pin that contract.

    #[test]
    fn parses_t_family_without_prefix_collision() {
        // `t` is a prefix of `tadj`, `trace`, `truncate`, and `twoqubitpaulierror`.
        assert_eq!(parse("t"), T);
        assert_eq!(parse("tadj"), TAdj);
        assert_eq!(parse("trace"), Trace);
        assert_eq!(parse("truncate"), Truncate);
        assert_eq!(parse("twoqubitpaulierror"), TwoQubitPauliError);
    }

    #[test]
    fn parses_s_family_without_prefix_collision() {
        // `s` is a prefix of `s_adj`, `sqrt_x`, `sqrt_y`, etc.
        assert_eq!(parse("s"), S);
        assert_eq!(parse("sadj"), SAdj);
        assert_eq!(parse("sqrtx"), SqrtX);
        assert_eq!(parse("sqrtxadj"), SqrtXAdj);
        assert_eq!(parse("sqrty"), SqrtY);
        assert_eq!(parse("sqrtyadj"), SqrtYAdj);
    }

    #[test]
    fn rejects_unknown_token() {
        assert!(CircuitSurfaceInstruction::parser()
            .parse("circuit.nope")
            .has_errors());
    }

    #[test]
    fn rejects_pascal_case_token() {
        // The parse token is lowercase; the Display form must not parse back.
        assert!(CircuitSurfaceInstruction::parser()
            .parse("circuit.CNOT")
            .has_errors());
    }

}
