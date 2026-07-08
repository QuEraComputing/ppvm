// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use smallvec::SmallVec;
use vihaco::Instruction;
use vihaco::Message;
use vihaco_parser::Parse;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Instruction, Parse)]
pub enum CircuitInstruction {
    // NOTE: longer tokens need to go first
    TwoQubitPauliError, // needs to go before T
    Truncate,           // needs to go before T
    Trace,              // needs to go before T

    // Single-Qubit Clifford gates
    X,
    Y,
    Z,
    H,

    #[token = "sqrt_x_adj"]
    SqrtXAdj,

    #[token = "sqrt_x"]
    SqrtX,

    #[token = "sqrt_y_adj"]
    SqrtYAdj,

    #[token = "sqrt_y"]
    SqrtY,

    #[token = "s_adj"]
    SAdj,
    S,

    // Controlled gates
    CNOT,
    CZ,

    // T gate
    TAdj,
    T,

    // Two-qubit rotations
    RXX,
    RYY,
    RZZ,

    // Single-qubit rotations
    RX,
    RY,
    RZ,

    // U3
    U3,

    // Measurement & Reset
    Measure,
    Reset,

    // RXY
    R,

    // Loss
    Loss,
    CorrelatedLoss,

    // Noise
    PauliError,
    Depolarize2,
    Depolarize,
}

impl std::fmt::Display for CircuitInstruction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use CircuitInstruction::*;
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
    use super::CircuitInstruction::*;
    use super::*;

    use chumsky::Parser as _;
    use vihaco::instruction::{FromBytes, OpCode, WriteBytes};
    use vihaco_parser_core::Parse as _;

    /// Every variant, in declaration order. Anything iterating over the full
    /// instruction set (round-trips, opcode uniqueness) goes through this so a
    /// newly added variant is automatically covered.
    const ALL: &[CircuitInstruction] = &[
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
    ];

    fn parse(src: &str) -> CircuitInstruction {
        CircuitInstruction::parser()
            .parse(src)
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
        assert_eq!(parse("s_adj"), SAdj);
        assert_eq!(parse("sqrt_x"), SqrtX);
        assert_eq!(parse("sqrt_x_adj"), SqrtXAdj);
        assert_eq!(parse("sqrt_y"), SqrtY);
        assert_eq!(parse("sqrt_y_adj"), SqrtYAdj);
    }

    #[test]
    fn rejects_unknown_token() {
        assert!(CircuitInstruction::parser().parse("nope").has_errors());
    }

    #[test]
    fn rejects_pascal_case_token() {
        // The parse token is lowercase; the Display form must not parse back.
        assert!(CircuitInstruction::parser().parse("CNOT").has_errors());
    }

    // ─── Display: PascalCase variant names ────────────────────────────────

    #[test]
    fn display_uses_pascal_case_names() {
        assert_eq!(H.to_string(), "H");
        assert_eq!(CNOT.to_string(), "CNOT");
        assert_eq!(TwoQubitPauliError.to_string(), "TwoQubitPauliError");
        assert_eq!(Trace.to_string(), "Trace");
        assert_eq!(Truncate.to_string(), "Truncate");
        // Custom-token variants display their Rust name, not the parse token.
        assert_eq!(SqrtXAdj.to_string(), "SqrtXAdj");
        assert_eq!(SAdj.to_string(), "SAdj");
    }

    // ─── Instruction codec (derived OpCode / WriteBytes / FromBytes) ──────

    #[test]
    fn opcodes_are_unit_width() {
        // All variants are field-less, so each encodes to a single byte.
        assert_eq!(CircuitInstruction::width(), 1);
    }

    #[test]
    fn opcodes_are_unique() {
        let mut seen = std::collections::HashSet::new();
        for inst in ALL {
            assert!(
                seen.insert(inst.opcode()),
                "duplicate opcode {} for {inst:?}",
                inst.opcode()
            );
        }
        assert_eq!(seen.len(), ALL.len());
    }

    #[test]
    fn opcodes_match_declaration_order() {
        // Opcodes default to the variant index, which is the on-disk contract
        // for bytecode. Reordering variants silently breaks old bytecode, so
        // pin the assignment here.
        for (index, inst) in ALL.iter().enumerate() {
            assert_eq!(inst.opcode() as usize, index, "{inst:?}");
        }
    }

    #[test]
    fn write_then_read_round_trips_every_variant() {
        for inst in ALL {
            let mut buf = Vec::new();
            inst.write_bytes(&mut buf).unwrap();
            assert_eq!(buf, [inst.opcode()], "{inst:?} should encode to one byte");

            let mut cursor = std::io::Cursor::new(buf);
            let back = CircuitInstruction::from_bytes(&mut cursor).unwrap();
            assert_eq!(back, *inst);
        }
    }

    #[test]
    fn from_bytes_rejects_unknown_opcode() {
        let mut cursor = std::io::Cursor::new([0xFFu8]);
        let err = CircuitInstruction::from_bytes(&mut cursor).unwrap_err();
        assert!(err.to_string().contains("invalid opcode"), "err: {err}");
    }
}
