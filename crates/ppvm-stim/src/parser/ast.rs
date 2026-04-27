//! Pure-Stim AST. Tags are preserved verbatim; dialect resolution
//! happens in `crate::normalize`.

use std::fmt;

/// Top-level program: a flat list of instructions (REPEAT bodies are nested).
#[derive(Debug, Clone, PartialEq)]
pub struct Program {
    pub instructions: Vec<RawInstruction>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum RawInstruction {
    Gate {
        name: GateName,
        tags: Vec<Tag>,
        args: Vec<f64>,
        targets: Vec<usize>,
        line: usize,
    },
    Noise {
        name: NoiseName,
        tags: Vec<Tag>,
        args: Vec<f64>,
        targets: Vec<usize>,
        line: usize,
    },
    Measure {
        name: MeasureName,
        tags: Vec<Tag>,
        args: Vec<f64>,
        targets: Vec<usize>,
        line: usize,
    },
    Annotation {
        kind: AnnotationKind,
        args: Vec<f64>,
        targets: Vec<usize>,
        line: usize,
    },
    Repeat {
        count: u64,
        body: Vec<RawInstruction>,
        line: usize,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct Tag {
    pub name: String,
    pub params: Vec<TagParam>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum TagParam {
    Positional(f64),
    Named { key: String, value: f64 },
}

/// Stim gate names — phase-1-supported and Stim-valid-but-unsupported.
/// Parser accepts every variant; normalizer rejects the unsupported ones.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GateName {
    // Reset (treated as a gate so it parses with no args)
    Reset, ResetZ,
    // Single-qubit Cliffords
    X, Y, Z, H, HXZ,
    S, SDag,
    SqrtZ, SqrtZDag,
    SqrtX, SqrtXDag,
    SqrtY, SqrtYDag,
    // Identity (carries dialect tags like S[T] is on S, but I[R_X(...)] is on Identity)
    Identity,
    // Two-qubit Cliffords
    CX, ZCX, CNot,
    CY, ZCY,
    CZ, ZCZ,
    // Phase-1-unsupported (parser only)
    Swap, ISwap, ISwapDag,
    SqrtXX, SqrtYY, SqrtZZ,
    CXSwap, SwapCX,
    XCX, XCY, XCZ, YCX, YCY, YCZ,
    CXYZ, CZYX,
    HXY, HYZ,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoiseName {
    // Supported
    Depolarize1, Depolarize2,
    PauliChannel1, PauliChannel2,
    XError, YError, ZError,
    /// `I_ERROR` — supported when tagged `[loss]` / `[correlated_loss]`,
    /// rejected otherwise (matches today's `stim.rs` behavior of dropping
    /// untagged `I_ERROR`).
    IError,
    // Unsupported
    HeraldedErase,
    HeraldedPauliChannel1,
    CorrelatedError,
    ElseCorrelatedError,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeasureName {
    // Supported
    M, MZ, MR,
    // Unsupported
    MX, MY, MRX, MRY,
    MXX, MYY, MZZ, MPP,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotationKind {
    Detector,
    /// `MPAD` is treated as an annotation, not a measurement (matches today's `stim.rs`).
    MPad,
    ObservableInclude,
    QubitCoords,
    ShiftCoords,
    Tick,
}

/// Required argument-count rule for a Stim instruction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArgCount {
    /// No args allowed. `(…)` parens must be absent.
    None,
    /// Exactly `n` args.
    Exact(usize),
    /// Args may be present or absent; if present, the count must equal `n`.
    /// Used for noise channels where the parens are mandatory but we don't
    /// special-case "no parens" vs "parens with wrong count" separately.
    /// (Not currently used; placeholder for future flexibility.)
    Optional(usize),
}

/// Required target-count rule for a Stim instruction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TargetArity {
    /// Any non-negative number of targets.
    Any,
    /// Targets must come in pairs (i.e. `len % 2 == 0`).
    Pairs,
    /// Targets must come in groups of four.
    Quadruples,
    /// At least one target required.
    AtLeastOne,
}

#[derive(Debug, thiserror::Error, PartialEq)]
#[non_exhaustive]
pub enum ParseError {
    #[error("syntax error at line {line}, col {col}: {message}")]
    Syntax { line: usize, col: usize, message: String },

    #[error("unknown instruction '{name}' at line {line}")]
    UnknownInstruction { name: String, line: usize },

    #[error("'{name}' at line {line} expected {expected} args, got {found}")]
    ArgCount {
        name: String,
        expected: usize,
        found: usize,
        line: usize,
    },

    #[error("'{name}' at line {line} expected target count divisible by {divisor}, got {found}")]
    TargetCount {
        name: String,
        divisor: usize,
        found: usize,
        line: usize,
    },
}

impl fmt::Display for GateName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Used in error messages — matches the canonical Stim spelling.
        f.write_str(self.canonical_name())
    }
}

impl GateName {
    pub fn canonical_name(&self) -> &'static str {
        match self {
            GateName::Reset => "R",
            GateName::ResetZ => "RZ",
            GateName::X => "X",
            GateName::Y => "Y",
            GateName::Z => "Z",
            GateName::H => "H",
            GateName::HXZ => "H_XZ",
            GateName::S => "S",
            GateName::SDag => "S_DAG",
            GateName::SqrtZ => "SQRT_Z",
            GateName::SqrtZDag => "SQRT_Z_DAG",
            GateName::SqrtX => "SQRT_X",
            GateName::SqrtXDag => "SQRT_X_DAG",
            GateName::SqrtY => "SQRT_Y",
            GateName::SqrtYDag => "SQRT_Y_DAG",
            GateName::Identity => "I",
            GateName::CX => "CX",
            GateName::ZCX => "ZCX",
            GateName::CNot => "CNOT",
            GateName::CY => "CY",
            GateName::ZCY => "ZCY",
            GateName::CZ => "CZ",
            GateName::ZCZ => "ZCZ",
            GateName::Swap => "SWAP",
            GateName::ISwap => "ISWAP",
            GateName::ISwapDag => "ISWAP_DAG",
            GateName::SqrtXX => "SQRT_XX",
            GateName::SqrtYY => "SQRT_YY",
            GateName::SqrtZZ => "SQRT_ZZ",
            GateName::CXSwap => "CXSWAP",
            GateName::SwapCX => "SWAPCX",
            GateName::XCX => "XCX",
            GateName::XCY => "XCY",
            GateName::XCZ => "XCZ",
            GateName::YCX => "YCX",
            GateName::YCY => "YCY",
            GateName::YCZ => "YCZ",
            GateName::CXYZ => "C_XYZ",
            GateName::CZYX => "C_ZYX",
            GateName::HXY => "H_XY",
            GateName::HYZ => "H_YZ",
        }
    }
}

impl NoiseName {
    pub fn canonical_name(&self) -> &'static str {
        match self {
            NoiseName::Depolarize1 => "DEPOLARIZE1",
            NoiseName::Depolarize2 => "DEPOLARIZE2",
            NoiseName::PauliChannel1 => "PAULI_CHANNEL_1",
            NoiseName::PauliChannel2 => "PAULI_CHANNEL_2",
            NoiseName::XError => "X_ERROR",
            NoiseName::YError => "Y_ERROR",
            NoiseName::ZError => "Z_ERROR",
            NoiseName::IError => "I_ERROR",
            NoiseName::HeraldedErase => "HERALDED_ERASE",
            NoiseName::HeraldedPauliChannel1 => "HERALDED_PAULI_CHANNEL_1",
            NoiseName::CorrelatedError => "CORRELATED_ERROR",
            NoiseName::ElseCorrelatedError => "ELSE_CORRELATED_ERROR",
        }
    }
}

impl MeasureName {
    pub fn canonical_name(&self) -> &'static str {
        match self {
            MeasureName::M => "M",
            MeasureName::MZ => "MZ",
            MeasureName::MR => "MR",
            MeasureName::MX => "MX",
            MeasureName::MY => "MY",
            MeasureName::MRX => "MRX",
            MeasureName::MRY => "MRY",
            MeasureName::MXX => "MXX",
            MeasureName::MYY => "MYY",
            MeasureName::MZZ => "MZZ",
            MeasureName::MPP => "MPP",
        }
    }
}

impl AnnotationKind {
    pub fn canonical_name(&self) -> &'static str {
        match self {
            AnnotationKind::Detector => "DETECTOR",
            AnnotationKind::MPad => "MPAD",
            AnnotationKind::ObservableInclude => "OBSERVABLE_INCLUDE",
            AnnotationKind::QubitCoords => "QUBIT_COORDS",
            AnnotationKind::ShiftCoords => "SHIFT_COORDS",
            AnnotationKind::Tick => "TICK",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TableEntry {
    Gate {
        name: GateName,
        args: ArgCount,
        targets: TargetArity,
    },
    Noise {
        name: NoiseName,
        args: ArgCount,
        targets: TargetArity,
    },
    Measure {
        name: MeasureName,
        args: ArgCount,
        targets: TargetArity,
    },
    Annotation {
        kind: AnnotationKind,
        args: ArgCount,
        targets: TargetArity,
    },
}

/// Master instruction table. The parser looks every instruction name up here.
/// Entries are sorted by family for readability; lookup is linear and runs at
/// the start of every instruction parse — performance is dwarfed by the rest
/// of the pipeline.
const TABLE: &[(&str, TableEntry)] = &[
    // --- Gates: reset ---
    ("R",  TableEntry::Gate { name: GateName::Reset,  args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("RZ", TableEntry::Gate { name: GateName::ResetZ, args: ArgCount::None, targets: TargetArity::AtLeastOne }),

    // --- Gates: single-qubit Clifford / paulis ---
    ("X",          TableEntry::Gate { name: GateName::X,         args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("Y",          TableEntry::Gate { name: GateName::Y,         args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("Z",          TableEntry::Gate { name: GateName::Z,         args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("H",          TableEntry::Gate { name: GateName::H,         args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("H_XZ",       TableEntry::Gate { name: GateName::HXZ,       args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("S",          TableEntry::Gate { name: GateName::S,         args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("S_DAG",      TableEntry::Gate { name: GateName::SDag,      args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("SQRT_Z",     TableEntry::Gate { name: GateName::SqrtZ,     args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("SQRT_Z_DAG", TableEntry::Gate { name: GateName::SqrtZDag,  args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("SQRT_X",     TableEntry::Gate { name: GateName::SqrtX,     args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("SQRT_X_DAG", TableEntry::Gate { name: GateName::SqrtXDag,  args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("SQRT_Y",     TableEntry::Gate { name: GateName::SqrtY,     args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("SQRT_Y_DAG", TableEntry::Gate { name: GateName::SqrtYDag,  args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("I",          TableEntry::Gate { name: GateName::Identity,  args: ArgCount::None, targets: TargetArity::AtLeastOne }),

    // --- Gates: two-qubit Clifford ---
    ("CX",   TableEntry::Gate { name: GateName::CX,   args: ArgCount::None, targets: TargetArity::Pairs }),
    ("ZCX",  TableEntry::Gate { name: GateName::ZCX,  args: ArgCount::None, targets: TargetArity::Pairs }),
    ("CNOT", TableEntry::Gate { name: GateName::CNot, args: ArgCount::None, targets: TargetArity::Pairs }),
    ("CY",   TableEntry::Gate { name: GateName::CY,   args: ArgCount::None, targets: TargetArity::Pairs }),
    ("ZCY",  TableEntry::Gate { name: GateName::ZCY,  args: ArgCount::None, targets: TargetArity::Pairs }),
    ("CZ",   TableEntry::Gate { name: GateName::CZ,   args: ArgCount::None, targets: TargetArity::Pairs }),
    ("ZCZ",  TableEntry::Gate { name: GateName::ZCZ,  args: ArgCount::None, targets: TargetArity::Pairs }),

    // --- Gates: phase-1-unsupported (parser accepts) ---
    ("SWAP",      TableEntry::Gate { name: GateName::Swap,     args: ArgCount::None, targets: TargetArity::Pairs }),
    ("ISWAP",     TableEntry::Gate { name: GateName::ISwap,    args: ArgCount::None, targets: TargetArity::Pairs }),
    ("ISWAP_DAG", TableEntry::Gate { name: GateName::ISwapDag, args: ArgCount::None, targets: TargetArity::Pairs }),
    ("SQRT_XX",   TableEntry::Gate { name: GateName::SqrtXX,   args: ArgCount::None, targets: TargetArity::Pairs }),
    ("SQRT_YY",   TableEntry::Gate { name: GateName::SqrtYY,   args: ArgCount::None, targets: TargetArity::Pairs }),
    ("SQRT_ZZ",   TableEntry::Gate { name: GateName::SqrtZZ,   args: ArgCount::None, targets: TargetArity::Pairs }),
    ("CXSWAP",    TableEntry::Gate { name: GateName::CXSwap,   args: ArgCount::None, targets: TargetArity::Pairs }),
    ("SWAPCX",    TableEntry::Gate { name: GateName::SwapCX,   args: ArgCount::None, targets: TargetArity::Pairs }),
    ("XCX",       TableEntry::Gate { name: GateName::XCX,      args: ArgCount::None, targets: TargetArity::Pairs }),
    ("XCY",       TableEntry::Gate { name: GateName::XCY,      args: ArgCount::None, targets: TargetArity::Pairs }),
    ("XCZ",       TableEntry::Gate { name: GateName::XCZ,      args: ArgCount::None, targets: TargetArity::Pairs }),
    ("YCX",       TableEntry::Gate { name: GateName::YCX,      args: ArgCount::None, targets: TargetArity::Pairs }),
    ("YCY",       TableEntry::Gate { name: GateName::YCY,      args: ArgCount::None, targets: TargetArity::Pairs }),
    ("YCZ",       TableEntry::Gate { name: GateName::YCZ,      args: ArgCount::None, targets: TargetArity::Pairs }),
    ("C_XYZ",     TableEntry::Gate { name: GateName::CXYZ,     args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("C_ZYX",     TableEntry::Gate { name: GateName::CZYX,     args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("H_XY",      TableEntry::Gate { name: GateName::HXY,      args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("H_YZ",      TableEntry::Gate { name: GateName::HYZ,      args: ArgCount::None, targets: TargetArity::AtLeastOne }),

    // --- Noise ---
    ("DEPOLARIZE1",     TableEntry::Noise { name: NoiseName::Depolarize1,    args: ArgCount::Exact(1),  targets: TargetArity::AtLeastOne }),
    ("DEPOLARIZE2",     TableEntry::Noise { name: NoiseName::Depolarize2,    args: ArgCount::Exact(1),  targets: TargetArity::Pairs }),
    ("PAULI_CHANNEL_1", TableEntry::Noise { name: NoiseName::PauliChannel1,  args: ArgCount::Exact(3),  targets: TargetArity::AtLeastOne }),
    ("PAULI_CHANNEL_2", TableEntry::Noise { name: NoiseName::PauliChannel2,  args: ArgCount::Exact(15), targets: TargetArity::Pairs }),
    ("X_ERROR",         TableEntry::Noise { name: NoiseName::XError,         args: ArgCount::Exact(1),  targets: TargetArity::AtLeastOne }),
    ("Y_ERROR",         TableEntry::Noise { name: NoiseName::YError,         args: ArgCount::Exact(1),  targets: TargetArity::AtLeastOne }),
    ("Z_ERROR",         TableEntry::Noise { name: NoiseName::ZError,         args: ArgCount::Exact(1),  targets: TargetArity::AtLeastOne }),
    // I_ERROR's arg count varies by tag (`[loss]` => 1, `[correlated_loss]` => 1 or 3).
    // Validate "args present" at parse time (any count); the normalizer enforces the
    // tag-specific count. We model this with ArgCount::None to skip parse-time arg
    // validation here, and let normalize emit InvalidTag for malformed tags.
    ("I_ERROR",                 TableEntry::Noise { name: NoiseName::IError,                 args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("HERALDED_ERASE",          TableEntry::Noise { name: NoiseName::HeraldedErase,          args: ArgCount::Exact(1), targets: TargetArity::AtLeastOne }),
    ("HERALDED_PAULI_CHANNEL_1",TableEntry::Noise { name: NoiseName::HeraldedPauliChannel1,  args: ArgCount::Exact(4), targets: TargetArity::AtLeastOne }),
    ("CORRELATED_ERROR",        TableEntry::Noise { name: NoiseName::CorrelatedError,        args: ArgCount::Exact(1), targets: TargetArity::AtLeastOne }),
    ("ELSE_CORRELATED_ERROR",   TableEntry::Noise { name: NoiseName::ElseCorrelatedError,    args: ArgCount::Exact(1), targets: TargetArity::AtLeastOne }),

    // --- Measurements ---
    ("M",   TableEntry::Measure { name: MeasureName::M,   args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("MZ",  TableEntry::Measure { name: MeasureName::MZ,  args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("MR",  TableEntry::Measure { name: MeasureName::MR,  args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("MX",  TableEntry::Measure { name: MeasureName::MX,  args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("MY",  TableEntry::Measure { name: MeasureName::MY,  args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("MRX", TableEntry::Measure { name: MeasureName::MRX, args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("MRY", TableEntry::Measure { name: MeasureName::MRY, args: ArgCount::None, targets: TargetArity::AtLeastOne }),
    ("MXX", TableEntry::Measure { name: MeasureName::MXX, args: ArgCount::None, targets: TargetArity::Pairs }),
    ("MYY", TableEntry::Measure { name: MeasureName::MYY, args: ArgCount::None, targets: TargetArity::Pairs }),
    ("MZZ", TableEntry::Measure { name: MeasureName::MZZ, args: ArgCount::None, targets: TargetArity::Pairs }),
    ("MPP", TableEntry::Measure { name: MeasureName::MPP, args: ArgCount::None, targets: TargetArity::AtLeastOne }),

    // --- Annotations ---
    ("DETECTOR",           TableEntry::Annotation { kind: AnnotationKind::Detector,          args: ArgCount::None, targets: TargetArity::Any }),
    ("MPAD",               TableEntry::Annotation { kind: AnnotationKind::MPad,              args: ArgCount::None, targets: TargetArity::Any }),
    ("OBSERVABLE_INCLUDE", TableEntry::Annotation { kind: AnnotationKind::ObservableInclude, args: ArgCount::None, targets: TargetArity::Any }),
    ("QUBIT_COORDS",       TableEntry::Annotation { kind: AnnotationKind::QubitCoords,       args: ArgCount::None, targets: TargetArity::Any }),
    ("SHIFT_COORDS",       TableEntry::Annotation { kind: AnnotationKind::ShiftCoords,       args: ArgCount::None, targets: TargetArity::Any }),
    ("TICK",               TableEntry::Annotation { kind: AnnotationKind::Tick,              args: ArgCount::None, targets: TargetArity::Any }),
];

/// Look up a Stim instruction name. `None` means unknown.
pub fn lookup(name: &str) -> Option<TableEntry> {
    TABLE.iter().find(|(n, _)| *n == name).map(|(_, e)| *e)
}

#[cfg(test)]
mod table_tests {
    use super::*;

    #[test]
    fn lookup_h_returns_gate_h_with_arity_any_no_args() {
        let entry = lookup("H").expect("H must be known");
        assert_eq!(
            entry,
            TableEntry::Gate {
                name: GateName::H,
                args: ArgCount::None,
                targets: TargetArity::AtLeastOne,
            }
        );
    }

    #[test]
    fn lookup_cx_requires_pairs() {
        let entry = lookup("CX").expect("CX must be known");
        assert!(matches!(
            entry,
            TableEntry::Gate { name: GateName::CX, args: ArgCount::None, targets: TargetArity::Pairs }
        ));
    }

    #[test]
    fn lookup_depolarize1_requires_one_arg_any_targets() {
        let entry = lookup("DEPOLARIZE1").expect("DEPOLARIZE1 must be known");
        assert!(matches!(
            entry,
            TableEntry::Noise { name: NoiseName::Depolarize1, args: ArgCount::Exact(1), targets: TargetArity::AtLeastOne }
        ));
    }

    #[test]
    fn lookup_pauli_channel_2_requires_15_args_pair_targets() {
        let entry = lookup("PAULI_CHANNEL_2").expect("PAULI_CHANNEL_2 must be known");
        assert!(matches!(
            entry,
            TableEntry::Noise { name: NoiseName::PauliChannel2, args: ArgCount::Exact(15), targets: TargetArity::Pairs }
        ));
    }

    #[test]
    fn lookup_m_returns_measure() {
        let entry = lookup("M").expect("M must be known");
        assert!(matches!(
            entry,
            TableEntry::Measure { name: MeasureName::M, .. }
        ));
    }

    #[test]
    fn lookup_detector_returns_annotation() {
        let entry = lookup("DETECTOR").expect("DETECTOR must be known");
        assert!(matches!(
            entry,
            TableEntry::Annotation { kind: AnnotationKind::Detector, .. }
        ));
    }

    #[test]
    fn lookup_unknown_returns_none() {
        assert!(lookup("FROBNICATE").is_none());
    }

    #[test]
    fn aliases_map_to_distinct_variants() {
        // CNOT, CX, ZCX are all spelled differently and produce distinct GateName
        // variants — normalize.rs will treat them as the same gate.
        assert!(matches!(lookup("CNOT").unwrap(),
            TableEntry::Gate { name: GateName::CNot, .. }));
        assert!(matches!(lookup("CX").unwrap(),
            TableEntry::Gate { name: GateName::CX, .. }));
        assert!(matches!(lookup("ZCX").unwrap(),
            TableEntry::Gate { name: GateName::ZCX, .. }));
    }
}
