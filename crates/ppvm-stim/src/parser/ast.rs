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

