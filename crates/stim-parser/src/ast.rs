//! Pure-Stim AST. Tags are preserved verbatim; the parser does not
//! resolve the Stim dialect — that is the consumer's responsibility.

use std::fmt;
use std::sync::Arc;

use chumsky::error::Rich;

use super::LineMap;

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
    MPad {
        tags: Vec<Tag>,
        prob: Option<f64>,
        bits: Vec<usize>,
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

/// All Stim gate names. The parser accepts every variant; consumers may
/// reject the ones their backend doesn't support.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GateName {
    // Reset (treated as a gate so it parses with no args)
    Reset,
    ResetZ,
    // Single-qubit Cliffords
    X,
    Y,
    Z,
    H,
    HXZ,
    S,
    SDag,
    SqrtZ,
    SqrtZDag,
    SqrtX,
    SqrtXDag,
    SqrtY,
    SqrtYDag,
    // Identity (carries dialect tags like S[T] is on S, but I[R_X(...)] is on Identity)
    Identity,
    // Two-qubit Cliffords
    CX,
    ZCX,
    CNot,
    CY,
    ZCY,
    CZ,
    ZCZ,
    // Phase-1-unsupported (parser only)
    Swap,
    ISwap,
    ISwapDag,
    SqrtXX,
    SqrtYY,
    SqrtZZ,
    CXSwap,
    SwapCX,
    XCX,
    XCY,
    XCZ,
    YCX,
    YCY,
    YCZ,
    CXYZ,
    CZYX,
    HXY,
    HYZ,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoiseName {
    // Supported
    Depolarize1,
    Depolarize2,
    PauliChannel1,
    PauliChannel2,
    XError,
    YError,
    ZError,
    /// `I_ERROR` — accepted at the raw-parse layer regardless of tag (its
    /// arg count is [`ArgCount::Deferred`]). The extended-dialect interpreter
    /// then promotes `I_ERROR[loss]` / `I_ERROR[correlated_loss]` to typed
    /// loss instructions and rejects every other tag combination (including
    /// untagged) as [`ExtendedParseError::InvalidTag`].
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
    M,
    MZ,
    MR,
    // Unsupported
    MX,
    MY,
    MRX,
    MRY,
    MXX,
    MYY,
    MZZ,
    MPP,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotationKind {
    Detector,
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
    /// Either no args or exactly `n` args. Used by Stim measurement
    /// instructions, where the optional single arg is the readout-flip
    /// probability.
    Optional(usize),
    /// Skip parse-time arg validation. The downstream layer (extended
    /// dialect) enforces an instruction-specific rule based on tags.
    Deferred,
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

/// Carries a chumsky 0.12 `Rich<char>` error plus a shared `LineMap`
/// so that `Display` formats `line:col` consistently with the typed
/// validation variants of [`ParseError`].
#[derive(Debug, Clone)]
pub struct SyntaxError {
    pub(crate) rich: Rich<'static, char>,
    pub(crate) line_map: Arc<LineMap>,
}

impl SyntaxError {
    /// Synthesise a `SyntaxError` from a (line, col) position. Used by
    /// the validator when it needs to emit a syntax error for an issue
    /// the grammar could not catch (e.g. an annotation-tolerated target
    /// that fails `usize` parsing for a non-annotation instruction).
    pub(crate) fn synth(
        line: usize,
        col: usize,
        message: impl Into<String>,
        line_map: Arc<LineMap>,
    ) -> Self {
        let line_idx = line.saturating_sub(1);
        let line_start = line_map.starts_at(line_idx).unwrap_or(0);
        let byte = line_start + col.saturating_sub(1);
        let span = chumsky::span::SimpleSpan::from(byte..byte);
        let rich = Rich::<char>::custom(span, message.into());
        SyntaxError { rich, line_map }
    }

    /// Construct from a chumsky `Rich` (with any borrow lifetime) and
    /// a shared `LineMap`. `into_owned` widens the lifetime to `'static`.
    pub(crate) fn from_rich<'src>(rich: Rich<'src, char>, line_map: Arc<LineMap>) -> Self {
        SyntaxError {
            rich: rich.into_owned(),
            line_map,
        }
    }
}

impl fmt::Display for SyntaxError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let span = self.rich.span();
        let (line, col) = self.line_map.line_col(span.start);
        write!(f, "syntax error at line {line}, col {col}: {}", self.rich)
    }
}

#[derive(Debug, thiserror::Error)]
pub enum ParseError {
    #[error("{0}")]
    Syntax(SyntaxError),

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
            AnnotationKind::ObservableInclude => "OBSERVABLE_INCLUDE",
            AnnotationKind::QubitCoords => "QUBIT_COORDS",
            AnnotationKind::ShiftCoords => "SHIFT_COORDS",
            AnnotationKind::Tick => "TICK",
        }
    }
}
