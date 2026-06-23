// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Single source of truth for Stim instruction identity and parse/validation
//! metadata: the name enums, arity rules, and the lookup table.

/// All Stim gate names. The parser accepts every variant; consumers may
/// reject the ones their backend doesn't support.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GateName {
    // Reset (treated as a gate so it parses with no args)
    Reset,
    ResetZ,
    // X/Y-basis resets (prepare |+> / |i>)
    ResetX,
    ResetY,
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
    // Non-Clifford T / T-dagger (also expressible as S[T] / S_DAG[T])
    T,
    TDag,
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
    /// Any arg count is accepted.
    Any,
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

/// Decoded instruction-table entry: family discriminant plus arity rules.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TableEntry {
    pub kind: EntryKind,
    pub args: ArgCount,
    pub targets: TargetArity,
}

/// Which AST family the instruction belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntryKind {
    Gate(GateName),
    Noise(NoiseName),
    Measure(MeasureName),
    Annotation(AnnotationKind),
    MPad,
}

#[expect(dead_code, reason = "used by TABLE in the next task")]
const fn gate(name: GateName, args: ArgCount, targets: TargetArity) -> TableEntry {
    TableEntry {
        kind: EntryKind::Gate(name),
        args,
        targets,
    }
}

#[expect(dead_code, reason = "used by TABLE in the next task")]
const fn noise(name: NoiseName, args: ArgCount, targets: TargetArity) -> TableEntry {
    TableEntry {
        kind: EntryKind::Noise(name),
        args,
        targets,
    }
}

#[expect(dead_code, reason = "used by TABLE in the next task")]
const fn measure(name: MeasureName) -> TableEntry {
    TableEntry {
        kind: EntryKind::Measure(name),
        args: ArgCount::Optional(1),
        targets: TargetArity::AtLeastOne,
    }
}

#[expect(dead_code, reason = "used by TABLE in the next task")]
const fn measure_pairs(name: MeasureName) -> TableEntry {
    TableEntry {
        kind: EntryKind::Measure(name),
        args: ArgCount::Optional(1),
        targets: TargetArity::Pairs,
    }
}

#[expect(dead_code, reason = "used by TABLE in the next task")]
const fn annotation(kind: AnnotationKind) -> TableEntry {
    TableEntry {
        kind: EntryKind::Annotation(kind),
        args: ArgCount::Any,
        targets: TargetArity::Any,
    }
}
