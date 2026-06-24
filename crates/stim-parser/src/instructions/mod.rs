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
    /// arg count is [`ArgCount::Deferred`]). The lowering pass then promotes
    /// `I_ERROR[loss]` / `I_ERROR[correlated_loss]` to typed loss instructions
    /// and rejects every other tag combination (including untagged) by emitting
    /// an `"invalid-tag"` diagnostic.
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
    /// Targets must come in groups of four. Reserved — no instruction in the
    /// current `TABLE` uses it, but the validator handles it for completeness.
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
    pub canonical: &'static str,
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

const fn gate(
    name: GateName,
    args: ArgCount,
    targets: TargetArity,
    canonical: &'static str,
) -> TableEntry {
    TableEntry {
        kind: EntryKind::Gate(name),
        args,
        targets,
        canonical,
    }
}

const fn noise(
    name: NoiseName,
    args: ArgCount,
    targets: TargetArity,
    canonical: &'static str,
) -> TableEntry {
    TableEntry {
        kind: EntryKind::Noise(name),
        args,
        targets,
        canonical,
    }
}

const fn measure(name: MeasureName, canonical: &'static str) -> TableEntry {
    TableEntry {
        kind: EntryKind::Measure(name),
        args: ArgCount::Optional(1),
        targets: TargetArity::AtLeastOne,
        canonical,
    }
}

const fn measure_pairs(name: MeasureName, canonical: &'static str) -> TableEntry {
    TableEntry {
        kind: EntryKind::Measure(name),
        args: ArgCount::Optional(1),
        targets: TargetArity::Pairs,
        canonical,
    }
}

const fn annotation(kind: AnnotationKind, canonical: &'static str) -> TableEntry {
    TableEntry {
        kind: EntryKind::Annotation(kind),
        args: ArgCount::Any,
        targets: TargetArity::Any,
        canonical,
    }
}

use ArgCount::{Exact, None as NoArgs};
use GateName as G;
use MeasureName as Me;
use NoiseName as N;
use TargetArity::{AtLeastOne, Pairs};

/// Master instruction table. Entries are sorted by family for readability;
/// lookup is linear and runs at the start of every instruction parse —
/// performance is dwarfed by the rest of the pipeline.
const TABLE: &[(&str, TableEntry)] = &[
    // --- Gates: reset ---
    ("R", gate(G::Reset, NoArgs, AtLeastOne, "R")),
    ("RZ", gate(G::ResetZ, NoArgs, AtLeastOne, "RZ")),
    ("RX", gate(G::ResetX, NoArgs, AtLeastOne, "RX")),
    ("RY", gate(G::ResetY, NoArgs, AtLeastOne, "RY")),
    // --- Gates: single-qubit Clifford / paulis ---
    ("X", gate(G::X, NoArgs, AtLeastOne, "X")),
    ("Y", gate(G::Y, NoArgs, AtLeastOne, "Y")),
    ("Z", gate(G::Z, NoArgs, AtLeastOne, "Z")),
    ("H", gate(G::H, NoArgs, AtLeastOne, "H")),
    ("H_XZ", gate(G::HXZ, NoArgs, AtLeastOne, "H_XZ")),
    ("S", gate(G::S, NoArgs, AtLeastOne, "S")),
    ("S_DAG", gate(G::SDag, NoArgs, AtLeastOne, "S_DAG")),
    ("SQRT_Z", gate(G::SqrtZ, NoArgs, AtLeastOne, "SQRT_Z")),
    (
        "SQRT_Z_DAG",
        gate(G::SqrtZDag, NoArgs, AtLeastOne, "SQRT_Z_DAG"),
    ),
    ("SQRT_X", gate(G::SqrtX, NoArgs, AtLeastOne, "SQRT_X")),
    (
        "SQRT_X_DAG",
        gate(G::SqrtXDag, NoArgs, AtLeastOne, "SQRT_X_DAG"),
    ),
    ("SQRT_Y", gate(G::SqrtY, NoArgs, AtLeastOne, "SQRT_Y")),
    (
        "SQRT_Y_DAG",
        gate(G::SqrtYDag, NoArgs, AtLeastOne, "SQRT_Y_DAG"),
    ),
    ("T", gate(G::T, NoArgs, AtLeastOne, "T")),
    ("T_DAG", gate(G::TDag, NoArgs, AtLeastOne, "T_DAG")),
    ("I", gate(G::Identity, NoArgs, AtLeastOne, "I")),
    // --- Gates: two-qubit Clifford ---
    ("CX", gate(G::CX, NoArgs, Pairs, "CX")),
    ("ZCX", gate(G::ZCX, NoArgs, Pairs, "ZCX")),
    ("CNOT", gate(G::CNot, NoArgs, Pairs, "CNOT")),
    ("CY", gate(G::CY, NoArgs, Pairs, "CY")),
    ("ZCY", gate(G::ZCY, NoArgs, Pairs, "ZCY")),
    ("CZ", gate(G::CZ, NoArgs, Pairs, "CZ")),
    ("ZCZ", gate(G::ZCZ, NoArgs, Pairs, "ZCZ")),
    // --- Gates: phase-1-unsupported (parser accepts) ---
    ("SWAP", gate(G::Swap, NoArgs, Pairs, "SWAP")),
    ("ISWAP", gate(G::ISwap, NoArgs, Pairs, "ISWAP")),
    ("ISWAP_DAG", gate(G::ISwapDag, NoArgs, Pairs, "ISWAP_DAG")),
    ("SQRT_XX", gate(G::SqrtXX, NoArgs, Pairs, "SQRT_XX")),
    ("SQRT_YY", gate(G::SqrtYY, NoArgs, Pairs, "SQRT_YY")),
    ("SQRT_ZZ", gate(G::SqrtZZ, NoArgs, Pairs, "SQRT_ZZ")),
    ("CXSWAP", gate(G::CXSwap, NoArgs, Pairs, "CXSWAP")),
    ("SWAPCX", gate(G::SwapCX, NoArgs, Pairs, "SWAPCX")),
    ("XCX", gate(G::XCX, NoArgs, Pairs, "XCX")),
    ("XCY", gate(G::XCY, NoArgs, Pairs, "XCY")),
    ("XCZ", gate(G::XCZ, NoArgs, Pairs, "XCZ")),
    ("YCX", gate(G::YCX, NoArgs, Pairs, "YCX")),
    ("YCY", gate(G::YCY, NoArgs, Pairs, "YCY")),
    ("YCZ", gate(G::YCZ, NoArgs, Pairs, "YCZ")),
    ("C_XYZ", gate(G::CXYZ, NoArgs, AtLeastOne, "C_XYZ")),
    ("C_ZYX", gate(G::CZYX, NoArgs, AtLeastOne, "C_ZYX")),
    ("H_XY", gate(G::HXY, NoArgs, AtLeastOne, "H_XY")),
    ("H_YZ", gate(G::HYZ, NoArgs, AtLeastOne, "H_YZ")),
    // --- Noise ---
    (
        "DEPOLARIZE1",
        noise(N::Depolarize1, Exact(1), AtLeastOne, "DEPOLARIZE1"),
    ),
    (
        "DEPOLARIZE2",
        noise(N::Depolarize2, Exact(1), Pairs, "DEPOLARIZE2"),
    ),
    (
        "PAULI_CHANNEL_1",
        noise(N::PauliChannel1, Exact(3), AtLeastOne, "PAULI_CHANNEL_1"),
    ),
    (
        "PAULI_CHANNEL_2",
        noise(N::PauliChannel2, Exact(15), Pairs, "PAULI_CHANNEL_2"),
    ),
    ("X_ERROR", noise(N::XError, Exact(1), AtLeastOne, "X_ERROR")),
    ("Y_ERROR", noise(N::YError, Exact(1), AtLeastOne, "Y_ERROR")),
    ("Z_ERROR", noise(N::ZError, Exact(1), AtLeastOne, "Z_ERROR")),
    // I_ERROR's arg count varies by tag (`[loss]` => 1, `[correlated_loss]` => 1 or 3).
    // The extended-dialect layer enforces the tag-specific count.
    (
        "I_ERROR",
        noise(N::IError, ArgCount::Deferred, AtLeastOne, "I_ERROR"),
    ),
    (
        "HERALDED_ERASE",
        noise(N::HeraldedErase, Exact(1), AtLeastOne, "HERALDED_ERASE"),
    ),
    (
        "HERALDED_PAULI_CHANNEL_1",
        noise(
            N::HeraldedPauliChannel1,
            Exact(4),
            AtLeastOne,
            "HERALDED_PAULI_CHANNEL_1",
        ),
    ),
    (
        "CORRELATED_ERROR",
        noise(N::CorrelatedError, Exact(1), AtLeastOne, "CORRELATED_ERROR"),
    ),
    (
        "ELSE_CORRELATED_ERROR",
        noise(
            N::ElseCorrelatedError,
            Exact(1),
            AtLeastOne,
            "ELSE_CORRELATED_ERROR",
        ),
    ),
    // --- Measurements (all share Optional(1) args) ---
    ("M", measure(Me::M, "M")),
    ("MZ", measure(Me::MZ, "MZ")),
    ("MR", measure(Me::MR, "MR")),
    ("MX", measure(Me::MX, "MX")),
    ("MY", measure(Me::MY, "MY")),
    ("MRX", measure(Me::MRX, "MRX")),
    ("MRY", measure(Me::MRY, "MRY")),
    ("MXX", measure_pairs(Me::MXX, "MXX")),
    ("MYY", measure_pairs(Me::MYY, "MYY")),
    ("MZZ", measure_pairs(Me::MZZ, "MZZ")),
    ("MPP", measure(Me::MPP, "MPP")),
    // --- Annotations ---
    ("DETECTOR", annotation(AnnotationKind::Detector, "DETECTOR")),
    (
        "MPAD",
        TableEntry {
            kind: EntryKind::MPad,
            args: ArgCount::Optional(1),
            targets: AtLeastOne,
            canonical: "MPAD",
        },
    ),
    (
        "OBSERVABLE_INCLUDE",
        annotation(AnnotationKind::ObservableInclude, "OBSERVABLE_INCLUDE"),
    ),
    (
        "QUBIT_COORDS",
        annotation(AnnotationKind::QubitCoords, "QUBIT_COORDS"),
    ),
    (
        "SHIFT_COORDS",
        annotation(AnnotationKind::ShiftCoords, "SHIFT_COORDS"),
    ),
    (
        "TICK",
        TableEntry {
            kind: EntryKind::Annotation(AnnotationKind::Tick),
            args: ArgCount::None,
            targets: TargetArity::Any,
            canonical: "TICK",
        },
    ),
];

use std::fmt;

impl GateName {
    pub fn canonical_name(self) -> &'static str {
        canonical_name(EntryKind::Gate(self))
    }
}

impl fmt::Display for GateName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.canonical_name())
    }
}

impl NoiseName {
    pub fn canonical_name(self) -> &'static str {
        canonical_name(EntryKind::Noise(self))
    }
}

impl MeasureName {
    pub fn canonical_name(self) -> &'static str {
        canonical_name(EntryKind::Measure(self))
    }
}

impl AnnotationKind {
    pub fn canonical_name(self) -> &'static str {
        canonical_name(EntryKind::Annotation(self))
    }
}

/// Look up a Stim instruction name. `None` means unknown.
pub fn lookup(name: &str) -> Option<TableEntry> {
    TABLE.iter().find(|(n, _)| *n == name).map(|(_, e)| *e)
}

/// Canonical spelling for a decoded instruction kind, derived from the
/// same TABLE rows that drive parsing (single source of truth).
pub fn canonical_name(kind: EntryKind) -> &'static str {
    TABLE
        .iter()
        .find(|(_, e)| e.kind == kind)
        .map(|(_, e)| e.canonical)
        .expect("every EntryKind has a TABLE row (enforced by completeness test)")
}

#[cfg(test)]
mod table_tests {
    use super::*;

    #[test]
    fn every_table_key_is_unique() {
        use std::collections::HashSet;
        let mut seen = HashSet::new();
        for (key, _) in TABLE {
            assert!(seen.insert(*key), "duplicate key {key:?} in TABLE");
        }
    }

    #[test]
    fn canonical_round_trips_through_lookup() {
        for (_, entry) in TABLE {
            let via_canonical = lookup(entry.canonical)
                .unwrap_or_else(|| panic!("canonical {:?} not in TABLE", entry.canonical));
            assert_eq!(
                via_canonical.kind, entry.kind,
                "canonical mismatch for {:?}",
                entry.canonical
            );
        }
    }

    #[test]
    fn every_row_kind_resolves_to_canonical() {
        for (_, entry) in TABLE {
            let _ = canonical_name(entry.kind); // must not panic
        }
    }
}
