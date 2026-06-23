// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Instruction lookup table for the Stim parser.
//!
//! [`lookup`] maps a raw instruction name (e.g. `"H"`) to a [`TableEntry`]
//! that describes its argument and target-arity rules.  The table is used at
//! the start of every instruction parse; linear scan cost is dwarfed by the
//! rest of the pipeline.

use super::ast::{AnnotationKind, ArgCount, GateName, MeasureName, NoiseName, TargetArity};

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

impl TableEntry {
    pub fn canonical(&self) -> &'static str {
        match self.kind {
            EntryKind::Gate(n) => n.canonical_name(),
            EntryKind::Noise(n) => n.canonical_name(),
            EntryKind::Measure(n) => n.canonical_name(),
            EntryKind::Annotation(k) => k.canonical_name(),
            EntryKind::MPad => "MPAD",
        }
    }
}

const fn gate(name: GateName, args: ArgCount, targets: TargetArity) -> TableEntry {
    TableEntry {
        kind: EntryKind::Gate(name),
        args,
        targets,
    }
}

const fn noise(name: NoiseName, args: ArgCount, targets: TargetArity) -> TableEntry {
    TableEntry {
        kind: EntryKind::Noise(name),
        args,
        targets,
    }
}

const fn measure(name: MeasureName) -> TableEntry {
    TableEntry {
        kind: EntryKind::Measure(name),
        args: ArgCount::Optional(1),
        targets: TargetArity::AtLeastOne,
    }
}

const fn measure_pairs(name: MeasureName) -> TableEntry {
    TableEntry {
        kind: EntryKind::Measure(name),
        args: ArgCount::Optional(1),
        targets: TargetArity::Pairs,
    }
}

const fn annotation(kind: AnnotationKind) -> TableEntry {
    TableEntry {
        kind: EntryKind::Annotation(kind),
        args: ArgCount::Any,
        targets: TargetArity::Any,
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
    ("R", gate(G::Reset, NoArgs, AtLeastOne)),
    ("RZ", gate(G::ResetZ, NoArgs, AtLeastOne)),
    ("RX", gate(G::ResetX, NoArgs, AtLeastOne)),
    ("RY", gate(G::ResetY, NoArgs, AtLeastOne)),
    // --- Gates: single-qubit Clifford / paulis ---
    ("X", gate(G::X, NoArgs, AtLeastOne)),
    ("Y", gate(G::Y, NoArgs, AtLeastOne)),
    ("Z", gate(G::Z, NoArgs, AtLeastOne)),
    ("H", gate(G::H, NoArgs, AtLeastOne)),
    ("H_XZ", gate(G::HXZ, NoArgs, AtLeastOne)),
    ("S", gate(G::S, NoArgs, AtLeastOne)),
    ("S_DAG", gate(G::SDag, NoArgs, AtLeastOne)),
    ("SQRT_Z", gate(G::SqrtZ, NoArgs, AtLeastOne)),
    ("SQRT_Z_DAG", gate(G::SqrtZDag, NoArgs, AtLeastOne)),
    ("SQRT_X", gate(G::SqrtX, NoArgs, AtLeastOne)),
    ("SQRT_X_DAG", gate(G::SqrtXDag, NoArgs, AtLeastOne)),
    ("SQRT_Y", gate(G::SqrtY, NoArgs, AtLeastOne)),
    ("SQRT_Y_DAG", gate(G::SqrtYDag, NoArgs, AtLeastOne)),
    ("I", gate(G::Identity, NoArgs, AtLeastOne)),
    // --- Gates: two-qubit Clifford ---
    ("CX", gate(G::CX, NoArgs, Pairs)),
    ("ZCX", gate(G::ZCX, NoArgs, Pairs)),
    ("CNOT", gate(G::CNot, NoArgs, Pairs)),
    ("CY", gate(G::CY, NoArgs, Pairs)),
    ("ZCY", gate(G::ZCY, NoArgs, Pairs)),
    ("CZ", gate(G::CZ, NoArgs, Pairs)),
    ("ZCZ", gate(G::ZCZ, NoArgs, Pairs)),
    // --- Gates: phase-1-unsupported (parser accepts) ---
    ("SWAP", gate(G::Swap, NoArgs, Pairs)),
    ("ISWAP", gate(G::ISwap, NoArgs, Pairs)),
    ("ISWAP_DAG", gate(G::ISwapDag, NoArgs, Pairs)),
    ("SQRT_XX", gate(G::SqrtXX, NoArgs, Pairs)),
    ("SQRT_YY", gate(G::SqrtYY, NoArgs, Pairs)),
    ("SQRT_ZZ", gate(G::SqrtZZ, NoArgs, Pairs)),
    ("CXSWAP", gate(G::CXSwap, NoArgs, Pairs)),
    ("SWAPCX", gate(G::SwapCX, NoArgs, Pairs)),
    ("XCX", gate(G::XCX, NoArgs, Pairs)),
    ("XCY", gate(G::XCY, NoArgs, Pairs)),
    ("XCZ", gate(G::XCZ, NoArgs, Pairs)),
    ("YCX", gate(G::YCX, NoArgs, Pairs)),
    ("YCY", gate(G::YCY, NoArgs, Pairs)),
    ("YCZ", gate(G::YCZ, NoArgs, Pairs)),
    ("C_XYZ", gate(G::CXYZ, NoArgs, AtLeastOne)),
    ("C_ZYX", gate(G::CZYX, NoArgs, AtLeastOne)),
    ("H_XY", gate(G::HXY, NoArgs, AtLeastOne)),
    ("H_YZ", gate(G::HYZ, NoArgs, AtLeastOne)),
    // --- Noise ---
    ("DEPOLARIZE1", noise(N::Depolarize1, Exact(1), AtLeastOne)),
    ("DEPOLARIZE2", noise(N::Depolarize2, Exact(1), Pairs)),
    (
        "PAULI_CHANNEL_1",
        noise(N::PauliChannel1, Exact(3), AtLeastOne),
    ),
    ("PAULI_CHANNEL_2", noise(N::PauliChannel2, Exact(15), Pairs)),
    ("X_ERROR", noise(N::XError, Exact(1), AtLeastOne)),
    ("Y_ERROR", noise(N::YError, Exact(1), AtLeastOne)),
    ("Z_ERROR", noise(N::ZError, Exact(1), AtLeastOne)),
    // I_ERROR's arg count varies by tag (`[loss]` => 1, `[correlated_loss]` => 1 or 3).
    // The extended-dialect layer enforces the tag-specific count.
    ("I_ERROR", noise(N::IError, ArgCount::Deferred, AtLeastOne)),
    (
        "HERALDED_ERASE",
        noise(N::HeraldedErase, Exact(1), AtLeastOne),
    ),
    (
        "HERALDED_PAULI_CHANNEL_1",
        noise(N::HeraldedPauliChannel1, Exact(4), AtLeastOne),
    ),
    (
        "CORRELATED_ERROR",
        noise(N::CorrelatedError, Exact(1), AtLeastOne),
    ),
    (
        "ELSE_CORRELATED_ERROR",
        noise(N::ElseCorrelatedError, Exact(1), AtLeastOne),
    ),
    // --- Measurements (all share Optional(1) args) ---
    ("M", measure(Me::M)),
    ("MZ", measure(Me::MZ)),
    ("MR", measure(Me::MR)),
    ("MX", measure(Me::MX)),
    ("MY", measure(Me::MY)),
    ("MRX", measure(Me::MRX)),
    ("MRY", measure(Me::MRY)),
    ("MXX", measure_pairs(Me::MXX)),
    ("MYY", measure_pairs(Me::MYY)),
    ("MZZ", measure_pairs(Me::MZZ)),
    ("MPP", measure(Me::MPP)),
    // --- Annotations ---
    ("DETECTOR", annotation(AnnotationKind::Detector)),
    (
        "MPAD",
        TableEntry {
            kind: EntryKind::MPad,
            args: ArgCount::Optional(1),
            targets: AtLeastOne,
        },
    ),
    (
        "OBSERVABLE_INCLUDE",
        annotation(AnnotationKind::ObservableInclude),
    ),
    ("QUBIT_COORDS", annotation(AnnotationKind::QubitCoords)),
    ("SHIFT_COORDS", annotation(AnnotationKind::ShiftCoords)),
    (
        "TICK",
        TableEntry {
            kind: EntryKind::Annotation(AnnotationKind::Tick),
            args: ArgCount::None,
            targets: TargetArity::Any,
        },
    ),
];

/// Look up a Stim instruction name. `None` means unknown.
pub fn lookup(name: &str) -> Option<TableEntry> {
    TABLE.iter().find(|(n, _)| *n == name).map(|(_, e)| *e)
}

#[cfg(test)]
mod table_tests {
    use super::*;
    use crate::ast::*;

    #[test]
    fn lookup_h_returns_gate_h_with_arity_at_least_one_no_args() {
        let entry = lookup("H").expect("H must be known");
        assert_eq!(
            entry,
            TableEntry {
                kind: EntryKind::Gate(GateName::H),
                args: ArgCount::None,
                targets: TargetArity::AtLeastOne,
            }
        );
    }

    #[test]
    fn lookup_cx_requires_pairs() {
        let entry = lookup("CX").expect("CX must be known");
        assert_eq!(entry.kind, EntryKind::Gate(GateName::CX));
        assert_eq!(entry.targets, TargetArity::Pairs);
    }

    #[test]
    fn lookup_depolarize1_requires_one_arg_any_targets() {
        let entry = lookup("DEPOLARIZE1").expect("DEPOLARIZE1 must be known");
        assert_eq!(entry.kind, EntryKind::Noise(NoiseName::Depolarize1));
        assert_eq!(entry.args, ArgCount::Exact(1));
        assert_eq!(entry.targets, TargetArity::AtLeastOne);
    }

    #[test]
    fn lookup_pauli_channel_2_requires_15_args_pair_targets() {
        let entry = lookup("PAULI_CHANNEL_2").expect("PAULI_CHANNEL_2 must be known");
        assert_eq!(entry.kind, EntryKind::Noise(NoiseName::PauliChannel2));
        assert_eq!(entry.args, ArgCount::Exact(15));
        assert_eq!(entry.targets, TargetArity::Pairs);
    }

    #[test]
    fn lookup_m_returns_measure() {
        let entry = lookup("M").expect("M must be known");
        assert_eq!(entry.kind, EntryKind::Measure(MeasureName::M));
    }

    #[test]
    fn lookup_detector_returns_annotation() {
        let entry = lookup("DETECTOR").expect("DETECTOR must be known");
        assert_eq!(entry.kind, EntryKind::Annotation(AnnotationKind::Detector));
    }

    #[test]
    fn lookup_unknown_returns_none() {
        assert!(lookup("FROBNICATE").is_none());
    }

    #[test]
    fn aliases_map_to_distinct_variants() {
        // CNOT, CX, ZCX are all spelled differently and produce distinct GateName
        // variants — downstream consumers may treat them as the same gate.
        assert_eq!(
            lookup("CNOT").unwrap().kind,
            EntryKind::Gate(GateName::CNot)
        );
        assert_eq!(lookup("CX").unwrap().kind, EntryKind::Gate(GateName::CX));
        assert_eq!(lookup("ZCX").unwrap().kind, EntryKind::Gate(GateName::ZCX));
    }

    #[test]
    fn every_table_key_is_unique() {
        use std::collections::HashSet;
        let mut seen: HashSet<&'static str> = HashSet::new();
        for (key, _) in TABLE {
            assert!(seen.insert(*key), "duplicate key {key:?} in TABLE");
        }
    }
}
