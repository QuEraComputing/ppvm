// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! AST-level property test: `parse(print(p))` reaches `p` (modulo
//! source-derived `line` metadata).
//!
//! Where `proptest_roundtrip.rs` starts from generated *strings* and
//! checks the printer fixpoint over what the parser happens to accept,
//! this file starts from generated *valid programs* and checks that
//! print → parse is an inverse pair. Catches printer bugs that the
//! parser-side roundtrip misses (e.g., dropping a tag, mis-grouping
//! correlated-loss pairs, mis-rendering a float).
//!
//! Line numbers are normalised to 0 before comparison: `print` doesn't
//! know the source position the AST originated from, so re-parsed line
//! numbers reflect the printed layout, not the generated one.

use proptest::prelude::*;
use stim_parser::ast::{
    AnnotationKind, GateName, MeasureName, NoiseName, Program, RawInstruction, Target,
};
use stim_parser::extended::parse_extended;
use stim_parser::extended::{Axis, ExtendedInstruction, ExtendedProgram, RawPassthrough};
use stim_parser::prelude::parse;

// ---- generators --------------------------------------------------------

/// Floats that round-trip exactly through the parser/printer (no NaN,
/// no infinity, no `pi` ambiguity, no exponents that look like
/// identifiers). All values printed in canonical decimal form.
fn float_lit() -> impl Strategy<Value = f64> {
    prop_oneof![
        Just(0.0),
        Just(0.25),
        Just(0.5),
        Just(0.75),
        Just(1.0),
        Just(1.5),
        Just(2.0),
        Just(-0.5),
        Just(-1.0),
        Just(0.001),
        Just(0.05),
        Just(0.1),
    ]
}

fn prob_lit() -> impl Strategy<Value = f64> {
    prop_oneof![Just(0.0), Just(0.01), Just(0.05), Just(0.1), Just(0.5)]
}

fn qubit() -> impl Strategy<Value = usize> {
    0usize..8
}

fn single_qubit_clifford() -> impl Strategy<Value = GateName> {
    prop_oneof![
        Just(GateName::H),
        Just(GateName::X),
        Just(GateName::Y),
        Just(GateName::Z),
        Just(GateName::S),
        Just(GateName::SDag),
        Just(GateName::SqrtX),
        Just(GateName::SqrtXDag),
        Just(GateName::SqrtY),
        Just(GateName::SqrtYDag),
        Just(GateName::Identity),
        Just(GateName::Reset),
    ]
}

fn two_qubit_clifford() -> impl Strategy<Value = GateName> {
    prop_oneof![
        Just(GateName::CX),
        Just(GateName::CY),
        Just(GateName::CZ),
        Just(GateName::CNot),
        Just(GateName::Swap),
    ]
}

fn one_q_targets() -> impl Strategy<Value = Vec<usize>> {
    prop::collection::vec(qubit(), 1..4)
}

fn two_q_pair_targets() -> impl Strategy<Value = Vec<usize>> {
    prop::collection::vec(qubit(), 1..3).prop_map(|qs| {
        let mut out = Vec::with_capacity(qs.len() * 2);
        for (i, q) in qs.iter().enumerate() {
            out.push(*q);
            // Force a distinct second qubit so a==b doesn't get rejected
            // by the executor; the parser/printer don't care.
            out.push((q + i + 1) % 8);
        }
        out
    })
}

fn gate_instr() -> impl Strategy<Value = RawInstruction> {
    prop_oneof![
        (single_qubit_clifford(), one_q_targets()).prop_map(|(name, targets)| {
            RawInstruction::Gate {
                name,
                tags: vec![],
                args: vec![],
                targets: targets.into_iter().map(Target::Qubit).collect(),
                line: 0,
            }
        }),
        (two_qubit_clifford(), two_q_pair_targets()).prop_map(|(name, targets)| {
            RawInstruction::Gate {
                name,
                tags: vec![],
                args: vec![],
                targets: targets.into_iter().map(Target::Qubit).collect(),
                line: 0,
            }
        }),
    ]
}

fn noise_instr() -> impl Strategy<Value = RawInstruction> {
    prop_oneof![
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| RawInstruction::Noise {
            name: NoiseName::Depolarize1,
            tags: vec![],
            args: vec![p],
            targets,
            line: 0,
        }),
        (prob_lit(), two_q_pair_targets()).prop_map(|(p, targets)| RawInstruction::Noise {
            name: NoiseName::Depolarize2,
            tags: vec![],
            args: vec![p],
            targets,
            line: 0,
        }),
        (prob_lit(), prob_lit(), prob_lit(), one_q_targets()).prop_map(|(a, b, c, targets)| {
            RawInstruction::Noise {
                name: NoiseName::PauliChannel1,
                tags: vec![],
                args: vec![a, b, c],
                targets,
                line: 0,
            }
        }),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| RawInstruction::Noise {
            name: NoiseName::XError,
            tags: vec![],
            args: vec![p],
            targets,
            line: 0,
        }),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| RawInstruction::Noise {
            name: NoiseName::YError,
            tags: vec![],
            args: vec![p],
            targets,
            line: 0,
        }),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| RawInstruction::Noise {
            name: NoiseName::ZError,
            tags: vec![],
            args: vec![p],
            targets,
            line: 0,
        }),
    ]
}

fn measure_instr() -> impl Strategy<Value = RawInstruction> {
    let name = prop_oneof![
        Just(MeasureName::M),
        Just(MeasureName::MZ),
        Just(MeasureName::MR),
    ];
    (name, proptest::option::of(prob_lit()), one_q_targets()).prop_map(|(name, args, targets)| {
        RawInstruction::Measure {
            name,
            tags: vec![],
            args: args.map(|p| vec![p]).unwrap_or_default(),
            targets,
            line: 0,
        }
    })
}

fn annotation_instr() -> impl Strategy<Value = RawInstruction> {
    prop_oneof![
        Just(RawInstruction::Annotation {
            kind: AnnotationKind::Tick,
            args: vec![],
            targets: vec![],
            line: 0,
        }),
        Just(RawInstruction::Annotation {
            kind: AnnotationKind::Detector,
            args: vec![],
            targets: vec![],
            line: 0,
        }),
    ]
}

fn mpad_instr() -> impl Strategy<Value = RawInstruction> {
    (
        proptest::option::of(prob_lit()),
        prop::collection::vec(0usize..2, 1..6),
    )
        .prop_map(|(prob, bits)| RawInstruction::MPad {
            tags: vec![],
            prob,
            bits,
            line: 0,
        })
}

fn flat_instr() -> impl Strategy<Value = RawInstruction> {
    prop_oneof![
        gate_instr(),
        noise_instr(),
        measure_instr(),
        annotation_instr(),
        mpad_instr(),
    ]
}

fn instr() -> impl Strategy<Value = RawInstruction> {
    // Depth 1, ≤6 nodes total, ≤2 branches expected. Generator stack
    // stays shallow; the body itself is `flat_instr()` (non-recursive).
    flat_instr().prop_recursive(1, 6, 2, |_inner| {
        (1u64..5, prop::collection::vec(flat_instr(), 1..3)).prop_map(|(count, body)| {
            RawInstruction::Repeat {
                count,
                body,
                line: 0,
            }
        })
    })
}

fn program() -> impl Strategy<Value = Program> {
    prop::collection::vec(instr(), 0..10).prop_map(|instructions| Program { instructions })
}

// ---- extended-AST generators ------------------------------------------

/// `RawPassthrough` generators mirror the raw-AST generators for the four
/// variants that survive the interpret pass unchanged.
fn ext_gate_instr() -> impl Strategy<Value = RawPassthrough> {
    prop_oneof![
        (single_qubit_clifford(), one_q_targets()).prop_map(|(name, targets)| {
            RawPassthrough::Gate {
                name,
                tags: vec![],
                args: vec![],
                targets: targets.into_iter().map(Target::Qubit).collect(),
                line: 0,
            }
        }),
        (two_qubit_clifford(), two_q_pair_targets()).prop_map(|(name, targets)| {
            RawPassthrough::Gate {
                name,
                tags: vec![],
                args: vec![],
                targets: targets.into_iter().map(Target::Qubit).collect(),
                line: 0,
            }
        }),
    ]
}

fn ext_noise_instr() -> impl Strategy<Value = RawPassthrough> {
    prop_oneof![
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| RawPassthrough::Noise {
            name: NoiseName::Depolarize1,
            tags: vec![],
            args: vec![p],
            targets,
            line: 0,
        }),
        (prob_lit(), two_q_pair_targets()).prop_map(|(p, targets)| RawPassthrough::Noise {
            name: NoiseName::Depolarize2,
            tags: vec![],
            args: vec![p],
            targets,
            line: 0,
        }),
        (prob_lit(), prob_lit(), prob_lit(), one_q_targets()).prop_map(|(a, b, c, targets)| {
            RawPassthrough::Noise {
                name: NoiseName::PauliChannel1,
                tags: vec![],
                args: vec![a, b, c],
                targets,
                line: 0,
            }
        }),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| RawPassthrough::Noise {
            name: NoiseName::XError,
            tags: vec![],
            args: vec![p],
            targets,
            line: 0,
        }),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| RawPassthrough::Noise {
            name: NoiseName::YError,
            tags: vec![],
            args: vec![p],
            targets,
            line: 0,
        }),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| RawPassthrough::Noise {
            name: NoiseName::ZError,
            tags: vec![],
            args: vec![p],
            targets,
            line: 0,
        }),
    ]
}

fn ext_measure_instr() -> impl Strategy<Value = RawPassthrough> {
    let name = prop_oneof![
        Just(MeasureName::M),
        Just(MeasureName::MZ),
        Just(MeasureName::MR),
    ];
    (name, proptest::option::of(prob_lit()), one_q_targets()).prop_map(|(name, args, targets)| {
        RawPassthrough::Measure {
            name,
            tags: vec![],
            args: args.map(|p| vec![p]).unwrap_or_default(),
            targets,
            line: 0,
        }
    })
}

fn ext_annotation_instr() -> impl Strategy<Value = RawPassthrough> {
    prop_oneof![
        Just(RawPassthrough::Annotation {
            kind: AnnotationKind::Tick,
            args: vec![],
            targets: vec![],
            line: 0,
        }),
        Just(RawPassthrough::Annotation {
            kind: AnnotationKind::Detector,
            args: vec![],
            targets: vec![],
            line: 0,
        }),
    ]
}

fn ext_raw() -> impl Strategy<Value = RawPassthrough> {
    prop_oneof![
        ext_gate_instr(),
        ext_noise_instr(),
        ext_measure_instr(),
        ext_annotation_instr(),
    ]
}

fn ext_flat() -> impl Strategy<Value = ExtendedInstruction> {
    prop_oneof![
        ext_raw().prop_map(ExtendedInstruction::Raw),
        one_q_targets().prop_map(|targets| ExtendedInstruction::T { targets, line: 0 }),
        one_q_targets().prop_map(|targets| ExtendedInstruction::TDag { targets, line: 0 }),
        (
            prop_oneof![Just(Axis::X), Just(Axis::Y), Just(Axis::Z)],
            float_lit(),
            one_q_targets(),
        )
            .prop_map(|(axis, theta, targets)| ExtendedInstruction::Rotation {
                axis,
                theta,
                targets,
                line: 0,
            }),
        (float_lit(), float_lit(), float_lit(), one_q_targets()).prop_map(
            |(theta, phi, lambda, targets)| ExtendedInstruction::U3 {
                theta,
                phi,
                lambda,
                targets,
                line: 0,
            }
        ),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| ExtendedInstruction::Loss {
            p,
            targets,
            line: 0,
        }),
        (
            prob_lit(),
            prob_lit(),
            prob_lit(),
            prop::collection::vec((qubit(), qubit()), 1..3),
        )
            .prop_map(
                |(p0, p1, p2, targets)| ExtendedInstruction::CorrelatedLoss {
                    ps: [p0, p1, p2],
                    targets,
                    line: 0,
                }
            ),
        (
            proptest::option::of(prob_lit()),
            prop::collection::vec(any::<bool>(), 1..6),
        )
            .prop_map(|(prob, bits)| ExtendedInstruction::MPad {
                tags: vec![],
                prob,
                bits,
                line: 0,
            }),
    ]
}

fn ext_instr() -> impl Strategy<Value = ExtendedInstruction> {
    ext_flat().prop_recursive(1, 6, 2, |_inner| {
        (1u64..5, prop::collection::vec(ext_flat(), 1..3)).prop_map(|(count, body)| {
            ExtendedInstruction::Repeat {
                count,
                body,
                line: 0,
            }
        })
    })
}

fn ext_program() -> impl Strategy<Value = ExtendedProgram> {
    prop::collection::vec(ext_instr(), 0..10)
        .prop_map(|instructions| ExtendedProgram { instructions })
}

// ---- line-number normalisers ------------------------------------------

fn zero_lines_raw(instrs: &mut [RawInstruction]) {
    for i in instrs {
        match i {
            RawInstruction::Gate { line, .. }
            | RawInstruction::Noise { line, .. }
            | RawInstruction::Measure { line, .. }
            | RawInstruction::Annotation { line, .. }
            | RawInstruction::MPad { line, .. }
            | RawInstruction::Mpp { line, .. } => *line = 0,
            RawInstruction::Repeat { body, line, .. } => {
                *line = 0;
                zero_lines_raw(body);
            }
        }
    }
}

fn zero_lines_raw_passthrough(r: &mut RawPassthrough) {
    match r {
        RawPassthrough::Gate { line, .. }
        | RawPassthrough::Noise { line, .. }
        | RawPassthrough::Measure { line, .. }
        | RawPassthrough::Annotation { line, .. } => *line = 0,
    }
}

fn zero_lines_ext(instrs: &mut [ExtendedInstruction]) {
    for i in instrs {
        match i {
            ExtendedInstruction::Raw(r) => zero_lines_raw_passthrough(r),
            ExtendedInstruction::T { line, .. }
            | ExtendedInstruction::TDag { line, .. }
            | ExtendedInstruction::Rotation { line, .. }
            | ExtendedInstruction::U3 { line, .. }
            | ExtendedInstruction::Loss { line, .. }
            | ExtendedInstruction::CorrelatedLoss { line, .. }
            | ExtendedInstruction::MPad { line, .. }
            | ExtendedInstruction::Mpp { line, .. } => *line = 0,
            ExtendedInstruction::Repeat { body, line, .. } => {
                *line = 0;
                zero_lines_ext(body);
            }
        }
    }
}

// ---- properties --------------------------------------------------------

proptest! {
    /// For every generated valid program, the printer/parser round-trip
    /// must reproduce the program (with line numbers zeroed).
    #[test]
    fn raw_print_parse_roundtrip(mut prog in program()) {
        zero_lines_raw(&mut prog.instructions);
        let printed = format!("{prog}");
        let mut reparsed = parse(&printed).unwrap_or_else(|e| {
            panic!("reparse failed for generated program: {e}\n--printed--\n{printed}")
        });
        zero_lines_raw(&mut reparsed.instructions);
        prop_assert_eq!(prog, reparsed, "raw print/parse not inverse on:\n{}", printed);
    }

    /// Same property for the extended AST: print → parse_extended must
    /// recover the original `ExtendedProgram`.
    #[test]
    fn extended_print_parse_roundtrip(mut prog in ext_program()) {
        zero_lines_ext(&mut prog.instructions);
        let printed = format!("{prog}");
        let mut reparsed = parse_extended(&printed).unwrap_or_else(|e| {
            panic!("reparse_extended failed for generated program: {e}\n--printed--\n{printed}")
        });
        zero_lines_ext(&mut reparsed.instructions);
        prop_assert_eq!(prog, reparsed, "extended print/parse not inverse on:\n{}", printed);
    }
}
