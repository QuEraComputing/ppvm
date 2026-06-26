// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! AST-level property test: `parse(print(p))` reaches `p` (modulo
//! source-derived `span` metadata).
//!
//! Where `proptest_roundtrip.rs` starts from generated *strings* and
//! checks the printer fixpoint over what the parser happens to accept,
//! this file starts from generated *valid programs* and checks that
//! print → parse is an inverse pair. Catches printer bugs that the
//! parser-side roundtrip misses (e.g., dropping a tag, mis-grouping
//! correlated-loss pairs, mis-rendering a float).
//!
//! Spans are zeroed before comparison: `print` doesn't know the source
//! position the AST originated from, so re-parsed spans reflect the
//! printed layout, not the generated one.

use std::sync::Arc;

use proptest::prelude::*;
use stim_parser::diagnostics::{LineMap, Span};
use stim_parser::prelude::{
    AnnotationKind, AnnotationOp, Axis, ExtendedInstruction, ExtendedProgram, GateName, GateOp,
    Instruction, MeasureName, MeasureOp, NoiseName, NoiseOp, Program, Target, parse,
    parse_extended,
};

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

fn span0() -> Span {
    Span::new(0, 0)
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

fn gate_instr() -> impl Strategy<Value = Instruction> {
    prop_oneof![
        (single_qubit_clifford(), one_q_targets()).prop_map(|(name, targets)| {
            Instruction::Gate(GateOp {
                name,
                tags: vec![],
                args: vec![],
                targets: targets.into_iter().map(Target::Qubit).collect(),
                span: span0(),
            })
        }),
        (two_qubit_clifford(), two_q_pair_targets()).prop_map(|(name, targets)| {
            Instruction::Gate(GateOp {
                name,
                tags: vec![],
                args: vec![],
                targets: targets.into_iter().map(Target::Qubit).collect(),
                span: span0(),
            })
        }),
    ]
}

fn noise_instr() -> impl Strategy<Value = Instruction> {
    prop_oneof![
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| Instruction::Noise(NoiseOp {
            name: NoiseName::Depolarize1,
            tags: vec![],
            args: vec![p],
            targets,
            span: span0(),
        })),
        (prob_lit(), two_q_pair_targets()).prop_map(|(p, targets)| Instruction::Noise(NoiseOp {
            name: NoiseName::Depolarize2,
            tags: vec![],
            args: vec![p],
            targets,
            span: span0(),
        })),
        (prob_lit(), prob_lit(), prob_lit(), one_q_targets()).prop_map(|(a, b, c, targets)| {
            Instruction::Noise(NoiseOp {
                name: NoiseName::PauliChannel1,
                tags: vec![],
                args: vec![a, b, c],
                targets,
                span: span0(),
            })
        }),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| Instruction::Noise(NoiseOp {
            name: NoiseName::XError,
            tags: vec![],
            args: vec![p],
            targets,
            span: span0(),
        })),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| Instruction::Noise(NoiseOp {
            name: NoiseName::YError,
            tags: vec![],
            args: vec![p],
            targets,
            span: span0(),
        })),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| Instruction::Noise(NoiseOp {
            name: NoiseName::ZError,
            tags: vec![],
            args: vec![p],
            targets,
            span: span0(),
        })),
    ]
}

fn measure_instr() -> impl Strategy<Value = Instruction> {
    let name = prop_oneof![
        Just(MeasureName::M),
        Just(MeasureName::MZ),
        Just(MeasureName::MR),
    ];
    (name, proptest::option::of(prob_lit()), one_q_targets()).prop_map(|(name, args, targets)| {
        Instruction::Measure(MeasureOp {
            name,
            tags: vec![],
            args: args.map(|p| vec![p]).unwrap_or_default(),
            targets,
            span: span0(),
        })
    })
}

fn annotation_instr() -> impl Strategy<Value = Instruction> {
    prop_oneof![
        Just(Instruction::Annotation(AnnotationOp {
            kind: AnnotationKind::Tick,
            args: vec![],
            targets: vec![],
            span: span0(),
        })),
        Just(Instruction::Annotation(AnnotationOp {
            kind: AnnotationKind::Detector,
            args: vec![],
            targets: vec![],
            span: span0(),
        })),
    ]
}

fn mpad_instr() -> impl Strategy<Value = Instruction> {
    (
        proptest::option::of(prob_lit()),
        prop::collection::vec(0usize..2, 1..6),
    )
        .prop_map(|(prob, bits)| Instruction::MPad {
            tags: vec![],
            prob,
            bits,
            span: span0(),
        })
}

fn flat_instr() -> impl Strategy<Value = Instruction> {
    prop_oneof![
        gate_instr(),
        noise_instr(),
        measure_instr(),
        annotation_instr(),
        mpad_instr(),
    ]
}

fn instr() -> impl Strategy<Value = Instruction> {
    // Depth 1, ≤6 nodes total, ≤2 branches expected. Generator stack
    // stays shallow; the body itself is `flat_instr()` (non-recursive).
    flat_instr().prop_recursive(1, 6, 2, |_inner| {
        (1u64..5, prop::collection::vec(flat_instr(), 1..3)).prop_map(|(count, body)| {
            Instruction::Repeat {
                count,
                body,
                span: span0(),
            }
        })
    })
}

fn program() -> impl Strategy<Value = Program> {
    prop::collection::vec(instr(), 0..10).prop_map(|instructions| Program {
        instructions,
        line_map: Arc::new(LineMap::new("")),
    })
}

// ---- extended-AST generators ------------------------------------------

fn ext_gate_instr() -> impl Strategy<Value = ExtendedInstruction> {
    prop_oneof![
        (single_qubit_clifford(), one_q_targets()).prop_map(|(name, targets)| {
            ExtendedInstruction::Gate(GateOp {
                name,
                tags: vec![],
                args: vec![],
                targets: targets.into_iter().map(Target::Qubit).collect(),
                span: span0(),
            })
        }),
        (two_qubit_clifford(), two_q_pair_targets()).prop_map(|(name, targets)| {
            ExtendedInstruction::Gate(GateOp {
                name,
                tags: vec![],
                args: vec![],
                targets: targets.into_iter().map(Target::Qubit).collect(),
                span: span0(),
            })
        }),
    ]
}

fn ext_noise_instr() -> impl Strategy<Value = ExtendedInstruction> {
    prop_oneof![
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| ExtendedInstruction::Noise(
            NoiseOp {
                name: NoiseName::Depolarize1,
                tags: vec![],
                args: vec![p],
                targets,
                span: span0(),
            }
        )),
        (prob_lit(), two_q_pair_targets()).prop_map(|(p, targets)| ExtendedInstruction::Noise(
            NoiseOp {
                name: NoiseName::Depolarize2,
                tags: vec![],
                args: vec![p],
                targets,
                span: span0(),
            }
        )),
        (prob_lit(), prob_lit(), prob_lit(), one_q_targets()).prop_map(|(a, b, c, targets)| {
            ExtendedInstruction::Noise(NoiseOp {
                name: NoiseName::PauliChannel1,
                tags: vec![],
                args: vec![a, b, c],
                targets,
                span: span0(),
            })
        }),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| ExtendedInstruction::Noise(
            NoiseOp {
                name: NoiseName::XError,
                tags: vec![],
                args: vec![p],
                targets,
                span: span0(),
            }
        )),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| ExtendedInstruction::Noise(
            NoiseOp {
                name: NoiseName::YError,
                tags: vec![],
                args: vec![p],
                targets,
                span: span0(),
            }
        )),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| ExtendedInstruction::Noise(
            NoiseOp {
                name: NoiseName::ZError,
                tags: vec![],
                args: vec![p],
                targets,
                span: span0(),
            }
        )),
    ]
}

fn ext_measure_instr() -> impl Strategy<Value = ExtendedInstruction> {
    let name = prop_oneof![
        Just(MeasureName::M),
        Just(MeasureName::MZ),
        Just(MeasureName::MR),
    ];
    (name, proptest::option::of(prob_lit()), one_q_targets()).prop_map(|(name, args, targets)| {
        ExtendedInstruction::Measure(MeasureOp {
            name,
            tags: vec![],
            args: args.map(|p| vec![p]).unwrap_or_default(),
            targets,
            span: span0(),
        })
    })
}

fn ext_annotation_instr() -> impl Strategy<Value = ExtendedInstruction> {
    prop_oneof![
        Just(ExtendedInstruction::Annotation(AnnotationOp {
            kind: AnnotationKind::Tick,
            args: vec![],
            targets: vec![],
            span: span0(),
        })),
        Just(ExtendedInstruction::Annotation(AnnotationOp {
            kind: AnnotationKind::Detector,
            args: vec![],
            targets: vec![],
            span: span0(),
        })),
    ]
}

fn ext_raw() -> impl Strategy<Value = ExtendedInstruction> {
    prop_oneof![
        ext_gate_instr(),
        ext_noise_instr(),
        ext_measure_instr(),
        ext_annotation_instr(),
    ]
}

fn ext_flat() -> impl Strategy<Value = ExtendedInstruction> {
    prop_oneof![
        ext_raw(),
        one_q_targets().prop_map(|targets| ExtendedInstruction::T {
            targets,
            span: span0()
        }),
        one_q_targets().prop_map(|targets| ExtendedInstruction::TDag {
            targets,
            span: span0()
        }),
        (
            prop_oneof![Just(Axis::X), Just(Axis::Y), Just(Axis::Z)],
            float_lit(),
            one_q_targets(),
        )
            .prop_map(|(axis, theta, targets)| ExtendedInstruction::Rotation {
                axis,
                // Rotation angles are stored in radians but always originate from
                // a `<n>*pi` tag, so only half-turn multiples are parser-producible.
                theta: theta * std::f64::consts::PI,
                targets,
                span: span0(),
            }),
        (float_lit(), float_lit(), float_lit(), one_q_targets()).prop_map(
            |(theta, phi, lambda, targets)| ExtendedInstruction::U3 {
                theta: theta * std::f64::consts::PI,
                phi: phi * std::f64::consts::PI,
                lambda: lambda * std::f64::consts::PI,
                targets,
                span: span0(),
            }
        ),
        (prob_lit(), one_q_targets()).prop_map(|(p, targets)| ExtendedInstruction::Loss {
            p,
            targets,
            span: span0(),
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
                    span: span0(),
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
                span: span0(),
            }),
    ]
}

fn ext_instr() -> impl Strategy<Value = ExtendedInstruction> {
    ext_flat().prop_recursive(1, 6, 2, |_inner| {
        (1u64..5, prop::collection::vec(ext_flat(), 1..3)).prop_map(|(count, body)| {
            ExtendedInstruction::Repeat {
                count,
                body,
                span: span0(),
            }
        })
    })
}

fn ext_program() -> impl Strategy<Value = ExtendedProgram> {
    prop::collection::vec(ext_instr(), 0..10).prop_map(|instructions| ExtendedProgram {
        instructions,
        line_map: Arc::new(LineMap::new("")),
    })
}

// ---- span normalisers ------------------------------------------

fn zero_spans_raw(instrs: &mut [Instruction]) {
    for i in instrs {
        match i {
            Instruction::Gate(op) => op.span = span0(),
            Instruction::Noise(op) => op.span = span0(),
            Instruction::Measure(op) => op.span = span0(),
            Instruction::Annotation(op) => op.span = span0(),
            Instruction::Mpp(op) => op.span = span0(),
            Instruction::MPad { span, .. } => *span = span0(),
            Instruction::Repeat { body, span, .. } => {
                *span = span0();
                zero_spans_raw(body);
            }
        }
    }
}

fn zero_spans_ext(instrs: &mut [ExtendedInstruction]) {
    for i in instrs {
        match i {
            ExtendedInstruction::Gate(op) => op.span = span0(),
            ExtendedInstruction::Noise(op) => op.span = span0(),
            ExtendedInstruction::Measure(op) => op.span = span0(),
            ExtendedInstruction::Annotation(op) => op.span = span0(),
            ExtendedInstruction::Mpp(op) => op.span = span0(),
            ExtendedInstruction::T { span, .. }
            | ExtendedInstruction::TDag { span, .. }
            | ExtendedInstruction::Rotation { span, .. }
            | ExtendedInstruction::U3 { span, .. }
            | ExtendedInstruction::Loss { span, .. }
            | ExtendedInstruction::CorrelatedLoss { span, .. }
            | ExtendedInstruction::MPad { span, .. } => *span = span0(),
            ExtendedInstruction::Repeat { body, span, .. } => {
                *span = span0();
                zero_spans_ext(body);
            }
        }
    }
}

// ---- properties --------------------------------------------------------

proptest! {
    /// For every generated valid program, the printer/parser round-trip
    /// must reproduce the program (with spans zeroed).
    #[test]
    fn raw_print_parse_roundtrip(mut prog in program()) {
        zero_spans_raw(&mut prog.instructions);
        let printed = format!("{prog}");
        let mut reparsed = parse(&printed).unwrap_or_else(|e| {
            panic!("reparse failed for generated program: {e}\n--printed--\n{printed}")
        });
        zero_spans_raw(&mut reparsed.instructions);
        prop_assert_eq!(prog, reparsed, "raw print/parse not inverse on:\n{}", printed);
    }

    /// Same property for the extended AST: print → parse_extended must
    /// recover the original `ExtendedProgram`.
    #[test]
    fn extended_print_parse_roundtrip(mut prog in ext_program()) {
        zero_spans_ext(&mut prog.instructions);
        let printed = format!("{prog}");
        let mut reparsed = parse_extended(&printed).unwrap_or_else(|e| {
            panic!("reparse_extended failed for generated program: {e}\n--printed--\n{printed}")
        });
        zero_spans_ext(&mut reparsed.instructions);
        prop_assert_eq!(prog, reparsed, "extended print/parse not inverse on:\n{}", printed);
    }
}
