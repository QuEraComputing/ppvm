// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

use stim_parser::prelude::*;

#[test]
fn parse_m_targets() {
    let p = parse("M 0 1 2").unwrap();
    match &p.instructions[0] {
        RawInstruction::Measure { name, targets, .. } => {
            assert_eq!(*name, MeasureName::M);
            assert_eq!(targets, &[0, 1, 2]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_mz_alias() {
    let p = parse("MZ 5").unwrap();
    let RawInstruction::Measure { name, .. } = &p.instructions[0] else {
        panic!()
    };
    assert_eq!(*name, MeasureName::MZ);
}

#[test]
fn parse_mr_alias() {
    let p = parse("MR 5").unwrap();
    let RawInstruction::Measure { name, .. } = &p.instructions[0] else {
        panic!()
    };
    assert_eq!(*name, MeasureName::MR);
}

#[test]
fn parse_annotation_tick() {
    let p = parse("TICK").unwrap();
    let RawInstruction::Annotation { kind, .. } = &p.instructions[0] else {
        panic!()
    };
    assert_eq!(*kind, AnnotationKind::Tick);
}

#[test]
fn parse_annotation_tick_with_args_rejected() {
    let err = parse("TICK(0.1)").unwrap_err();
    match err {
        ParseError::ArgCount {
            name,
            expected,
            found,
            ..
        } => {
            assert_eq!(name, "TICK");
            assert_eq!(expected, 0);
            assert_eq!(found, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_annotation_with_args_and_targets() {
    let p = parse("QUBIT_COORDS(0, 0) 0").unwrap();
    let RawInstruction::Annotation {
        kind,
        args,
        targets,
        ..
    } = &p.instructions[0]
    else {
        panic!()
    };
    assert_eq!(*kind, AnnotationKind::QubitCoords);
    assert_eq!(args.len(), 2);
    assert_eq!(targets, &[0]);
}

#[test]
fn parse_detector_no_targets() {
    let p = parse("DETECTOR").unwrap();
    let RawInstruction::Annotation { kind, targets, .. } = &p.instructions[0] else {
        panic!()
    };
    assert_eq!(*kind, AnnotationKind::Detector);
    assert!(targets.is_empty());
}

#[test]
fn parse_observable_include_with_paren_arg() {
    let p = parse("OBSERVABLE_INCLUDE(0)").unwrap();
    let RawInstruction::Annotation { kind, args, .. } = &p.instructions[0] else {
        panic!()
    };
    assert_eq!(*kind, AnnotationKind::ObservableInclude);
    assert_eq!(args.len(), 1);
}

#[test]
fn parse_m_with_noise_arg() {
    let p = parse("M(0.001) 0 1 2").unwrap();
    match &p.instructions[0] {
        RawInstruction::Measure {
            name,
            args,
            targets,
            ..
        } => {
            assert_eq!(*name, MeasureName::M);
            assert_eq!(args.len(), 1);
            assert!((args[0] - 0.001).abs() < 1e-12);
            assert_eq!(targets, &[0, 1, 2]);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_mz_with_noise_arg() {
    let p = parse("MZ(0.5) 5").unwrap();
    let RawInstruction::Measure { name, args, .. } = &p.instructions[0] else {
        panic!()
    };
    assert_eq!(*name, MeasureName::MZ);
    assert_eq!(args.len(), 1);
    assert!((args[0] - 0.5).abs() < 1e-12);
}

#[test]
fn parse_mr_with_noise_arg() {
    let p = parse("MR(0.01) 0").unwrap();
    let RawInstruction::Measure { name, args, .. } = &p.instructions[0] else {
        panic!()
    };
    assert_eq!(*name, MeasureName::MR);
    assert_eq!(args[0], 0.01);
}

#[test]
fn parse_m_without_noise_still_works() {
    let p = parse("M 0").unwrap();
    let RawInstruction::Measure { args, .. } = &p.instructions[0] else {
        panic!()
    };
    assert!(args.is_empty());
}

#[test]
fn parse_m_with_two_args_rejected() {
    let err = parse("M(0.1, 0.2) 0").unwrap_err();
    match err {
        ParseError::ArgCount {
            name,
            expected,
            found,
            ..
        } => {
            assert_eq!(name, "M");
            assert_eq!(expected, 1);
            assert_eq!(found, 2);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_mpad_no_args_no_tags() {
    let p = parse("MPAD 0 1 0").unwrap();
    let RawInstruction::MPad {
        tags, prob, bits, ..
    } = &p.instructions[0]
    else {
        panic!("{:?}", p.instructions[0]);
    };
    assert!(tags.is_empty());
    assert_eq!(*prob, None);
    assert_eq!(bits, &[0usize, 1, 0]);
}

#[test]
fn parse_mpad_with_prob() {
    let p = parse("MPAD(0.25) 0 1").unwrap();
    let RawInstruction::MPad { prob, bits, .. } = &p.instructions[0] else {
        panic!()
    };
    assert_eq!(*prob, Some(0.25));
    assert_eq!(bits, &[0usize, 1]);
}
