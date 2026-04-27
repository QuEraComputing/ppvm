use ppvm_stim::{parse, RawInstruction, MeasureName, AnnotationKind};

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
    let RawInstruction::Measure { name, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(*name, MeasureName::MZ);
}

#[test]
fn parse_mr_alias() {
    let p = parse("MR 5").unwrap();
    let RawInstruction::Measure { name, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(*name, MeasureName::MR);
}

#[test]
fn parse_annotation_tick() {
    let p = parse("TICK").unwrap();
    let RawInstruction::Annotation { kind, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(*kind, AnnotationKind::Tick);
}

#[test]
fn parse_annotation_with_args_and_targets() {
    let p = parse("QUBIT_COORDS(0, 0) 0").unwrap();
    let RawInstruction::Annotation { kind, args, targets, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(*kind, AnnotationKind::QubitCoords);
    assert_eq!(args.len(), 2);
    assert_eq!(targets, &[0]);
}

#[test]
fn parse_detector_no_targets() {
    let p = parse("DETECTOR").unwrap();
    let RawInstruction::Annotation { kind, targets, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(*kind, AnnotationKind::Detector);
    assert!(targets.is_empty());
}

#[test]
fn parse_observable_include_with_paren_arg() {
    let p = parse("OBSERVABLE_INCLUDE(0)").unwrap();
    let RawInstruction::Annotation { kind, args, .. } = &p.instructions[0] else { panic!() };
    assert_eq!(*kind, AnnotationKind::ObservableInclude);
    assert_eq!(args.len(), 1);
}
