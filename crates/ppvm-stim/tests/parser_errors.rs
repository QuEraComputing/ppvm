use ppvm_stim::{parse, ParseError};

#[test]
fn unknown_instruction_returns_error() {
    let err = parse("FROBNICATE 0").unwrap_err();
    match err {
        ParseError::UnknownInstruction { name, line } => {
            assert_eq!(name, "FROBNICATE");
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn arg_count_mismatch() {
    let err = parse("DEPOLARIZE1(0.1, 0.2) 0").unwrap_err();
    match err {
        ParseError::ArgCount { name, expected, found, line } => {
            assert_eq!(name, "DEPOLARIZE1");
            assert_eq!(expected, 1);
            assert_eq!(found, 2);
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn target_pair_violation() {
    let err = parse("CX 0 1 2").unwrap_err();
    match err {
        ParseError::TargetCount { name, divisor, found, line } => {
            assert_eq!(name, "CX");
            assert_eq!(divisor, 2);
            assert_eq!(found, 3);
            assert_eq!(line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn at_least_one_target_required_for_h() {
    let err = parse("H").unwrap_err();
    assert!(matches!(err, ParseError::TargetCount { .. } | ParseError::Syntax { .. }));
}

#[test]
fn invalid_target_yields_syntax_error() {
    let err = parse("H abc").unwrap_err();
    assert!(matches!(err, ParseError::Syntax { .. }));
}

#[test]
fn unclosed_bracket_yields_syntax_error() {
    let err = parse("S[T 0").unwrap_err();
    assert!(matches!(err, ParseError::Syntax { .. }));
}

#[test]
fn line_numbers_in_errors_are_correct() {
    let err = parse("X 0\nY 0\nFROBNICATE 0").unwrap_err();
    match err {
        ParseError::UnknownInstruction { line, .. } => assert_eq!(line, 3),
        other => panic!("{other:?}"),
    }
}
