#[allow(unused_imports)]
use ppvm_stim::{parse, RawInstruction, GateName, AnnotationKind};

#[test]
fn empty_program_parses_to_empty_instructions() {
    let p = parse("").unwrap();
    assert!(p.instructions.is_empty());
}

#[test]
fn whitespace_only_program() {
    let p = parse("   \n\n\t\n").unwrap();
    assert!(p.instructions.is_empty());
}

#[test]
fn comments_and_blank_lines_skipped() {
    let src = "
# header comment

X 0
# mid comment
M 0
";
    let p = parse(src).unwrap();
    assert_eq!(p.instructions.len(), 2);
    let RawInstruction::Gate { name, line, .. } = &p.instructions[0] else {
        panic!()
    };
    assert_eq!(*name, GateName::X);
    assert_eq!(*line, 4); // counting the leading blank line
}

#[test]
fn line_numbers_are_one_indexed_and_count_blank_and_comment_lines() {
    let src = "X 0\n# c\nY 0\n\nZ 0";
    let p = parse(src).unwrap();
    assert_eq!(p.instructions.len(), 3);
    let lines: Vec<usize> = p
        .instructions
        .iter()
        .map(|i| match i {
            RawInstruction::Gate { line, .. } => *line,
            _ => panic!(),
        })
        .collect();
    assert_eq!(lines, vec![1, 3, 5]);
}

#[test]
fn leading_indentation_tolerated() {
    let p = parse("    H 0\n\tH 1").unwrap();
    assert_eq!(p.instructions.len(), 2);
}

#[test]
fn parse_simple_repeat() {
    let src = "REPEAT 3 {\n    X 0\n    M 0\n}";
    let p = parse(src).unwrap();
    assert_eq!(p.instructions.len(), 1);
    match &p.instructions[0] {
        RawInstruction::Repeat { count, body, line } => {
            assert_eq!(*count, 3);
            assert_eq!(body.len(), 2);
            assert_eq!(*line, 1);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn parse_nested_repeat() {
    let src = "REPEAT 2 {\n    REPEAT 3 {\n        H 0\n    }\n}";
    let p = parse(src).unwrap();
    let RawInstruction::Repeat { body, .. } = &p.instructions[0] else { panic!() };
    assert!(matches!(body[0], RawInstruction::Repeat { count: 3, .. }));
}

#[test]
fn parse_repeat_then_following_instruction() {
    let src = "REPEAT 2 { H 0 }\nM 0";
    let p = parse(src).unwrap();
    assert_eq!(p.instructions.len(), 2);
    assert!(matches!(p.instructions[0], RawInstruction::Repeat { .. }));
    assert!(matches!(p.instructions[1], RawInstruction::Measure { .. }));
}

#[test]
fn parse_repeat_one_line() {
    // Single-line REPEAT works because the tokenizer splits `{` and `}`
    // out wherever they appear, not just at line boundaries.
    let p = parse("REPEAT 5 { H 0 }").unwrap();
    let RawInstruction::Repeat { count, body, .. } = &p.instructions[0] else {
        panic!("expected Repeat, got {:?}", &p.instructions[0]);
    };
    assert_eq!(*count, 5);
    assert_eq!(body.len(), 1);
    assert!(matches!(body[0], RawInstruction::Gate { name: GateName::H, .. }));
}
