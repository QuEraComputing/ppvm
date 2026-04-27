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
