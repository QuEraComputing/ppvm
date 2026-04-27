pub mod ast;
pub mod table;

use ast::{
    AnnotationKind, ArgCount, GateName, MeasureName, NoiseName,
    ParseError, Program, RawInstruction, TargetArity,
};
use table::{TableEntry, lookup};

/// Map a byte offset in `src` to a 1-indexed line number.
struct LineMap {
    /// `starts[i]` = byte offset of the start of line (i+1).
    starts: Vec<usize>,
}

impl LineMap {
    fn new(src: &str) -> Self {
        let mut starts = vec![0];
        for (i, b) in src.bytes().enumerate() {
            if b == b'\n' {
                starts.push(i + 1);
            }
        }
        Self { starts }
    }

    #[allow(dead_code)] // used by chumsky integration in a future refactor
    fn line_of(&self, byte_offset: usize) -> usize {
        // Largest start <= byte_offset.
        match self.starts.binary_search(&byte_offset) {
            Ok(i) => i + 1,
            Err(i) => i, // i is the insertion index; start of line `i` is at starts[i-1].
        }
    }
}

/// Parse Stim source into a [`Program`].
pub fn parse(src: &str) -> Result<Program, ParseError> {
    let line_map = LineMap::new(src);
    let mut instructions = Vec::new();

    for (line_idx, raw_line) in src.lines().enumerate() {
        let line_no = line_idx + 1;
        let trimmed = raw_line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let instr = parse_line(trimmed, line_no, &line_map)?;
        instructions.push(instr);
    }

    Ok(Program { instructions })
}

fn parse_line(line: &str, line_no: usize, _line_map: &LineMap) -> Result<RawInstruction, ParseError> {
    // Split off the instruction token: name (possibly followed by `[…]` and/or `(…)`)
    // separated from targets by the first whitespace at bracket-depth 0.
    let (head, targets_part) = split_head_and_targets(line);
    let (name_str, _tags_src, _args_src) = parse_head(head, line_no)?;

    let entry = lookup(name_str).ok_or(ParseError::UnknownInstruction {
        name: name_str.to_string(),
        line: line_no,
    })?;

    let targets: Vec<usize> = targets_part
        .split_whitespace()
        .map(|t| t.parse::<usize>().map_err(|_| ParseError::Syntax {
            line: line_no,
            col: 1,
            message: format!("invalid target {t:?}"),
        }))
        .collect::<Result<_, _>>()?;

    Ok(build_instruction(entry, name_str, vec![], vec![], targets, line_no))
}

fn split_head_and_targets(line: &str) -> (&str, &str) {
    let mut depth: usize = 0;
    for (i, c) in line.char_indices() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => depth = depth.saturating_sub(1),
            ' ' | '\t' if depth == 0 => return (line[..i].trim(), line[i + 1..].trim()),
            _ => {}
        }
    }
    (line.trim(), "")
}

/// Split `head` (e.g. `S[T](0.5)`) into `(name, Option<tags_src>, Option<args_src>)`.
/// Both `tags_src` and `args_src` are returned without their delimiters.
fn parse_head<'a>(
    head: &'a str,
    _line_no: usize,
) -> Result<(&'a str, Option<&'a str>, Option<&'a str>), ParseError> {
    // Tasks 5–6 fill in tag and arg parsing. For Task 4 we only need to support
    // a bare instruction name with no `[…]` and no `(…)`.
    Ok((head.trim(), None, None))
}

fn build_instruction(
    entry: TableEntry,
    _name_str: &str,
    _tags: Vec<crate::parser::ast::Tag>,
    args: Vec<f64>,
    targets: Vec<usize>,
    line: usize,
) -> RawInstruction {
    match entry {
        TableEntry::Gate { name, .. } => RawInstruction::Gate {
            name,
            tags: vec![],
            args,
            targets,
            line,
        },
        TableEntry::Noise { name, .. } => RawInstruction::Noise {
            name,
            tags: vec![],
            args,
            targets,
            line,
        },
        TableEntry::Measure { name, .. } => RawInstruction::Measure {
            name,
            tags: vec![],
            args,
            targets,
            line,
        },
        TableEntry::Annotation { kind, .. } => RawInstruction::Annotation {
            kind,
            args,
            targets,
            line,
        },
    }
}

// Suppress dead_code warnings for entries used in later tasks.
#[allow(dead_code)]
const _: (TargetArity, ArgCount, NoiseName, MeasureName, AnnotationKind, GateName) = (
    TargetArity::Any,
    ArgCount::None,
    NoiseName::Depolarize1,
    MeasureName::M,
    AnnotationKind::Detector,
    GateName::H,
);
