pub mod ast;
pub mod table;

use ast::{
    AnnotationKind, ArgCount, GateName, MeasureName, NoiseName,
    ParseError, Program, RawInstruction, TargetArity,
};
use table::{TableEntry, lookup};

/// Parse `<coeff>*pi`, `pi`, or a plain f64.
fn parse_pi_expr(s: &str, line_no: usize) -> Result<f64, ParseError> {
    let s = s.trim();
    if s == "pi" {
        return Ok(std::f64::consts::PI);
    }
    if let Some(coeff) = s.strip_suffix("*pi") {
        return coeff
            .trim()
            .parse::<f64>()
            .map(|c| c * std::f64::consts::PI)
            .map_err(|_| ParseError::Syntax {
                line: line_no,
                col: 1,
                message: format!("invalid pi-expression {s:?}"),
            });
    }
    s.parse::<f64>().map_err(|_| ParseError::Syntax {
        line: line_no,
        col: 1,
        message: format!("invalid number {s:?}"),
    })
}

/// Split `s` by commas that are not inside parentheses or brackets.
fn split_commas_shallow(s: &str) -> Vec<&str> {
    let mut depth = 0usize;
    let mut start = 0;
    let mut result = Vec::new();
    for (i, c) in s.char_indices() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => depth = depth.saturating_sub(1),
            ',' if depth == 0 => {
                result.push(s[start..i].trim());
                start = i + 1;
            }
            _ => {}
        }
    }
    let last = s[start..].trim();
    if !last.is_empty() {
        result.push(last);
    }
    result
}

fn parse_tags(tags_src: &str, line_no: usize) -> Result<Vec<crate::parser::ast::Tag>, ParseError> {
    use crate::parser::ast::{Tag, TagParam};
    let mut out = Vec::new();
    for tag in split_commas_shallow(tags_src) {
        let tag = tag.trim();
        if tag.is_empty() {
            continue;
        }
        if let Some(paren_start) = tag.find('(') {
            let name = tag[..paren_start].trim().to_string();
            let inner = tag[paren_start + 1..]
                .strip_suffix(')')
                .ok_or(ParseError::Syntax {
                    line: line_no,
                    col: paren_start + 1,
                    message: "unclosed '(' in tag".into(),
                })?;
            let mut params = Vec::new();
            for raw in split_commas_shallow(inner) {
                let raw = raw.trim();
                if raw.is_empty() {
                    continue;
                }
                if let Some(eq_idx) = raw.find('=') {
                    let key = raw[..eq_idx].trim().to_string();
                    let value = parse_pi_expr(raw[eq_idx + 1..].trim(), line_no)?;
                    params.push(TagParam::Named { key, value });
                } else {
                    let value = parse_pi_expr(raw, line_no)?;
                    params.push(TagParam::Positional(value));
                }
            }
            out.push(Tag { name, params });
        } else {
            out.push(Tag {
                name: tag.to_string(),
                params: vec![],
            });
        }
    }
    Ok(out)
}

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
    let (head, targets_part) = split_head_and_targets(line);
    let (name_str, tags_src, args_src) = parse_head(head, line_no)?;

    let entry = lookup(name_str).ok_or(ParseError::UnknownInstruction {
        name: name_str.to_string(),
        line: line_no,
    })?;

    let tags = match tags_src {
        Some(s) => parse_tags(s, line_no)?,
        None => Vec::new(),
    };

    let args: Vec<f64> = match args_src {
        Some(s) if !s.trim().is_empty() => split_commas_shallow(s)
            .into_iter()
            .map(|p| parse_pi_expr(p, line_no))
            .collect::<Result<_, _>>()?,
        _ => Vec::new(),
    };

    let targets: Vec<usize> = targets_part
        .split_whitespace()
        .map(|t| t.parse::<usize>().map_err(|_| ParseError::Syntax {
            line: line_no,
            col: 1,
            message: format!("invalid target {t:?}"),
        }))
        .collect::<Result<_, _>>()?;

    Ok(build_instruction(entry, name_str, tags, args, targets, line_no))
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
    line_no: usize,
) -> Result<(&'a str, Option<&'a str>, Option<&'a str>), ParseError> {
    if let Some(bracket_start) = head.find('[') {
        let name = head[..bracket_start].trim();
        let after_open = &head[bracket_start + 1..];
        let bracket_end = after_open.find(']').ok_or(ParseError::Syntax {
            line: line_no,
            col: bracket_start + 1,
            message: "unclosed '['".into(),
        })?;
        let tags_src = &after_open[..bracket_end];
        let after_bracket = after_open[bracket_end + 1..].trim();
        let args_src = match after_bracket.strip_prefix('(') {
            Some(rest) => Some(rest.strip_suffix(')').ok_or(ParseError::Syntax {
                line: line_no,
                col: bracket_start + bracket_end + 2,
                message: "unclosed '('".into(),
            })?),
            None => None,
        };
        Ok((name, Some(tags_src), args_src))
    } else if let Some(paren_start) = head.find('(') {
        let name = head[..paren_start].trim();
        let inner = head[paren_start + 1..]
            .strip_suffix(')')
            .ok_or(ParseError::Syntax {
                line: line_no,
                col: paren_start + 1,
                message: "unclosed '('".into(),
            })?;
        Ok((name, None, Some(inner)))
    } else {
        Ok((head.trim(), None, None))
    }
}

fn build_instruction(
    entry: TableEntry,
    _name_str: &str,
    tags: Vec<crate::parser::ast::Tag>,
    args: Vec<f64>,
    targets: Vec<usize>,
    line: usize,
) -> RawInstruction {
    match entry {
        TableEntry::Gate { name, .. } => RawInstruction::Gate { name, tags, args, targets, line },
        TableEntry::Noise { name, .. } => RawInstruction::Noise { name, tags, args, targets, line },
        TableEntry::Measure { name, .. } => RawInstruction::Measure { name, tags, args, targets, line },
        TableEntry::Annotation { kind, .. } => RawInstruction::Annotation { kind, args, targets, line },
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
