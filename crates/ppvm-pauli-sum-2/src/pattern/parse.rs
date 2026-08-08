// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::fmt;
use std::iter::Peekable;
use std::str::Chars;

use ppvm_traits_2::Pauli;

use super::data::{Decorated, OpPattern, PauliPattern, SiteSet};

/// A typed Pauli-pattern parse failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PatternParseError {
    ExpectedAtom,
    InvalidAlternation,
    EmptyAlternation,
    ExpectedDecoration,
    ExpectedCount,
    InvalidCount,
    InvalidCountCharacter(char),
}

impl fmt::Display for PatternParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ExpectedAtom => f.write_str("expected X, Y, Z, '_', or '['"),
            Self::InvalidAlternation => f.write_str("expected X, Y, Z, or ']' in alternation"),
            Self::EmptyAlternation => f.write_str("empty alternation"),
            Self::ExpectedDecoration => {
                f.write_str("expected '*', '{count}', or an absolute position")
            }
            Self::ExpectedCount => f.write_str("expected a decimal count"),
            Self::InvalidCount => f.write_str("pattern number does not fit usize"),
            Self::InvalidCountCharacter(ch) => {
                write!(f, "expected a digit or '}}' after '{{', found '{ch}'")
            }
        }
    }
}

impl std::error::Error for PatternParseError {}

impl PauliPattern {
    /// Parse the original Pauli-pattern grammar.
    pub fn parse(input: impl AsRef<str>) -> Result<Self, PatternParseError> {
        let mut chars = input.as_ref().trim().chars().peekable();
        let mut patterns = Vec::new();
        while chars.peek().is_some() {
            let mut op = parse_atom(&mut chars)?;
            if chars.next_if_eq(&'?').is_some() {
                op = op.optional();
            }
            match chars.peek().copied() {
                Some('*') => {
                    chars.next();
                    patterns.push(Decorated::Star(op));
                }
                Some('{') => {
                    chars.next();
                    patterns.push(Decorated::Repeat(op, parse_count(&mut chars, true)?));
                }
                Some(ch) if ch.is_ascii_digit() => {
                    patterns.push(Decorated::Position(op, parse_count(&mut chars, false)?));
                }
                _ => return Err(PatternParseError::ExpectedDecoration),
            }
        }
        Ok(Self(patterns))
    }
}

fn parse_atom(chars: &mut Peekable<Chars<'_>>) -> Result<OpPattern, PatternParseError> {
    match chars.next() {
        Some('X') => Ok(OpPattern::Single(Pauli::X)),
        Some('Y') => Ok(OpPattern::Single(Pauli::Y)),
        Some('Z') => Ok(OpPattern::Single(Pauli::Z)),
        Some('_') => Ok(OpPattern::AnyPauliOrIdentity),
        Some('[') => parse_alternation(chars),
        _ => Err(PatternParseError::ExpectedAtom),
    }
}

fn parse_alternation(chars: &mut Peekable<Chars<'_>>) -> Result<OpPattern, PatternParseError> {
    let mut set = SiteSet(0);
    while let Some(ch) = chars.peek().copied() {
        match ch {
            ']' => {
                chars.next();
                break;
            }
            'X' => {
                chars.next();
                set = set.union(SiteSet::X);
            }
            'Y' => {
                chars.next();
                set = set.union(SiteSet::Y);
            }
            'Z' => {
                chars.next();
                set = set.union(SiteSet::Z);
            }
            _ => return Err(PatternParseError::InvalidAlternation),
        }
    }
    op_from_nonidentity_set(set)
}

fn op_from_nonidentity_set(set: SiteSet) -> Result<OpPattern, PatternParseError> {
    match set.0 {
        0 => Err(PatternParseError::EmptyAlternation),
        bits if bits == SiteSet::X.0 => Ok(OpPattern::Single(Pauli::X)),
        bits if bits == SiteSet::Z.0 => Ok(OpPattern::Single(Pauli::Z)),
        bits if bits == SiteSet::Y.0 => Ok(OpPattern::Single(Pauli::Y)),
        bits if bits == (SiteSet::X.0 | SiteSet::Z.0) => Ok(OpPattern::Double(Pauli::X, Pauli::Z)),
        bits if bits == (SiteSet::X.0 | SiteSet::Y.0) => Ok(OpPattern::Double(Pauli::X, Pauli::Y)),
        bits if bits == (SiteSet::Z.0 | SiteSet::Y.0) => Ok(OpPattern::Double(Pauli::Z, Pauli::Y)),
        _ => Ok(OpPattern::AnyNonIdentity),
    }
}

fn parse_count(chars: &mut Peekable<Chars<'_>>, braced: bool) -> Result<usize, PatternParseError> {
    let mut number = String::new();
    while let Some(ch) = chars.peek().copied() {
        if ch.is_ascii_digit() {
            number.push(ch);
            chars.next();
        } else if braced && ch == '}' {
            chars.next();
            break;
        } else if braced {
            return Err(PatternParseError::InvalidCountCharacter(ch));
        } else {
            break;
        }
    }
    if number.is_empty() {
        return Err(PatternParseError::ExpectedCount);
    }
    number.parse().map_err(|_| PatternParseError::InvalidCount)
}
