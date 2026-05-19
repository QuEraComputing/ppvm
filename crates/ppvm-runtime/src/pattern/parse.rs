// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use anyhow::Result;
use std::collections::BTreeSet;

use super::data::{Decorated, NotIdentity, OpPattern, PauliPattern};

impl PauliPattern {
    /// Parse a pattern from its textual form.
    ///
    /// The grammar supports literal Paulis (`X`, `Y`, `Z`), the
    /// don't-care symbol `_`, alternation in brackets (`[XY]`),
    /// optional-identity suffix (`X?`), star repetition (`X*`),
    /// counted repetition (`X{3}`), and absolute-position anchoring
    /// (`X3` matches an `X` at qubit 3).
    ///
    /// # Examples
    ///
    /// ```
    /// use ppvm_runtime::pattern::{Contains, PauliPattern};
    /// use ppvm_runtime::word::PauliWord;
    ///
    /// let pat = PauliPattern::parse("X0Y1Z2").unwrap();
    /// let w: PauliWord<u64> = "XYZ".into();
    /// assert!(pat.contains(&w));
    ///
    /// let alt = PauliPattern::parse("[XY]0Y1Z2").unwrap();
    /// let w2: PauliWord<u64> = "YYZ".into();
    /// assert!(alt.contains(&w2));
    /// ```
    pub fn parse(input: impl AsRef<str>) -> Result<Self> {
        let mut patterns = Vec::new();
        let mut chars = input.as_ref().trim().chars().peekable();
        while let Some(c) = chars.next() {
            let op = match c {
                '[' => {
                    let mut alter = BTreeSet::new();
                    while let Some(ch) = chars.peek() {
                        let op = match ch {
                            ']' => {
                                chars.next(); // consume ']'
                                break;
                            }
                            'X' => {
                                chars.next();
                                NotIdentity::X
                            }
                            'Y' => {
                                chars.next();
                                NotIdentity::Y
                            }
                            'Z' => {
                                chars.next();
                                NotIdentity::Z
                            }
                            _ => Err(anyhow::anyhow!("Expected X, Y, or Z in []"))?,
                        };
                        alter.insert(op);
                    } // while

                    let mut alter_iter = alter.iter();
                    match alter.len() {
                        1 => OpPattern::Single(*alter_iter.next().unwrap()),
                        2 => OpPattern::Double(
                            *alter_iter.next().unwrap(),
                            *alter_iter.next().unwrap(),
                        ),
                        3 => OpPattern::AnyNonIdentity,
                        _ => Err(anyhow::anyhow!("Too many Pauli characters in pattern"))?,
                    }
                }
                '_' => OpPattern::AnyPauliOrIdentity,
                'X' => OpPattern::Single(NotIdentity::X),
                'Y' => OpPattern::Single(NotIdentity::Y),
                'Z' => OpPattern::Single(NotIdentity::Z),
                _ => Err(anyhow::anyhow!("Expected X, Y, or Z in Pauli pattern"))?,
            }; // op

            let op = match chars.peek() {
                Some('?') => {
                    chars.next();
                    op.add_identity()
                }
                _ => op,
            };

            match chars.peek() {
                Some('*') => {
                    chars.next(); // consume '*'
                    patterns.push(Decorated::Star(op));
                }
                Some('{') => {
                    chars.next(); // consume '{'
                    let mut num = String::new();
                    while let Some(ch) = chars.peek() {
                        if ch.is_ascii_digit() {
                            num.push(*ch);
                            chars.next();
                        } else if *ch == '}' {
                            chars.next(); // consume '}'
                            break;
                        } else {
                            Err(anyhow::anyhow!("Expected digit or '}}' after '{{'"))?;
                        }
                    }
                    let count = num.parse::<usize>()?;
                    patterns.push(Decorated::Repeat(op, count));
                }
                Some(x) if x.is_ascii_digit() => {
                    let mut num = String::new();
                    while let Some(ch) = chars.peek() {
                        if ch.is_ascii_digit() {
                            num.push(*ch);
                            chars.next();
                        } else {
                            break;
                        }
                    }
                    let position = num.parse::<usize>()?;
                    patterns.push(Decorated::Position(op, position));
                }
                _ => Err(anyhow::anyhow!("Expected '*' or digit after Pauli pattern"))?,
            }
        }
        Ok(PauliPattern(patterns))
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = &Decorated> {
        self.0.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use Decorated::*;
    use NotIdentity::*;
    use OpPattern::*;

    #[test]
    fn test_parse_patterns() {
        let answer = PauliPattern::parse("X1Y2Z3").unwrap();
        let expect = PauliPattern(vec![
            Position(Single(X), 1),
            Position(Single(Y), 2),
            Position(Single(Z), 3),
        ]);
        assert_eq!(answer, expect);

        let answer = PauliPattern::parse("[XY]1Y2_3").unwrap();
        let expect = PauliPattern(vec![
            Position(Double(X, Y), 1),
            Position(Single(Y), 2),
            Position(AnyPauliOrIdentity, 3),
        ]);
        assert_eq!(answer, expect);

        let answer = PauliPattern::parse("X?1Y2Z3").unwrap();
        let expect = PauliPattern(vec![
            Position(SingleOrIdentity(X), 1),
            Position(Single(Y), 2),
            Position(Single(Z), 3),
        ]);
        assert_eq!(answer, expect);

        let answer = PauliPattern::parse("[XY]*Z5").unwrap();
        let expect = PauliPattern(vec![Star(Double(X, Y)), Position(Single(Z), 5)]);
        assert_eq!(answer, expect);

        let answer = PauliPattern::parse("[XY]{3}Z5").unwrap();
        let expect = PauliPattern(vec![Repeat(Double(X, Y), 3), Position(Single(Z), 5)]);
        assert_eq!(answer, expect);
    }

    #[test]
    #[cfg(feature = "bincode")]
    fn test_bincode() {
        use bincode;
        let pat = PauliPattern::parse("[XY]{2}Z3").unwrap();
        let mut encoded = [0u8; 100];
        bincode::encode_into_slice(&pat, &mut encoded[..], bincode::config::standard()).unwrap();
        let (decoded, _) =
            bincode::decode_from_slice(&encoded, bincode::config::standard()).unwrap();
        assert_eq!(pat, decoded);
    }
}
