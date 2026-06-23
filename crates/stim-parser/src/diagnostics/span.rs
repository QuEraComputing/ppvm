// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Byte spans and the line/column map shared by every diagnostic.

/// Half-open byte range `[start, end)` into the source string.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}

impl Span {
    pub fn new(start: usize, end: usize) -> Self {
        Span { start, end }
    }

    /// 1-indexed `(line, col)` of the span start.
    pub fn line_col(&self, line_map: &LineMap) -> (usize, usize) {
        line_map.line_col(self.start)
    }

    /// 1-indexed line of the span start.
    pub fn line(&self, line_map: &LineMap) -> usize {
        line_map.line_of(self.start)
    }
}

impl From<chumsky::span::SimpleSpan<usize>> for Span {
    fn from(s: chumsky::span::SimpleSpan<usize>) -> Self {
        Span::new(s.start, s.end)
    }
}

/// Maps byte offsets in source to 1-indexed line/column positions.
pub struct LineMap {
    /// `starts[i]` = byte offset of the start of line (i+1).
    starts: Vec<usize>,
}

impl std::fmt::Debug for LineMap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LineMap")
            .field("lines", &self.starts.len())
            .finish()
    }
}

impl LineMap {
    /// Build a `LineMap` for `src`.
    pub fn new(src: &str) -> Self {
        let mut starts = vec![0];
        for (i, b) in src.bytes().enumerate() {
            if b == b'\n' {
                starts.push(i + 1);
            }
        }
        Self { starts }
    }

    /// 1-indexed line number for a byte offset.
    pub fn line_of(&self, byte_offset: usize) -> usize {
        match self.starts.binary_search(&byte_offset) {
            Ok(i) => i + 1,
            Err(i) => i, // i is the insertion index; start of line `i` is at starts[i-1].
        }
    }

    /// 1-indexed `(line, col)` for a byte offset.
    pub fn line_col(&self, byte_offset: usize) -> (usize, usize) {
        let line = self.line_of(byte_offset);
        let line_start = self.starts[line - 1];
        let col = byte_offset - line_start + 1;
        (line, col)
    }

    /// Byte offset of the start of line `(line_idx + 1)`. `None` for
    /// out-of-range indices.
    pub fn starts_at(&self, line_idx: usize) -> Option<usize> {
        self.starts.get(line_idx).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_col_at_line_start() {
        let m = LineMap::new("abc\ndef\nghi");
        assert_eq!(m.line_col(0), (1, 1));
        assert_eq!(m.line_col(4), (2, 1));
        assert_eq!(m.line_col(8), (3, 1));
    }

    #[test]
    fn line_col_mid_line() {
        let m = LineMap::new("abc\ndef\nghi");
        assert_eq!(m.line_col(2), (1, 3));
        assert_eq!(m.line_col(6), (2, 3));
    }

    #[test]
    fn span_resolves_against_line_map() {
        let m = LineMap::new("X 0\nH 0");
        let span = Span::new(4, 5);
        assert_eq!(span.line_col(&m), (2, 1));
        assert_eq!(span.line(&m), 2);
    }
}
