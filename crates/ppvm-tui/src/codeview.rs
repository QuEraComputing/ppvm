// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! A minimal scrollable list of lines with one optionally-highlighted row.
//! Backs the Program panel (highlight = the program counter) and the REPL log
//! (no highlight). Deliberately mirrors the shape of stellarscope's `CodeView`
//! so the two stay interchangeable, but has no external dependencies.

/// A list of displayable lines plus an optional cursor (the highlighted row).
#[derive(Debug, Clone, Default)]
pub struct CodeView<T> {
    lines: Vec<T>,
    cursor: Option<usize>,
}

impl<T> CodeView<T> {
    pub fn new() -> Self {
        Self {
            lines: Vec::new(),
            cursor: None,
        }
    }

    /// Drop all lines and clear the cursor.
    pub fn clear(&mut self) {
        self.lines.clear();
        self.cursor = None;
    }

    /// Append one line.
    pub fn push(&mut self, line: T) {
        self.lines.push(line);
    }

    /// Highlight row `idx` (or none).
    pub fn set_cursor(&mut self, idx: Option<usize>) {
        self.cursor = idx;
    }

    pub fn cursor(&self) -> Option<usize> {
        self.cursor
    }

    pub fn lines(&self) -> &[T] {
        &self.lines
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn push_and_cursor_round_trip() {
        let mut cv: CodeView<String> = CodeView::new();
        assert!(cv.lines().is_empty());
        assert_eq!(cv.cursor(), None);

        cv.push("0000: h".to_string());
        cv.push("0001: measure".to_string());
        cv.set_cursor(Some(1));

        assert_eq!(cv.lines().len(), 2);
        assert_eq!(cv.cursor(), Some(1));

        cv.clear();
        assert!(cv.lines().is_empty());
        assert_eq!(cv.cursor(), None);
    }
}
