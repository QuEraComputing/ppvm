// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! A readline-style single-line editor with command history. Owns the command
//! buffer, edit cursor, and history; knows nothing about the PPVM or the
//! debugger. The host handles submit (Enter) and quit keys, then delegates the
//! remaining editing keys here via [`LineEditor::handle_key`].

use crossterm::event::{KeyCode, KeyEvent};

/// The command line's editable buffer plus its command history.
#[derive(Debug, Default)]
pub struct LineEditor {
    /// The command-line buffer.
    input: String,
    /// Char index of the edit cursor within `input` (0..=input char count).
    cursor: usize,
    /// Submitted command lines, oldest first (for Up/Down recall).
    history: Vec<String>,
    /// Position within `history` while recalling; `None` = editing the live line.
    history_pos: Option<usize>,
    /// The live line stashed when entering history, restored on the way back.
    draft: String,
}

impl LineEditor {
    pub fn new() -> Self {
        Self::default()
    }

    // ─── read-only accessors ─────────────────────────────────────────────

    pub fn input(&self) -> &str {
        &self.input
    }

    /// Char index of the edit cursor within the input line. Used to place the
    /// terminal cursor.
    pub fn cursor(&self) -> usize {
        self.cursor
    }

    pub fn is_empty(&self) -> bool {
        self.input.is_empty()
    }

    /// Clear the input line and reset the cursor.
    pub fn clear(&mut self) {
        self.input.clear();
        self.cursor = 0;
    }

    // ─── key handling ────────────────────────────────────────────────────

    /// Apply one editing key (text, cursor movement, or history recall).
    /// Returns whether it was consumed. Submit and quit keys are the host's
    /// responsibility and are not handled here.
    pub fn handle_key(&mut self, key: KeyEvent) -> bool {
        match key.code {
            KeyCode::Char(c) => {
                self.insert_char(c);
                true
            }
            KeyCode::Backspace => {
                self.backspace();
                true
            }
            KeyCode::Delete => {
                self.delete();
                true
            }
            KeyCode::Left => {
                self.cursor = self.cursor.saturating_sub(1);
                true
            }
            KeyCode::Right => {
                self.cursor = (self.cursor + 1).min(self.input_len());
                true
            }
            KeyCode::Home => {
                self.cursor = 0;
                true
            }
            KeyCode::End => {
                self.cursor = self.input_len();
                true
            }
            KeyCode::Up => {
                self.history_prev();
                true
            }
            KeyCode::Down => {
                self.history_next();
                true
            }
            _ => false,
        }
    }

    /// Take the current line: record it in history and reset editor state,
    /// returning the (untrimmed) line for the host to dispatch.
    pub fn submit(&mut self) -> String {
        let line = std::mem::take(&mut self.input);
        self.cursor = 0;
        self.history_pos = None;
        self.draft.clear();
        let trimmed = line.trim();
        // Skip blanks and consecutive duplicates, matching a shell's history.
        if !trimmed.is_empty() && self.history.last().map(String::as_str) != Some(trimmed) {
            self.history.push(trimmed.to_string());
        }
        line
    }

    // ─── line editing ────────────────────────────────────────────────────

    /// Char count of the current input line.
    fn input_len(&self) -> usize {
        self.input.chars().count()
    }

    /// Byte offset of char index `i` within `input` (or the end of the string).
    fn byte_index(&self, i: usize) -> usize {
        self.input
            .char_indices()
            .nth(i)
            .map(|(b, _)| b)
            .unwrap_or(self.input.len())
    }

    /// Insert `c` at the cursor and step past it.
    fn insert_char(&mut self, c: char) {
        let b = self.byte_index(self.cursor);
        self.input.insert(b, c);
        self.cursor += 1;
    }

    /// Delete the char before the cursor (Backspace).
    fn backspace(&mut self) {
        if self.cursor > 0 {
            let b = self.byte_index(self.cursor - 1);
            self.input.remove(b);
            self.cursor -= 1;
        }
    }

    /// Delete the char at the cursor (Delete).
    fn delete(&mut self) {
        if self.cursor < self.input_len() {
            let b = self.byte_index(self.cursor);
            self.input.remove(b);
        }
    }

    /// Replace the input line, moving the cursor to its end.
    fn set_input(&mut self, line: String) {
        self.cursor = line.chars().count();
        self.input = line;
    }

    // ─── command history ─────────────────────────────────────────────────

    /// Recall an older history entry (Up), stashing the live line on first entry.
    fn history_prev(&mut self) {
        if self.history.is_empty() {
            return;
        }
        let i = match self.history_pos {
            None => {
                self.draft = std::mem::take(&mut self.input);
                self.history.len() - 1
            }
            Some(0) => return, // already at the oldest
            Some(i) => i - 1,
        };
        self.history_pos = Some(i);
        self.set_input(self.history[i].clone());
    }

    /// Move toward newer entries (Down); past the newest, restore the stashed
    /// live line.
    fn history_next(&mut self) {
        match self.history_pos {
            None => {}
            Some(i) if i + 1 < self.history.len() => {
                self.history_pos = Some(i + 1);
                self.set_input(self.history[i + 1].clone());
            }
            Some(_) => {
                self.history_pos = None;
                let draft = std::mem::take(&mut self.draft);
                self.set_input(draft);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crossterm::event::{KeyCode, KeyModifiers};

    fn key(code: KeyCode) -> KeyEvent {
        KeyEvent::new(code, KeyModifiers::NONE)
    }

    /// Type each char, then submit (recording it in history).
    fn type_line(ed: &mut LineEditor, line: &str) {
        for c in line.chars() {
            ed.handle_key(key(KeyCode::Char(c)));
        }
        ed.submit();
    }

    #[test]
    fn typing_inserts_and_advances_the_cursor() {
        let mut ed = LineEditor::new();
        for c in "hz".chars() {
            ed.handle_key(key(KeyCode::Char(c)));
        }
        assert_eq!(ed.input(), "hz");
        assert_eq!(ed.cursor(), 2);
    }

    #[test]
    fn backspace_edits_the_buffer() {
        let mut ed = LineEditor::new();
        ed.handle_key(key(KeyCode::Char('h')));
        ed.handle_key(key(KeyCode::Char('i')));
        ed.handle_key(key(KeyCode::Backspace));
        assert_eq!(ed.input(), "h");
    }

    #[test]
    fn left_right_home_end_move_the_cursor() {
        let mut ed = LineEditor::new();
        for c in "abc".chars() {
            ed.handle_key(key(KeyCode::Char(c)));
        }
        assert_eq!(ed.cursor(), 3);
        ed.handle_key(key(KeyCode::Left));
        assert_eq!(ed.cursor(), 2);
        ed.handle_key(key(KeyCode::Home));
        assert_eq!(ed.cursor(), 0);
        ed.handle_key(key(KeyCode::Left)); // saturates at 0
        assert_eq!(ed.cursor(), 0);
        ed.handle_key(key(KeyCode::End));
        assert_eq!(ed.cursor(), 3);
        ed.handle_key(key(KeyCode::Right)); // saturates at len
        assert_eq!(ed.cursor(), 3);
    }

    #[test]
    fn insert_and_delete_at_the_cursor() {
        let mut ed = LineEditor::new();
        for c in "ac".chars() {
            ed.handle_key(key(KeyCode::Char(c)));
        }
        ed.handle_key(key(KeyCode::Left)); // cursor between 'a' and 'c'
        ed.handle_key(key(KeyCode::Char('b')));
        assert_eq!(ed.input(), "abc");
        assert_eq!(ed.cursor(), 2);
        ed.handle_key(key(KeyCode::Backspace)); // removes 'b' before cursor
        assert_eq!(ed.input(), "ac");
        assert_eq!(ed.cursor(), 1);
        ed.handle_key(key(KeyCode::Delete)); // removes 'c' at cursor
        assert_eq!(ed.input(), "a");
        assert_eq!(ed.cursor(), 1);
    }

    #[test]
    fn submit_records_history_and_clears_the_line() {
        let mut ed = LineEditor::new();
        for c in "x 0".chars() {
            ed.handle_key(key(KeyCode::Char(c)));
        }
        let line = ed.submit();
        assert_eq!(line, "x 0", "submit returns the entered line");
        assert!(ed.is_empty(), "buffer clears on submit");
        assert_eq!(ed.cursor(), 0);
    }

    #[test]
    fn up_arrow_recalls_previous_commands() {
        let mut ed = LineEditor::new();
        type_line(&mut ed, "device 1");
        type_line(&mut ed, "x 0");
        ed.handle_key(key(KeyCode::Up)); // newest first
        assert_eq!(ed.input(), "x 0");
        assert_eq!(ed.cursor(), 3, "cursor lands at end of the recalled line");
        ed.handle_key(key(KeyCode::Up));
        assert_eq!(ed.input(), "device 1");
        ed.handle_key(key(KeyCode::Up)); // already oldest, stays put
        assert_eq!(ed.input(), "device 1");
    }

    #[test]
    fn down_arrow_returns_toward_the_live_line() {
        let mut ed = LineEditor::new();
        type_line(&mut ed, "device 1");
        type_line(&mut ed, "x 0");
        for c in "meas".chars() {
            ed.handle_key(key(KeyCode::Char(c)));
        }
        ed.handle_key(key(KeyCode::Up)); // stash "meas", show "x 0"
        ed.handle_key(key(KeyCode::Up)); // "device 1"
        ed.handle_key(key(KeyCode::Down)); // "x 0"
        assert_eq!(ed.input(), "x 0");
        ed.handle_key(key(KeyCode::Down)); // past newest -> restore draft
        assert_eq!(ed.input(), "meas");
        assert_eq!(ed.cursor(), 4);
    }

    #[test]
    fn history_skips_consecutive_duplicates() {
        let mut ed = LineEditor::new();
        type_line(&mut ed, "device 1");
        type_line(&mut ed, "device 1"); // same command twice -> one entry
        ed.handle_key(key(KeyCode::Up));
        assert_eq!(ed.input(), "device 1");
        ed.handle_key(key(KeyCode::Up)); // only one entry -> stays
        assert_eq!(ed.input(), "device 1");
    }

    #[test]
    fn clear_empties_the_line() {
        let mut ed = LineEditor::new();
        for c in "abc".chars() {
            ed.handle_key(key(KeyCode::Char(c)));
        }
        ed.clear();
        assert!(ed.is_empty());
        assert_eq!(ed.cursor(), 0);
    }

    #[test]
    fn non_editing_keys_are_not_consumed() {
        let mut ed = LineEditor::new();
        assert!(!ed.handle_key(key(KeyCode::Enter)));
        assert!(!ed.handle_key(key(KeyCode::Esc)));
    }
}
