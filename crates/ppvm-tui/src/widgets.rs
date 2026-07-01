// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! ratatui `Widget` components for the ppvm TUI. Each borrows `&AppState` and
//! only reads it, so a host app can render any of them into its own layout.

use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::text::{Line, Text};
use ratatui::widgets::{Block, Paragraph, Widget, Wrap};

use crate::app::AppState;

/// The left panel: the loaded program's listing (with a `▶` at the pc) or, in a
/// REPL session, the command/result log.
pub struct ProgramView<'a>(pub &'a AppState);

impl Widget for ProgramView<'_> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        let (title, view) = self.0.active_listing();
        let mut text = Text::default();
        for (i, line) in view.lines().iter().enumerate() {
            let marked = if view.cursor() == Some(i) {
                format!("▶ {line}")
            } else {
                format!("  {line}")
            };
            text.push_line(Line::from(marked));
        }
        Paragraph::new(text)
            .block(Block::bordered().title(title))
            .render(area, buf);
    }
}

/// The right panel: the tableau state (`PPVM::state_string`).
pub struct StateView<'a>(pub &'a AppState);

impl Widget for StateView<'_> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        Paragraph::new(self.0.state_text())
            .block(Block::bordered().title("State"))
            .wrap(Wrap { trim: false })
            .render(area, buf);
    }
}

/// The measurement-record band.
pub struct RecordView<'a>(pub &'a AppState);

impl Widget for RecordView<'_> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        Paragraph::new(self.0.measurement_bits())
            .block(Block::bordered().title("Measurement record"))
            .render(area, buf);
    }
}

/// The footer: prompt + input, then the hint and status line.
pub struct CommandLine<'a>(pub &'a AppState);

impl CommandLine<'_> {
    /// The command prompt. Also used to offset the terminal cursor past it.
    pub const PROMPT: &'static str = "ppvm> ";
}

impl Widget for CommandLine<'_> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        let text = Text::from(vec![
            Line::from(format!("{}{}", Self::PROMPT, self.0.input())),
            Line::from(format!("{}    {}", self.0.hint(), self.0.status())),
        ]);
        Paragraph::new(text).render(area, buf);
    }
}

#[cfg(test)]
mod tests {
    use crate::AppState;
    use ratatui::Terminal;
    use ratatui::backend::TestBackend;

    const BP_PROGRAM: &str = "device circuit.n_qubits 1;\n\
                              fn @main() { breakpoint\n const.u64 0\n circuit.measure\n ret }\n";

    #[test]
    fn renders_all_panels_without_panic() {
        let mut app = AppState::new();
        app.load_source(BP_PROGRAM).unwrap();

        let backend = TestBackend::new(100, 30);
        let mut terminal = Terminal::new(backend).unwrap();
        terminal.draw(|f| app.render(f)).unwrap();

        let content: String = terminal
            .backend()
            .buffer()
            .content
            .iter()
            .map(|c| c.symbol())
            .collect();
        assert!(content.contains("Program"), "missing Program panel");
        assert!(content.contains("State"), "missing State panel");
        assert!(
            content.contains("Measurement record"),
            "missing record panel"
        );
        assert!(content.contains("ppvm>"), "missing command prompt");
    }
}
