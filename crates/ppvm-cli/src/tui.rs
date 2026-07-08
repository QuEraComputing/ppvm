// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Terminal ownership for the ppvm TUI: raw mode + alternate screen behind an
//! RAII guard (restored even on panic), plus a blocking event loop that drives
//! the terminal-agnostic `ppvm_tui::AppState`.

use std::io;

use crossterm::event::{self, Event};
use crossterm::execute;
use crossterm::terminal::{
    EnterAlternateScreen, LeaveAlternateScreen, disable_raw_mode, enable_raw_mode,
};
use eyre::Result;
use ppvm_tui::AppState;
use ratatui::Terminal;
use ratatui::backend::CrosstermBackend;

/// Restores the terminal on drop — including when the app panics mid-loop.
struct TerminalGuard;

impl Drop for TerminalGuard {
    fn drop(&mut self) {
        let _ = disable_raw_mode();
        let _ = execute!(io::stdout(), LeaveAlternateScreen);
    }
}

/// Launch the TUI. With `file`, open it loaded and paused at pc 0; without,
/// start an empty REPL session.
pub fn run(file: Option<&str>) -> Result<()> {
    let mut app = match file {
        Some(path) => AppState::from_file(path)?,
        None => AppState::new(),
    };

    // Guard immediately after raw mode is on, so any later setup error
    // (EnterAlternateScreen, Terminal::new) still restores the terminal.
    enable_raw_mode()?;
    let _guard = TerminalGuard;
    execute!(io::stdout(), EnterAlternateScreen)?;

    let mut terminal = Terminal::new(CrosstermBackend::new(io::stdout()))?;
    terminal.clear()?;

    while !app.should_exit {
        terminal.draw(|frame| app.render(frame))?;
        if let Event::Key(key) = event::read()? {
            app.handle_key(key);
        }
    }
    Ok(())
}
