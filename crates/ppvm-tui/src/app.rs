// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `AppState` — the terminal-agnostic state of the ppvm TUI. Owns an optional
//! [`PPVM`] plus the command buffer, status line, program listing, and REPL
//! log. `dispatch` runs one command string; `handle_key` edits the buffer and
//! submits on Enter. Nothing here touches a terminal or runs a loop.

use crossterm::event::{KeyCode, KeyEvent, KeyEventKind, KeyModifiers};
use eyre::{Result, eyre};
use ppvm_vihaco::composite::{PPVM, StepOutcome};
use ppvm_vihaco::measurements::MeasurementResult;
use ppvm_vihaco::{CircuitInstruction, PPVMModule, compile_program, load_module_file};
use ratatui::Frame;
use ratatui::layout::{Constraint, Layout};

use crate::codeview::CodeView;
use crate::command::{Command, parse_command};
use crate::widgets::{CommandLine, ProgramView, RecordView, StateView};

/// Terminal-agnostic state for the ppvm TUI.
pub struct AppState {
    /// The live machine. `None` until `device N` or a program is loaded.
    machine: Option<PPVM>,
    /// The loaded module (kept for `:reset`).
    module: Option<PPVMModule>,
    /// Program instruction listing (populated when a program is loaded).
    program: CodeView<String>,
    /// REPL scrollback: entered commands and inline results.
    log: CodeView<String>,
    /// Qubit count of the current REPL device (for `:reset`).
    n_qubits: usize,
    /// True while a program is loaded (Program panel) vs a REPL session (Log).
    has_program: bool,
    /// True while the debugger is paused (at start or a breakpoint).
    paused: bool,
    /// True once the loaded program has run to Return/Halt.
    finished: bool,
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
    /// The status/error line.
    status: String,
    /// Set to leave the event loop.
    pub should_exit: bool,
}

impl Default for AppState {
    fn default() -> Self {
        Self::new()
    }
}

impl AppState {
    pub fn new() -> Self {
        Self {
            machine: None,
            module: None,
            program: CodeView::new(),
            log: CodeView::new(),
            n_qubits: 0,
            has_program: false,
            paused: false,
            finished: false,
            input: String::new(),
            cursor: 0,
            history: Vec::new(),
            history_pos: None,
            draft: String::new(),
            status: String::new(),
            should_exit: false,
        }
    }

    // ─── command dispatch ────────────────────────────────────────────────

    /// Run one command-line string. Command-level errors are non-fatal: they
    /// are written to the status line and the app keeps running.
    pub fn dispatch(&mut self, line: &str) {
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            self.log.push(format!("ppvm> {trimmed}"));
        }
        let result = match parse_command(line) {
            Ok(cmd) => self.run_command(cmd),
            Err(e) => Err(e),
        };
        if let Err(e) = result {
            self.set_status(format!("error: {e}"));
        }
    }

    fn run_command(&mut self, cmd: Command) -> Result<()> {
        match cmd {
            Command::Quit => {
                self.should_exit = true;
                Ok(())
            }
            Command::Device(n) => self.new_device(n),
            Command::Gate {
                inst,
                qubits,
                params,
            } => self.apply_gate(inst, &qubits, &params),
            Command::Step => self.step(),
            Command::Continue => self.cont(),
            Command::Reset => self.reset(),
            Command::Load(path) => self.load_file(&path),
        }
    }

    fn new_device(&mut self, n: usize) -> Result<()> {
        self.machine = Some(PPVM::with_qubits(n)?);
        self.n_qubits = n;
        self.has_program = false;
        self.paused = false;
        self.finished = false;
        self.program.clear();
        self.set_status(format!("fresh {n}-qubit device"));
        Ok(())
    }

    // ─── program loading ─────────────────────────────────────────────────

    /// Build an `AppState` with `path` loaded and paused at pc 0.
    pub fn from_file(path: &str) -> Result<Self> {
        let mut app = Self::new();
        app.load_file(path)?;
        Ok(app)
    }

    /// Compile `.sst` source and load it, paused at pc 0. (Test/embedding entry
    /// that avoids touching the filesystem.)
    pub fn load_source(&mut self, src: &str) -> Result<()> {
        let module = compile_program(src)?;
        self.load_module(module);
        self.set_status("loaded program");
        Ok(())
    }

    fn load_file(&mut self, path: &str) -> Result<()> {
        let module = load_module_file(path).map_err(|e| eyre!("failed to load {path}: {e}"))?;
        self.load_module(module);
        self.set_status(format!("loaded {path}"));
        Ok(())
    }

    /// Core loader: rebuild the machine from `module` and pause at pc 0.
    fn load_module(&mut self, module: PPVMModule) {
        let mut m = PPVM::default();
        // A fresh machine + load + init gives clean tableau/record state; these
        // only fail on malformed modules, which `compile_program` already
        // rejects, so surface as a status rather than unwinding the UI.
        if let Err(e) = m.load(&module).and_then(|()| m.init()) {
            self.set_status(format!("error: {e}"));
            return;
        }
        self.program.clear();
        for (i, inst) in module.code.iter().enumerate() {
            self.program.push(format!("{i:04}: {inst}"));
        }
        self.machine = Some(m);
        self.module = Some(module);
        self.has_program = true;
        self.paused = true;
        self.finished = false;
        self.refresh_cursor();
    }

    fn refresh_cursor(&mut self) {
        let pc = self.machine.as_ref().map(|m| m.current_pc() as usize);
        if let Some(pc) = pc {
            self.program.set_cursor(Some(pc));
        }
    }

    // ─── stepping ────────────────────────────────────────────────────────

    fn step(&mut self) -> Result<()> {
        if !self.has_program {
            self.set_status("nothing to step — load a program with :load");
            return Ok(());
        }
        if self.finished {
            self.set_status("program finished — :reset to run again");
            return Ok(());
        }
        let outcome = self.machine.as_mut().unwrap().step_once()?;
        self.apply_outcome(outcome);
        self.refresh_cursor();
        Ok(())
    }

    fn cont(&mut self) -> Result<()> {
        if !self.has_program {
            self.set_status("nothing to continue — load a program with :load");
            return Ok(());
        }
        while !self.finished {
            let outcome = self.machine.as_mut().unwrap().step_once()?;
            match outcome {
                StepOutcome::Continue => {}
                StepOutcome::Breakpoint => {
                    self.apply_outcome(outcome);
                    self.refresh_cursor();
                    return Ok(());
                }
                StepOutcome::Return | StepOutcome::Halt => {
                    self.apply_outcome(outcome);
                    break;
                }
            }
        }
        self.refresh_cursor();
        Ok(())
    }

    /// Fold a single step outcome into the app's paused/finished/status state.
    fn apply_outcome(&mut self, outcome: StepOutcome) {
        match outcome {
            StepOutcome::Continue => self.set_status(""),
            StepOutcome::Breakpoint => {
                self.paused = true;
                self.set_status("-- breakpoint hit --");
            }
            StepOutcome::Return | StepOutcome::Halt => {
                self.finished = true;
                self.set_status("program finished");
            }
        }
    }

    fn reset(&mut self) -> Result<()> {
        if let Some(module) = self.module.clone() {
            self.load_module(module);
            self.set_status("reset");
        } else if self.n_qubits > 0 {
            self.machine = Some(PPVM::with_qubits(self.n_qubits)?);
            self.set_status("reset device");
        } else {
            self.set_status("nothing to reset");
        }
        Ok(())
    }

    pub fn has_program(&self) -> bool {
        self.has_program
    }

    pub fn paused(&self) -> bool {
        self.paused
    }

    fn apply_gate(
        &mut self,
        inst: CircuitInstruction,
        qubits: &[usize],
        params: &[f64],
    ) -> Result<()> {
        let m = self
            .machine
            .as_mut()
            .ok_or_else(|| eyre!("no device — run `device N` or :load a file first"))?;
        let before = m.measurement_record().len();
        m.apply_circuit_instruction(inst, qubits, params)?;
        // Any new record entries are this gate's measurement outcomes.
        let new: Vec<MeasurementResult> = m.measurement_record()[before..].to_vec();
        if new.is_empty() {
            self.set_status("");
        } else {
            let bits = format_record(&new);
            self.log.push(format!("  => {bits}"));
            self.set_status(format!("=> {bits}"));
        }
        Ok(())
    }

    fn set_status(&mut self, s: impl Into<String>) {
        self.status = s.into();
    }

    // ─── key handling ────────────────────────────────────────────────────

    /// Apply one key event. Returns whether it was consumed.
    pub fn handle_key(&mut self, key: KeyEvent) -> bool {
        if key.kind != KeyEventKind::Press {
            return false;
        }
        match key.code {
            KeyCode::Enter => {
                self.submit();
                true
            }
            // Ctrl-C: clear a non-empty buffer, else quit (shell-like).
            KeyCode::Char('c') if key.modifiers.contains(KeyModifiers::CONTROL) => {
                if self.input.is_empty() {
                    self.should_exit = true;
                } else {
                    self.clear_input();
                }
                true
            }
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
            KeyCode::Esc => {
                if self.input.is_empty() {
                    self.should_exit = true;
                } else {
                    self.clear_input();
                }
                true
            }
            _ => false,
        }
    }

    // ─── line editing & command history ──────────────────────────────────

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

    /// Clear the input line and reset the cursor.
    fn clear_input(&mut self) {
        self.input.clear();
        self.cursor = 0;
    }

    /// Replace the input line, moving the cursor to its end.
    fn set_input(&mut self, line: String) {
        self.cursor = line.chars().count();
        self.input = line;
    }

    /// Submit the current line: record it in history, then dispatch it.
    fn submit(&mut self) {
        let line = std::mem::take(&mut self.input);
        self.cursor = 0;
        self.history_pos = None;
        self.draft.clear();
        let trimmed = line.trim();
        // Skip blanks and consecutive duplicates, matching a shell's history.
        if !trimmed.is_empty() && self.history.last().map(String::as_str) != Some(trimmed) {
            self.history.push(trimmed.to_string());
        }
        self.dispatch(&line);
    }

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

    // ─── read-only accessors (used by the widgets in Task 6) ──────────────

    pub fn input(&self) -> &str {
        &self.input
    }

    /// Char index of the edit cursor within the input line. Used to place the
    /// terminal cursor; a host embedding `CommandLine` uses this to position it.
    pub fn cursor(&self) -> usize {
        self.cursor
    }

    pub fn status(&self) -> &str {
        &self.status
    }

    /// Which listing the Program panel shows: the loaded program, or the log.
    pub fn active_listing(&self) -> (&'static str, &CodeView<String>) {
        if self.has_program {
            ("Program", &self.program)
        } else {
            ("Log", &self.log)
        }
    }

    /// The tableau rendering for the State panel.
    pub fn state_text(&self) -> String {
        match &self.machine {
            Some(m) => m.state_string(),
            None => "(no device — type `device N` or :load <file>)".to_string(),
        }
    }

    /// The measurement record as flat bits, events separated by spaces.
    pub fn measurement_bits(&self) -> String {
        match &self.machine {
            Some(m) => {
                let rec = m.measurement_record();
                if rec.is_empty() {
                    "(none)".to_string()
                } else {
                    format_record(&rec)
                }
            }
            None => "(none)".to_string(),
        }
    }

    /// A contextual footer hint.
    pub fn hint(&self) -> &'static str {
        if self.has_program && self.paused {
            "Enter=step  :c=continue  :reset  :q=quit"
        } else if self.machine.is_some() {
            "type a gate, or :load <file>   :q=quit"
        } else {
            ":load <file>  or  device N  to begin   :q=quit"
        }
    }

    /// Convenience full-screen composer for the standalone `ppvm` TUI. A host
    /// app (e.g. stellarscope) ignores this and lays out the individual
    /// `…View` widgets itself.
    pub fn render(&self, frame: &mut Frame) {
        let root = Layout::vertical([
            Constraint::Min(6),    // Program | State
            Constraint::Length(3), // measurement record
            Constraint::Length(2), // command line
        ])
        .split(frame.area());

        let top = Layout::horizontal([Constraint::Percentage(50), Constraint::Percentage(50)])
            .split(root[0]);

        frame.render_widget(ProgramView(self), top[0]);
        frame.render_widget(StateView(self), top[1]);
        frame.render_widget(RecordView(self), root[1]);
        frame.render_widget(CommandLine(self), root[2]);

        // Place the terminal cursor in the input line, just after the prompt,
        // clamped to the command area so a long line can't run off-panel.
        let cmd = root[2];
        let col = cmd.x + CommandLine::PROMPT.len() as u16 + self.cursor as u16;
        let x = col.min(cmd.x + cmd.width.saturating_sub(1));
        frame.set_cursor_position((x, cmd.y));
    }
}

/// Render a measurement record as flat bits: `Zero`→`0`, `One`→`1`, `Lost`→`2`
/// (the outcome's own enum value), events joined by spaces.
fn format_record(record: &[MeasurementResult]) -> String {
    record
        .iter()
        .map(|event| {
            event
                .iter()
                .map(|o| char::from(b'0' + *o as u8))
                .collect::<String>()
        })
        .collect::<Vec<_>>()
        .join(" ")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(code: KeyCode) -> KeyEvent {
        KeyEvent::new(code, KeyModifiers::NONE)
    }

    #[test]
    fn device_then_x_then_measure_records_one() {
        let mut app = AppState::new();
        app.dispatch("device 1");
        app.dispatch("x 0");
        app.dispatch("measure 0");
        assert_eq!(app.measurement_bits(), "1");
        assert!(app.status().contains("=> 1"), "status: {}", app.status());
    }

    #[test]
    fn fresh_measure_is_zero() {
        let mut app = AppState::new();
        app.dispatch("device 1");
        app.dispatch("measure 0");
        assert_eq!(app.measurement_bits(), "0");
    }

    #[test]
    fn gate_without_device_is_a_nonfatal_error() {
        let mut app = AppState::new();
        app.dispatch("x 0");
        assert!(app.status().contains("no device"));
        // Still usable afterwards.
        app.dispatch("device 1");
        app.dispatch("measure 0");
        assert_eq!(app.measurement_bits(), "0");
    }

    #[test]
    fn cnot_respects_control_target_order() {
        let mut app = AppState::new();
        app.dispatch("device 2");
        app.dispatch("x 0");
        app.dispatch("cnot 0 1");
        app.dispatch("measure 0");
        app.dispatch("measure 1");
        // Two separate measurement events, so two space-separated bits.
        assert_eq!(app.measurement_bits(), "1 1");
    }

    #[test]
    fn out_of_range_qubit_errors_not_panics() {
        let mut app = AppState::new();
        app.dispatch("device 1");
        app.dispatch("x 3");
        assert!(
            app.status().contains("out of range"),
            "status: {}",
            app.status()
        );
    }

    #[test]
    fn enter_key_dispatches_the_buffered_line() {
        let mut app = AppState::new();
        for c in "device 1".chars() {
            app.handle_key(key(KeyCode::Char(c)));
        }
        assert_eq!(app.input(), "device 1");
        app.handle_key(key(KeyCode::Enter));
        assert!(app.status().contains("1-qubit device"));
        assert!(app.input().is_empty(), "buffer should clear on submit");
    }

    #[test]
    fn backspace_edits_the_buffer() {
        let mut app = AppState::new();
        app.handle_key(key(KeyCode::Char('h')));
        app.handle_key(key(KeyCode::Char('i')));
        app.handle_key(key(KeyCode::Backspace));
        assert_eq!(app.input(), "h");
    }

    #[test]
    fn ctrl_c_on_empty_buffer_exits() {
        let mut app = AppState::new();
        app.handle_key(KeyEvent::new(KeyCode::Char('c'), KeyModifiers::CONTROL));
        assert!(app.should_exit);
    }

    #[test]
    fn ctrl_c_on_nonempty_buffer_clears_it() {
        let mut app = AppState::new();
        app.handle_key(key(KeyCode::Char('x')));
        app.handle_key(KeyEvent::new(KeyCode::Char('c'), KeyModifiers::CONTROL));
        assert!(!app.should_exit);
        assert!(app.input().is_empty());
    }

    #[test]
    fn esc_on_empty_buffer_exits() {
        let mut app = AppState::new();
        app.handle_key(key(KeyCode::Esc));
        assert!(app.should_exit);
    }

    #[test]
    fn esc_on_nonempty_buffer_clears_it() {
        let mut app = AppState::new();
        app.handle_key(key(KeyCode::Char('x')));
        app.handle_key(key(KeyCode::Esc));
        assert!(!app.should_exit);
        assert!(app.input().is_empty());
    }

    #[test]
    fn quit_command_sets_should_exit() {
        let mut app = AppState::new();
        app.dispatch(":q");
        assert!(app.should_exit);
    }

    /// A 1-qubit program with a breakpoint before measuring q0 (|0> -> 0).
    const BP_PROGRAM: &str = "device circuit.n_qubits 1;\n\
                              fn @main() { breakpoint\n const.u64 0\n circuit.measure\n ret }\n";

    #[test]
    fn load_source_starts_paused_with_a_listing() {
        let mut app = AppState::new();
        app.load_source(BP_PROGRAM).unwrap();
        assert!(app.has_program());
        assert!(app.paused());
        let (title, view) = app.active_listing();
        assert_eq!(title, "Program");
        assert!(!view.lines().is_empty());
        assert_eq!(view.cursor(), Some(0), "cursor starts at pc 0");
    }

    #[test]
    fn continue_pauses_at_breakpoint_then_finishes() {
        let mut app = AppState::new();
        app.load_source(BP_PROGRAM).unwrap();
        app.dispatch(":c");
        assert!(
            app.status().contains("breakpoint"),
            "status: {}",
            app.status()
        );
        app.dispatch(":c");
        assert!(
            app.status().contains("finished"),
            "status: {}",
            app.status()
        );
        // |0> measured is 0.
        assert_eq!(app.measurement_bits(), "0");
    }

    #[test]
    fn empty_line_steps_and_advances_cursor() {
        let mut app = AppState::new();
        app.load_source(BP_PROGRAM).unwrap();
        let start = app.active_listing().1.cursor();
        app.dispatch(""); // empty line == step
        let after = app.active_listing().1.cursor();
        assert_ne!(start, after, "stepping should move the cursor");
    }

    #[test]
    fn inject_gate_at_breakpoint_then_resume() {
        // At the breakpoint, inject X on q0; resuming, the program measures |1>.
        let mut app = AppState::new();
        app.load_source(BP_PROGRAM).unwrap();
        app.dispatch(":c"); // run to the breakpoint
        assert!(app.status().contains("breakpoint"));
        app.dispatch("x 0"); // inject while paused
        app.dispatch(":c"); // resume; program measures q0
        assert!(
            app.status().contains("finished"),
            "status: {}",
            app.status()
        );
        assert_eq!(
            app.measurement_bits(),
            "1",
            "injected X should flip the result"
        );
    }

    #[test]
    fn reset_returns_a_program_to_the_start() {
        let mut app = AppState::new();
        app.load_source(BP_PROGRAM).unwrap();
        app.dispatch(":c");
        app.dispatch(":c"); // finished
        app.dispatch(":reset");
        assert!(app.paused());
        assert_eq!(app.active_listing().1.cursor(), Some(0));
        assert_eq!(app.measurement_bits(), "(none)");
    }

    // ─── line editing & history ──────────────────────────────────────────

    /// Type each char, then press Enter (submitting to history + dispatch).
    fn type_line(app: &mut AppState, line: &str) {
        for c in line.chars() {
            app.handle_key(key(KeyCode::Char(c)));
        }
        app.handle_key(key(KeyCode::Enter));
    }

    #[test]
    fn typing_inserts_and_advances_the_cursor() {
        let mut app = AppState::new();
        for c in "hz".chars() {
            app.handle_key(key(KeyCode::Char(c)));
        }
        assert_eq!(app.input(), "hz");
        assert_eq!(app.cursor(), 2);
    }

    #[test]
    fn left_right_home_end_move_the_cursor() {
        let mut app = AppState::new();
        for c in "abc".chars() {
            app.handle_key(key(KeyCode::Char(c)));
        }
        assert_eq!(app.cursor(), 3);
        app.handle_key(key(KeyCode::Left));
        assert_eq!(app.cursor(), 2);
        app.handle_key(key(KeyCode::Home));
        assert_eq!(app.cursor(), 0);
        app.handle_key(key(KeyCode::Left)); // saturates at 0
        assert_eq!(app.cursor(), 0);
        app.handle_key(key(KeyCode::End));
        assert_eq!(app.cursor(), 3);
        app.handle_key(key(KeyCode::Right)); // saturates at len
        assert_eq!(app.cursor(), 3);
    }

    #[test]
    fn insert_and_delete_at_the_cursor() {
        let mut app = AppState::new();
        for c in "ac".chars() {
            app.handle_key(key(KeyCode::Char(c)));
        }
        app.handle_key(key(KeyCode::Left)); // cursor between 'a' and 'c'
        app.handle_key(key(KeyCode::Char('b')));
        assert_eq!(app.input(), "abc");
        assert_eq!(app.cursor(), 2);
        app.handle_key(key(KeyCode::Backspace)); // removes 'b' before cursor
        assert_eq!(app.input(), "ac");
        assert_eq!(app.cursor(), 1);
        app.handle_key(key(KeyCode::Delete)); // removes 'c' at cursor
        assert_eq!(app.input(), "a");
        assert_eq!(app.cursor(), 1);
    }

    #[test]
    fn up_arrow_recalls_previous_commands() {
        let mut app = AppState::new();
        type_line(&mut app, "device 1");
        type_line(&mut app, "x 0");
        app.handle_key(key(KeyCode::Up)); // newest first
        assert_eq!(app.input(), "x 0");
        assert_eq!(app.cursor(), 3, "cursor lands at end of the recalled line");
        app.handle_key(key(KeyCode::Up));
        assert_eq!(app.input(), "device 1");
        app.handle_key(key(KeyCode::Up)); // already oldest, stays put
        assert_eq!(app.input(), "device 1");
    }

    #[test]
    fn down_arrow_returns_toward_the_live_line() {
        let mut app = AppState::new();
        type_line(&mut app, "device 1");
        type_line(&mut app, "x 0");
        for c in "meas".chars() {
            app.handle_key(key(KeyCode::Char(c)));
        }
        app.handle_key(key(KeyCode::Up)); // stash "meas", show "x 0"
        app.handle_key(key(KeyCode::Up)); // "device 1"
        app.handle_key(key(KeyCode::Down)); // "x 0"
        assert_eq!(app.input(), "x 0");
        app.handle_key(key(KeyCode::Down)); // past newest -> restore draft
        assert_eq!(app.input(), "meas");
        assert_eq!(app.cursor(), 4);
    }

    #[test]
    fn history_skips_consecutive_duplicates() {
        let mut app = AppState::new();
        type_line(&mut app, "device 1");
        type_line(&mut app, "device 1"); // same command twice -> one entry
        app.handle_key(key(KeyCode::Up));
        assert_eq!(app.input(), "device 1");
        app.handle_key(key(KeyCode::Up)); // only one entry -> stays
        assert_eq!(app.input(), "device 1");
    }

    #[test]
    fn render_places_the_terminal_cursor_after_the_prompt() {
        use crate::widgets::CommandLine;
        use ratatui::Terminal;
        use ratatui::backend::TestBackend;

        let mut app = AppState::new();
        for c in "ab".chars() {
            app.handle_key(key(KeyCode::Char(c)));
        }
        let mut terminal = Terminal::new(TestBackend::new(80, 24)).unwrap();
        terminal.draw(|f| app.render(f)).unwrap();

        // Command area starts at x=0, so the cursor sits at prompt width + 2.
        let pos = terminal.get_cursor_position().unwrap();
        assert_eq!(pos.x, (CommandLine::PROMPT.len() + 2) as u16);
    }
}
