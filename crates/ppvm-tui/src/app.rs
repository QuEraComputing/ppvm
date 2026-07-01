// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `AppState` — the terminal-agnostic state of the ppvm TUI. Owns an optional
//! [`PPVM`] plus the command buffer, status line, program listing, and REPL
//! log. `dispatch` runs one command string; `handle_key` edits the buffer and
//! submits on Enter. Nothing here touches a terminal or runs a loop.

use crossterm::event::{KeyCode, KeyEvent, KeyEventKind, KeyModifiers};
use eyre::{Result, eyre};
use ppvm_vihaco::CircuitInstruction;
use ppvm_vihaco::composite::PPVM;
use ppvm_vihaco::measurements::MeasurementResult;

use crate::codeview::CodeView;
use crate::command::{Command, parse_command};

/// Terminal-agnostic state for the ppvm TUI.
pub struct AppState {
    /// The live machine. `None` until `device N` or a program is loaded.
    machine: Option<PPVM>,
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
            program: CodeView::new(),
            log: CodeView::new(),
            n_qubits: 0,
            has_program: false,
            paused: false,
            finished: false,
            input: String::new(),
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
            // Step/Continue/Reset/Load are implemented in Task 5.
            Command::Step | Command::Continue | Command::Reset | Command::Load(_) => {
                self.set_status("not a REPL command (load a program first)");
                Ok(())
            }
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
                let line = std::mem::take(&mut self.input);
                self.dispatch(&line);
                true
            }
            // Ctrl-C: clear a non-empty buffer, else quit (shell-like).
            KeyCode::Char('c') if key.modifiers.contains(KeyModifiers::CONTROL) => {
                if self.input.is_empty() {
                    self.should_exit = true;
                } else {
                    self.input.clear();
                }
                true
            }
            KeyCode::Char(c) => {
                self.input.push(c);
                true
            }
            KeyCode::Backspace => {
                self.input.pop();
                true
            }
            KeyCode::Esc => {
                if self.input.is_empty() {
                    self.should_exit = true;
                } else {
                    self.input.clear();
                }
                true
            }
            _ => false,
        }
    }

    // ─── read-only accessors (used by the widgets in Task 6) ──────────────

    pub fn input(&self) -> &str {
        &self.input
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
    fn quit_command_sets_should_exit() {
        let mut app = AppState::new();
        app.dispatch(":q");
        assert!(app.should_exit);
    }
}
