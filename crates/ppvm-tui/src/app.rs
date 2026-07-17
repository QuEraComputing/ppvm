// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `AppState` — the terminal-agnostic state of the ppvm TUI. Owns an optional
//! [`PPVM`] plus a [`LineEditor`], status line, program listing, and REPL log.
//! `dispatch` runs one command string; `handle_key` routes app-level keys
//! (submit, quit) and delegates editing to the [`LineEditor`]. Nothing here
//! touches a terminal or runs a loop.

use eyre::{Result, eyre};
use ppvm_vihaco::component::{ComplexityMetric, ComplexityMetricKind};
use ppvm_vihaco::composite::{PPVM, StepOutcome};
use ppvm_vihaco::measurements::MeasurementResult;
use ppvm_vihaco::{CircuitInstruction, PPVMModule, compile_program, load_module_file};
use ratatui::Frame;
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyEventKind, KeyModifiers};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::widgets::Clear;

use crate::codeview::CodeView;
use crate::command::{Command, parse_command};
use crate::editor::LineEditor;
use crate::widgets::{
    CommandLine, ComplexityTreeView, HelpOverlay, ProgramView, RecordView, StackView, StateView,
};

const STATE_DETAIL_QUBIT_LIMIT: usize = 16;
const STATE_DETAIL_COMPLEXITY_LIMIT: usize = 128;
const TREE_LAYER_STRIDE: usize = 6;
const DEFAULT_TREE_HEIGHT: u16 = 8;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ComplexityLayer {
    kind: ComplexityMetricKind,
    count: usize,
}

impl From<ComplexityMetric> for ComplexityLayer {
    fn from(metric: ComplexityMetric) -> Self {
        Self {
            kind: metric.kind,
            count: metric.count,
        }
    }
}

#[derive(Debug, Default)]
struct ComplexityHistory {
    layers: Vec<ComplexityLayer>,
}

impl ComplexityHistory {
    fn clear(&mut self) {
        self.layers.clear();
    }

    fn push_if_changed(&mut self, metric: ComplexityMetric) {
        let layer = ComplexityLayer::from(metric);
        if self.layers.last() != Some(&layer) {
            self.layers.push(layer);
        }
    }

    fn counts(&self) -> Vec<usize> {
        self.layers.iter().map(|layer| layer.count).collect()
    }

    fn visible_layer_count(&self, width: usize) -> usize {
        (width.div_ceil(TREE_LAYER_STRIDE))
            .max(1)
            .min(self.layers.len())
    }

    fn render(&self, width: u16, height: u16) -> String {
        let Some(current) = self.layers.last() else {
            return "(no device)".to_string();
        };

        let width = usize::from(width).max(1);
        let graph_height = usize::from(height.saturating_sub(1)).max(1);
        let visible_layers = self.visible_layer_count(width);
        let start = self.layers.len() - visible_layers;
        let end = self.layers.len() - 1;
        let layers = &self.layers[start..];

        let mut canvas = vec![vec![' '; width]; graph_height];
        let layer_rows: Vec<Vec<usize>> = layers
            .iter()
            .map(|layer| spread_positions(layer.count, graph_height))
            .collect();
        let layer_x: Vec<usize> = (0..layers.len())
            .map(|idx| idx * TREE_LAYER_STRIDE)
            .filter(|&x| x < width)
            .collect();

        for idx in 0..layer_x.len().saturating_sub(1) {
            draw_transposed_connector(
                &mut canvas,
                layer_x[idx],
                layer_x[idx + 1],
                &layer_rows[idx],
                &layer_rows[idx + 1],
            );
        }

        for (&x, rows) in layer_x.iter().zip(layer_rows.iter()) {
            for &row in rows {
                canvas[row][x] = '●';
            }
        }

        let mut lines = vec![format!(
            "layers {start:03}..{end:03} | current: {} {}",
            current.count,
            current.kind.noun(current.count)
        )];
        lines.extend(
            canvas
                .into_iter()
                .map(|row| row.into_iter().collect::<String>()),
        );
        lines.join("\n")
    }
}

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
    /// True while the debugger is stopped and waiting for input (at start, after
    /// a step, or at a breakpoint).
    paused: bool,
    /// True once the loaded program has run to Return/Halt.
    finished: bool,
    /// The command line: buffer, edit cursor, and history.
    editor: LineEditor,
    /// Whether the State panel should render full state details when safe.
    show_state_details: bool,
    /// Branching-complexity history shown in the tree panel.
    complexity: ComplexityHistory,
    /// The status/error line.
    status: String,
    /// Whether the help overlay is currently shown.
    show_help: bool,
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
            editor: LineEditor::new(),
            show_state_details: true,
            complexity: ComplexityHistory::default(),
            status: String::new(),
            show_help: false,
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
            Command::Help => {
                self.show_help = !self.show_help;
                Ok(())
            }
            Command::ToggleState => {
                self.show_state_details = !self.show_state_details;
                let mode = if self.show_state_details {
                    "enabled"
                } else {
                    "hidden"
                };
                self.set_status(format!("state details {mode}"));
                Ok(())
            }
        }
    }

    fn new_device(&mut self, n: usize) -> Result<()> {
        self.machine = Some(PPVM::with_qubits(n)?);
        // Forget any previously loaded program so `:reset` resets this device
        // rather than resurrecting the old program.
        self.module = None;
        self.n_qubits = n;
        self.has_program = false;
        self.paused = false;
        self.finished = false;
        self.program.clear();
        self.reset_complexity_history();
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
        if self.load_module(module) {
            self.set_status("loaded program");
        }
        Ok(())
    }

    fn load_file(&mut self, path: &str) -> Result<()> {
        let module = load_module_file(path).map_err(|e| eyre!("failed to load {path}: {e}"))?;
        if self.load_module(module) {
            self.set_status(format!("loaded {path}"));
        }
        Ok(())
    }

    /// Core loader: rebuild the machine from `module` and pause at pc 0.
    /// Returns `false` if load/init fails (error is written to the status line).
    fn load_module(&mut self, module: PPVMModule) -> bool {
        let mut m = PPVM::default();
        // A fresh machine + load + init gives clean tableau/record state; these
        // only fail on malformed modules, which `compile_program` already
        // rejects, so surface as a status rather than unwinding the UI.
        if let Err(e) = m.load(&module).and_then(|()| m.init()) {
            self.set_status(format!("error: {e}"));
            return false;
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
        self.reset_complexity_history();
        self.refresh_cursor();
        true
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
        self.record_complexity();
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
            self.record_complexity();
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
            StepOutcome::Continue => {
                self.set_status("");
            }
            StepOutcome::Breakpoint => {
                self.paused = true;
                self.set_status("-- breakpoint hit --");
            }
            StepOutcome::Return | StepOutcome::Halt => {
                self.paused = false;
                self.finished = true;
                self.set_status("program finished");
            }
        }
    }

    fn reset(&mut self) -> Result<()> {
        if let Some(module) = self.module.clone() {
            if self.load_module(module) {
                self.set_status("reset");
            }
        } else if self.n_qubits > 0 {
            self.machine = Some(PPVM::with_qubits(self.n_qubits)?);
            self.reset_complexity_history();
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
        self.record_complexity();
        Ok(())
    }

    fn set_status(&mut self, s: impl Into<String>) {
        self.status = s.into();
    }

    fn reset_complexity_history(&mut self) {
        self.complexity.clear();
        self.record_complexity();
    }

    fn record_complexity(&mut self) {
        if let Some(machine) = &self.machine {
            self.complexity.push_if_changed(machine.complexity_metric());
        }
    }

    // ─── key handling ────────────────────────────────────────────────────

    /// Apply one key event. Returns whether it was consumed. App-level keys
    /// (submit, quit) are handled here; the rest go to the line editor.
    pub fn handle_key(&mut self, key: KeyEvent) -> bool {
        if key.kind != KeyEventKind::Press {
            return false;
        }
        match key.code {
            KeyCode::Enter => {
                self.submit();
                true
            }
            // Ctrl-C / Esc: clear a non-empty buffer, else quit (shell-like).
            KeyCode::Char('c') if key.modifiers.contains(KeyModifiers::CONTROL) => {
                self.clear_or_quit();
                true
            }
            KeyCode::Esc => {
                self.clear_or_quit();
                true
            }
            _ => self.editor.handle_key(key),
        }
    }

    /// Clear a non-empty command line, or quit when it is already empty.
    fn clear_or_quit(&mut self) {
        if self.editor.is_empty() {
            self.should_exit = true;
        } else {
            self.editor.clear();
        }
    }

    /// Submit the current line: record it in history (via the editor), then
    /// dispatch it.
    fn submit(&mut self) {
        let line = self.editor.submit();
        self.dispatch(&line);
    }

    // ─── read-only accessors (used by the widgets in Task 6) ──────────────

    pub fn input(&self) -> &str {
        self.editor.input()
    }

    /// Char index of the edit cursor within the input line. Used to place the
    /// terminal cursor; a host embedding `CommandLine` uses this to position it.
    pub fn cursor(&self) -> usize {
        self.editor.cursor()
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
            Some(m) => {
                let summary = state_summary(m);
                let metric = m.complexity_metric();
                if !self.show_state_details {
                    return format!("state details hidden (:state to show)\n{summary}");
                }
                if m.n_qubits() > STATE_DETAIL_QUBIT_LIMIT
                    || metric.count > STATE_DETAIL_COMPLEXITY_LIMIT
                {
                    return format!(
                        "state details suppressed (limit: {STATE_DETAIL_QUBIT_LIMIT} qubits, {STATE_DETAIL_COMPLEXITY_LIMIT} {})\n{summary}",
                        metric.kind.noun(STATE_DETAIL_COMPLEXITY_LIMIT),
                    );
                }
                format!("{summary}\n{}", m.compact_state_string())
            }
            None => "(no device — type `device N` or :load <file>)".to_string(),
        }
    }

    /// Current CPU operand stack, top entry first.
    pub fn stack_text(&self) -> String {
        let Some(machine) = &self.machine else {
            return "(no device)".to_string();
        };
        let stack = machine.stack_snapshot();
        if stack.is_empty() {
            return "(empty)".to_string();
        }
        stack
            .iter()
            .enumerate()
            .rev()
            .map(|(idx, value)| {
                let marker = if idx + 1 == stack.len() { "top" } else { "   " };
                format!("{idx:04} {marker} {value}")
            })
            .collect::<Vec<_>>()
            .join("\n")
    }

    /// Exact complexity counts recorded in the tree, exposed for tests and
    /// embedders that want to render the history themselves.
    pub fn complexity_counts(&self) -> Vec<usize> {
        self.complexity.counts()
    }

    pub fn complexity_graph_text(&self, width: u16) -> String {
        self.complexity.render(width, DEFAULT_TREE_HEIGHT)
    }

    pub fn complexity_graph_text_for_area(&self, width: u16, height: u16) -> String {
        self.complexity.render(width, height)
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
        if self.has_program && self.finished {
            ":reset=run again  ·  :load <file>  ·  :help  ·  :q=quit"
        } else if self.has_program && self.paused {
            "Enter=step  :c=continue  :reset  :help  :q=quit"
        } else if self.machine.is_some() {
            "type a gate, or :load <file>   :help  :q=quit"
        } else {
            ":load <file>  ·  device N  ·  :help  ·  :q=quit"
        }
    }

    /// Whether the help overlay should be drawn.
    pub fn show_help(&self) -> bool {
        self.show_help
    }

    /// The command reference shown in the help overlay.
    pub fn help_text(&self) -> &'static str {
        HELP_TEXT
    }

    /// Convenience full-screen composer for the standalone `ppvm` TUI. A host
    /// app (e.g. stellarscope) ignores this and lays out the individual
    /// `…View` widgets itself.
    pub fn render(&self, frame: &mut Frame) {
        let root = Layout::vertical([
            Constraint::Ratio(3, 5), // Program | complexity tree
            Constraint::Ratio(2, 5), // Stack | State
            Constraint::Length(3),   // measurement record
            Constraint::Length(2),   // command line
        ])
        .split(frame.area());

        let top = Layout::horizontal([Constraint::Percentage(45), Constraint::Percentage(55)])
            .split(root[0]);
        let middle = Layout::horizontal([Constraint::Percentage(45), Constraint::Percentage(55)])
            .split(root[1]);

        frame.render_widget(ProgramView(self), top[0]);
        frame.render_widget(ComplexityTreeView(self), top[1]);
        frame.render_widget(StackView(self), middle[0]);
        frame.render_widget(StateView(self), middle[1]);
        frame.render_widget(RecordView(self), root[2]);
        frame.render_widget(CommandLine(self), root[3]);

        // Place the terminal cursor in the input line, just after the prompt,
        // clamped to the command area so a long line can't run off-panel.
        let cmd = root[3];
        let col = cmd.x + CommandLine::PROMPT.len() as u16 + self.editor.cursor() as u16;
        let x = col.min(cmd.x + cmd.width.saturating_sub(1));
        frame.set_cursor_position((x, cmd.y));

        // The help overlay floats above everything else when toggled on.
        if self.show_help {
            let full = frame.area();
            let w = 76.min(full.width.saturating_sub(2));
            let h = 20.min(full.height.saturating_sub(2));
            let popup = Rect {
                x: full.x + full.width.saturating_sub(w) / 2,
                y: full.y + full.height.saturating_sub(h) / 2,
                width: w,
                height: h,
            };
            frame.render_widget(Clear, popup);
            frame.render_widget(HelpOverlay(self), popup);
        }
    }
}

/// The command reference shown by `:help`. Mirrors the command grammar in
/// [`crate::command`]: bare tokens are gate ops, `:`-prefixed are meta/debug.
const HELP_TEXT: &str = "\
Meta / debug
  device N              create a fresh N-qubit tableau device
  :load <file>          load a .sst / .ssb program (paused at start)
  Enter (empty)  :s     step one instruction
  :continue  :c         run to the next breakpoint or the end
  :reset                restart the loaded program / device
  :state                toggle detailed state rendering
  :help  :h             toggle this help
  :quit  :q  (Ctrl-C)   leave

Gates  (q = qubit index; angles / probabilities are floats)
  x y z h s sadj sqrtx sqrty sqrtxadj sqrtyadj t tadj reset measure <q>
  cnot <c> <t>     cz <a> <b>
  rx ry rz <q> <angle>     r <q> <axis> <angle>
  u3 <q> <theta> <phi> <lam>
  rxx ryy rzz <a> <b> <angle>
  depolarize loss <q> <p>     depolarize2 <a> <b> <p>
  paulierror <q> <px> <py> <pz>     correlatedloss <a> <b> <p0> <p1> <p2>

Line editing: ←/→ move · Home/End · Backspace/Del · ↑/↓ history";

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

fn state_summary(machine: &PPVM) -> String {
    let metric = machine.complexity_metric();
    format!(
        "{} backend | {} {} | {} {}",
        machine.backend_name(),
        machine.n_qubits(),
        plural(machine.n_qubits(), "qubit", "qubits"),
        metric.count,
        metric.kind.noun(metric.count)
    )
}

fn plural(count: usize, one: &'static str, many: &'static str) -> &'static str {
    if count == 1 { one } else { many }
}

fn spread_positions(count: usize, height: usize) -> Vec<usize> {
    if count == 0 {
        return Vec::new();
    }
    let visible = count.min(height.max(1));
    if visible == 1 {
        return vec![height.saturating_sub(1) / 2];
    }

    let full_span = height.saturating_sub(1);
    let inner_span = height.saturating_sub(3);
    let inner_capacity = height.saturating_sub(2);
    let desired_span = if visible <= inner_capacity && inner_span > 0 {
        if visible == 2 {
            2.min(inner_span)
        } else {
            (visible - 1).min(inner_span)
        }
    } else {
        full_span
    };
    let top = full_span.saturating_sub(desired_span) / 2;
    let mut out = Vec::with_capacity(visible);
    for idx in 0..visible {
        let pos = top + idx * desired_span / (visible - 1);
        if out.last() != Some(&pos) {
            out.push(pos);
        }
    }
    out
}

fn draw_transposed_connector(
    canvas: &mut [Vec<char>],
    from_x: usize,
    to_x: usize,
    from_rows: &[usize],
    to_rows: &[usize],
) {
    if from_rows.is_empty() || to_rows.is_empty() || from_x >= to_x {
        return;
    }
    let bus_x = (from_x + to_x) / 2;
    let min_row = from_rows
        .iter()
        .chain(to_rows.iter())
        .min()
        .copied()
        .unwrap_or(0);
    let max_row = from_rows
        .iter()
        .chain(to_rows.iter())
        .max()
        .copied()
        .unwrap_or(min_row);

    for &row in from_rows {
        draw_canvas_horizontal(canvas, row, from_x + 1, bus_x);
    }
    for &row in to_rows {
        draw_canvas_horizontal(canvas, row, bus_x, to_x.saturating_sub(1));
    }
    for row in min_row..=max_row {
        let glyph = if row == min_row {
            '╭'
        } else if row == max_row {
            '╰'
        } else {
            '│'
        };
        put_canvas_connector(canvas, row, bus_x, glyph);
    }
    for &row in from_rows {
        put_canvas_connector(canvas, row, bus_x, '┤');
    }
    for &row in to_rows {
        put_canvas_connector(canvas, row, bus_x, '├');
    }
}

fn draw_canvas_horizontal(canvas: &mut [Vec<char>], row: usize, start: usize, end: usize) {
    if start > end {
        return;
    }
    for idx in start..=end {
        put_canvas_connector(canvas, row, idx, '─');
    }
}

fn put_canvas_connector(canvas: &mut [Vec<char>], row: usize, col: usize, glyph: char) {
    let Some(line) = canvas.get_mut(row) else {
        return;
    };
    let Some(slot) = line.get_mut(col) else {
        return;
    };
    *slot = merge_connector(*slot, glyph);
}

fn merge_connector(existing: char, incoming: char) -> char {
    match (existing, incoming) {
        (' ', glyph) => glyph,
        ('─', glyph) | (glyph, '─') => glyph,
        (same, glyph) if same == glyph => same,
        ('╭', '├') | ('├', '╭') => '┌',
        ('╰', '├') | ('├', '╰') => '└',
        ('╭', '┤') | ('┤', '╭') => '┐',
        ('╰', '┤') | ('┤', '╰') => '┘',
        ('│', '├') | ('├', '│') => '├',
        ('│', '┤') | ('┤', '│') => '┤',
        ('├', '┤') | ('┤', '├') => '┼',
        ('│', '┴') | ('┴', '│') => '┴',
        ('│', '┬') | ('┬', '│') => '┬',
        ('│', _) | (_, '│') => '┼',
        ('╭' | '╮' | '╰' | '╯', '┴' | '┬' | '┼')
        | ('┴' | '┬' | '┼', '╭' | '╮' | '╰' | '╯')
        | ('┴', '┬')
        | ('┬', '┴') => '┼',
        (_, glyph) => glyph,
    }
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
    fn state_command_toggles_state_details() {
        let mut app = AppState::new();
        app.dispatch("device 1");
        assert!(app.state_text().contains("Tableau"));

        app.dispatch(":state");
        assert!(
            app.state_text().contains("state details hidden"),
            "state text: {}",
            app.state_text()
        );

        app.dispatch(":state");
        assert!(
            !app.state_text().contains("state details hidden"),
            "state text: {}",
            app.state_text()
        );
    }

    #[test]
    fn large_device_state_is_summarized_without_detail_formatting() {
        let mut app = AppState::new();
        app.dispatch("device 17");
        let state = app.state_text();
        assert!(
            state.contains("state details suppressed"),
            "state text: {state}"
        );
        assert!(state.contains("17 qubits"), "state text: {state}");
        assert!(state.contains("coefficient"), "state text: {state}");
    }

    #[test]
    fn ten_qubit_state_uses_compact_side_by_side_tableau_format() {
        let mut app = AppState::new();
        app.dispatch("device 10");

        let state = app.state_text();
        assert!(
            state.lines().count() <= 12,
            "state text should be compact enough for a 10-qubit tableau:\n{state}"
        );
        assert!(
            state.lines().any(|line| {
                line.contains("q00")
                    && line.contains(" D ")
                    && line.contains(" S ")
                    && line.contains("Index 0:")
            }),
            "coefficient column should share a row with paired tableau rows:\n{state}"
        );
        assert!(
            !state
                .lines()
                .any(|line| line.trim_start().starts_with("Index 0:")),
            "coefficients should not be stacked below the tableau:\n{state}"
        );
    }

    #[test]
    fn stack_text_tracks_current_cpu_stack() {
        const STACK_PROGRAM: &str = "device circuit.n_qubits 1;\n\
                                    fn @main() { const.u64 7\n breakpoint\n ret }\n";
        let mut app = AppState::new();
        app.load_source(STACK_PROGRAM).unwrap();

        assert_eq!(app.stack_text(), "(empty)");
        app.dispatch(""); // const.u64 7

        let stack = app.stack_text();
        assert!(stack.contains("top"), "stack text: {stack}");
        assert!(stack.contains("7"), "stack text: {stack}");
    }

    #[test]
    fn complexity_history_records_only_count_changes() {
        let mut app = AppState::new();
        app.dispatch("device 1");
        assert_eq!(app.complexity_counts(), vec![1]);

        app.dispatch("h 0");
        assert_eq!(
            app.complexity_counts(),
            vec![1],
            "Clifford H should not add a layer when coefficient count is unchanged"
        );

        app.dispatch("t 0");
        assert_eq!(app.complexity_counts(), vec![1, 2]);

        app.dispatch("measure 0");
        assert_eq!(
            app.complexity_counts(),
            vec![1, 2, 1],
            "measurement should record the shrink back to one coefficient"
        );
    }

    #[test]
    fn complexity_tree_uses_unicode_connectors() {
        let mut app = AppState::new();
        app.dispatch("device 2");
        app.dispatch("h 0");
        app.dispatch("h 1");
        app.dispatch("t 0");
        app.dispatch("t 1");
        assert_eq!(app.complexity_counts(), vec![1, 2, 4]);

        let graph = app.complexity_graph_text(72);
        assert!(graph.contains('●'), "graph text: {graph}");
        assert!(graph.contains('┌'), "graph text: {graph}");
        assert!(graph.contains('└'), "graph text: {graph}");
        assert!(graph.contains('├'), "graph text: {graph}");
        assert!(
            graph.contains('┤') || graph.contains('┼'),
            "graph text: {graph}"
        );
        assert!(
            !graph.lines().skip(1).any(|line| line.contains('o')),
            "graph text: {graph}"
        );
        assert!(!graph.contains('/'), "graph text: {graph}");
        assert!(!graph.contains('\\'), "graph text: {graph}");
    }

    #[test]
    fn complexity_tree_is_transposed_and_autoscrolls_to_latest_layers() {
        let mut app = AppState::new();
        app.complexity.clear();
        for count in [1, 2, 4, 2, 5, 3, 6, 1] {
            app.complexity.push_if_changed(ComplexityMetric {
                kind: ComplexityMetricKind::Coefficients,
                count,
            });
        }

        let graph = app.complexity_graph_text(24);
        assert!(graph.contains("layers 004..007"), "graph text: {graph}");
        assert!(!graph.contains("000"), "graph text: {graph}");
        assert!(!graph.contains("001"), "graph text: {graph}");
        assert!(
            !graph.lines().any(|line| line.starts_with("000      1")),
            "graph text: {graph}"
        );
        assert!(graph.contains("004"), "graph text: {graph}");
        assert!(graph.contains("007"), "graph text: {graph}");
    }

    #[test]
    fn complexity_tree_small_branching_stays_near_center() {
        let mut app = AppState::new();
        app.complexity.clear();
        for count in [1, 2] {
            app.complexity.push_if_changed(ComplexityMetric {
                kind: ComplexityMetricKind::Coefficients,
                count,
            });
        }

        let graph_height = 7;
        let graph = app.complexity_graph_text_for_area(24, graph_height + 2);
        let node_rows: Vec<usize> = graph
            .lines()
            .skip(1)
            .take(graph_height as usize)
            .enumerate()
            .filter_map(|(row, line)| line.contains('●').then_some(row))
            .collect();

        assert!(
            node_rows.contains(&(graph_height as usize / 2)),
            "initial node should stay centered:\n{graph}"
        );
        assert!(
            !node_rows.contains(&0) && !node_rows.contains(&(graph_height as usize - 1)),
            "two-way branching should not jump straight to the borders:\n{graph}"
        );
    }

    #[test]
    fn complexity_tree_centers_initial_node_in_short_even_viewport() {
        let mut app = AppState::new();
        app.dispatch("device 1");

        let graph_height = 5;
        let graph = app.complexity_graph_text_for_area(24, graph_height + 1);
        let node_row = graph
            .lines()
            .skip(1)
            .take(graph_height as usize)
            .position(|line| line.contains('●'))
            .expect("initial node should be visible");

        assert!(
            node_row <= (graph_height as usize - 1) / 2,
            "initial node should use the upper center row in a short even viewport:\n{graph}"
        );
    }

    #[test]
    fn complexity_tree_initial_view_does_not_start_with_bottom_label() {
        let mut app = AppState::new();
        app.dispatch("device 1");

        let graph = app.complexity_graph_text_for_area(24, 8);
        assert!(
            !graph.lines().last().unwrap_or("").contains("000:1"),
            "initial tree should not draw its layer label at bottom-left:\n{graph}"
        );
    }

    #[test]
    fn rendered_complexity_tree_initial_node_is_vertically_centered() {
        use ratatui::Terminal;
        use ratatui::backend::TestBackend;

        let mut app = AppState::new();
        app.dispatch("device 1");

        let mut terminal = Terminal::new(TestBackend::new(100, 18)).unwrap();
        terminal.draw(|f| app.render(f)).unwrap();
        let lines = buffer_lines(terminal.backend().buffer(), 100);

        let tree_top = row_containing(&lines, "Complexity tree").unwrap();
        let tree_bottom = lines
            .iter()
            .enumerate()
            .skip(tree_top + 1)
            .find_map(|(row, line)| line.contains('┘').then_some(row))
            .unwrap();
        let node_row = lines
            .iter()
            .enumerate()
            .skip(tree_top + 1)
            .take(tree_bottom - tree_top)
            .find_map(|(row, line)| line.contains('●').then_some(row))
            .unwrap();

        let panel_mid = tree_top + (tree_bottom - tree_top) / 2;
        let lower_slack = (tree_bottom - tree_top) / 4;
        assert!(
            node_row <= panel_mid + lower_slack,
            "initial node should not render near the bottom of the tree panel\n{}",
            lines.join("\n")
        );
    }

    #[test]
    fn standalone_layout_puts_tree_above_state_on_the_right() {
        use ratatui::Terminal;
        use ratatui::backend::TestBackend;

        let mut app = AppState::new();
        app.dispatch("device 2");
        app.dispatch("h 0");
        app.dispatch("t 0");

        let mut terminal = Terminal::new(TestBackend::new(120, 36)).unwrap();
        terminal.draw(|f| app.render(f)).unwrap();
        let lines = buffer_lines(terminal.backend().buffer(), 120);

        let tree_row = row_containing(&lines, "Complexity tree").unwrap();
        let state_row = row_containing(&lines, "State").unwrap();
        assert!(
            tree_row < state_row,
            "tree should be above state\n{}",
            lines.join("\n")
        );
        assert!(
            lines[tree_row].find("Complexity tree").unwrap() > 50,
            "tree should be in the right pane\n{}",
            lines.join("\n")
        );
        assert!(
            lines[state_row].find("State").unwrap() > 50,
            "state should remain in the right pane\n{}",
            lines.join("\n")
        );
    }

    #[test]
    fn complexity_history_uses_paulisum_term_counts() {
        const PAULISUM_PROGRAM: &str = "device circuit.n_qubits 1;\n\
                                       device circuit.backend paulisum;\n\
                                       device circuit.observable Z;\n\
                                       fn @main() {\n\
                                         const.u64 0\n\
                                         const.f64 0.7\n\
                                         circuit.ry\n\
                                         ret\n\
                                       }\n";
        let mut app = AppState::new();
        app.load_source(PAULISUM_PROGRAM).unwrap();
        assert_eq!(app.complexity_counts(), vec![1]);

        app.dispatch(""); // const.u64 0
        app.dispatch(""); // const.f64 0.7
        app.dispatch(""); // circuit.ry: Z -> cos(theta) Z + sin(theta) X

        assert_eq!(app.complexity_counts(), vec![1, 2]);
        assert!(
            app.state_text().contains("2 terms"),
            "state text: {}",
            app.state_text()
        );
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
    fn finishing_clears_paused() {
        let mut app = AppState::new();
        app.load_source(BP_PROGRAM).unwrap();
        app.dispatch(":c"); // pause at breakpoint
        assert!(app.paused());
        app.dispatch(":c"); // run to Return
        assert!(!app.paused(), "paused must clear once the program finishes");
    }

    #[test]
    fn finished_hint_does_not_suggest_stepping() {
        let mut app = AppState::new();
        app.load_source(BP_PROGRAM).unwrap();
        app.dispatch(":c"); // pause at breakpoint
        app.dispatch(":c"); // run to Return
        assert!(
            !app.hint().contains("step") && !app.hint().contains("continue"),
            "finished-program hint should not suggest stepping, got: {}",
            app.hint()
        );
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
    fn step_keeps_paused_and_stepping_hint() {
        let mut app = AppState::new();
        app.load_source(BP_PROGRAM).unwrap();
        app.dispatch(""); // empty line == step
        assert!(app.paused(), "manual step should leave the debugger paused");
        assert!(
            app.hint().contains("step"),
            "hint should still advertise stepping, got: {}",
            app.hint()
        );
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

    #[test]
    fn device_after_a_program_then_reset_resets_the_device() {
        // Switching to a REPL device must forget the previously loaded program,
        // so `:reset` resets the device rather than resurrecting the program.
        let mut app = AppState::new();
        app.load_source(BP_PROGRAM).unwrap();
        assert!(app.has_program());
        app.dispatch("device 2");
        assert!(!app.has_program(), "device should switch to a REPL session");
        app.dispatch(":reset");
        assert!(
            !app.has_program(),
            "reset must reset the device, not reload the old program"
        );
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

    #[test]
    fn help_command_toggles_the_overlay() {
        let mut app = AppState::new();
        assert!(!app.show_help());
        app.dispatch(":help");
        assert!(app.show_help(), ":help should open the overlay");
        app.dispatch(":help");
        assert!(!app.show_help(), ":help again should close it");
    }

    #[test]
    fn help_overlay_renders_the_command_reference() {
        use ratatui::Terminal;
        use ratatui::backend::TestBackend;

        let mut app = AppState::new();
        app.dispatch(":help");
        let mut terminal = Terminal::new(TestBackend::new(90, 30)).unwrap();
        terminal.draw(|f| app.render(f)).unwrap();

        let content: String = terminal
            .backend()
            .buffer()
            .content
            .iter()
            .map(|c| c.symbol())
            .collect();
        assert!(content.contains("Help"), "overlay title missing");
        assert!(
            content.contains("cnot"),
            "gate reference missing from overlay"
        );
    }

    fn buffer_lines(buffer: &ratatui::buffer::Buffer, width: u16) -> Vec<String> {
        buffer
            .content
            .chunks(width as usize)
            .map(|row| row.iter().map(|cell| cell.symbol()).collect())
            .collect()
    }

    fn row_containing(lines: &[String], needle: &str) -> Option<usize> {
        lines.iter().position(|line| line.contains(needle))
    }
}
