# Design: `ppvm` TUI — a composable ratatui debugger + REPL for `ppvm-cli`

- **Date:** 2026-07-01
- **Status:** Approved (design); pending implementation plan
- **Branch:** `david/cli-tui`
- **Scope:** Add a terminal UI to the `ppvm` command. The bare `ppvm`
  (no subcommand) launches the TUI. The TUI unifies two capabilities on one
  screen: (1) **step controls for debugging** a loaded `.sst`/`.ssb` program,
  pausing at breakpoints; and (2) a **small interactive REPL** to initialize a
  tableau and apply circuit ops. State display shows the measurement record and
  the tableau. The TUI is built to be **composable** with the ratatui TUI in
  `~/git/stellarscope/` — i.e. its view components must be embeddable into
  another ratatui app later. The existing non-interactive subcommands
  (`run`/`debug`/`parse`/`dump`) are left untouched.

## 1. Motivation

`ppvm-cli` today is a batch tool: `run` executes a program and prints the
measurement record; `debug` offers a line-oriented step loop over stdin.
A bare `ppvm` used to drop into a rustyline REPL, since removed (commit
`8beb96ac`, "Remove the REPL"). Both the removed REPL and the surviving
line-oriented `debug` loop (`commands::debug` / `debug_loop`) are proven
blueprints — they already lower gate commands, step the engine, honor
breakpoints, and render measurements. This design re-homes that behavior into
an interactive, panelled TUI so a user can watch tableau state, the measurement
record, and program position update live while stepping or poking at a device.

The engine already exposes everything the TUI needs as public API on
`ppvm_vihaco::composite::PPVM` — no engine changes are required for v1.

## 2. Locked decisions

| Decision | Choice |
| --- | --- |
| Screen model | **Unified single screen**: persistent panels + one command line. Not tabs, not separate binaries. |
| Program loading | `ppvm <file>` launches the TUI with the program loaded and **paused at pc 0**; bare `ppvm` starts an empty REPL session; `:load <path>` (re)loads from inside the TUI. |
| Existing subcommands | `run` / `debug` / `parse` / `dump` unchanged (still scriptable + tested). |
| Crate boundary | **New library crate `crates/ppvm-tui`** holds app state + `Widget` components; `ppvm-cli` (bin) owns only terminal setup + the event loop. |
| TUI deps | `ratatui = "0.29.0"`, `crossterm = "0.29.0"` — pinned to match stellarscope so panels can dock together. |
| Program panel | A **local** `CodeView<T: Display>`-alike, mirroring stellarscope's shape/`Widget` signature. **Not** a dependency on `stellarscope-inspect` (see §9). |
| Tableau display | Rendered via the already-public `PPVM::state_string()`. `GeneralizedTableau`'s `Display` is left as-is for v1. |
| PauliSum backend | **Deferred.** REPL `device N` creates a Tableau-backed machine only. Loaded programs that declare other backends still run and render. |

## 3. Architecture overview

Two crates, following stellarscope's own lib(inspect/optics)+bin(loop) split:

```
crates/ppvm-tui/            # NEW lib crate — reusable, terminal-agnostic
  src/lib.rs
  src/app.rs               # AppState: owns a PPVM + UI state; handle_key; dispatch
  src/command.rs           # command grammar: parse a line -> Command; gate_spec map
  src/codeview.rs          # local CodeView<T> ring-buffer + cursor (data)
  src/widgets/mod.rs
  src/widgets/program.rs   # impl Widget for &ProgramView   (the CodeView render)
  src/widgets/state.rs     # impl Widget for &StateView     (tableau via state_string)
  src/widgets/record.rs    # impl Widget for &RecordView    (measurement record)
  src/widgets/command.rs   # impl Widget for &CommandLine   (prompt + hint + status)

crates/ppvm-cli/           # bin — owns terminal + loop only
  src/main.rs              # `command: Option<Commands>`; None => tui::run()
  src/tui.rs               # terminal setup (raw mode + alt screen, RAII guard) + event loop
  src/commands.rs          # unchanged
```

**Separation of concerns.** `ppvm-tui` never touches the terminal, never runs a
loop, and never blocks on input. It exposes:

- `AppState` — holds a `PPVM`, the command buffer, a status/error string, the
  program `CodeView`, a REPL scrollback `CodeView`, a `paused: bool`, and scroll
  offsets. Constructed empty (`AppState::new()`) or from a file
  (`AppState::from_file(path)`).
- `AppState::handle_key(&mut self, key: KeyEvent) -> bool` — apply one key
  event; returns whether it was consumed. Pure w.r.t. the terminal.
- `AppState::dispatch(&mut self, line: &str)` — run one command-line string
  (used by Enter and directly by tests).
- The four `Widget` impls on `&…View` newtypes that borrow from `AppState` —
  these are the units a **host app embeds** (render whichever it wants, where it
  wants).
- `AppState::render(&self, frame: &mut Frame)` — a **convenience** full-screen
  composer that lays out all four panels (the standalone layout in §4). A host
  like stellarscope ignores this and lays out the individual `…View` widgets
  itself; `ppvm-cli` just calls it.

`ppvm-cli`'s `tui::run()` owns a `Terminal<CrosstermBackend<Stdout>>`, enables
raw mode + the alternate screen behind an RAII guard that restores them on drop
(even on panic — mirroring stellarscope's `TerminalGuard`), then loops:
poll → `handle_key` → `terminal.draw(|f| app.render(f))` → break on
`app.should_exit`.

## 4. Layout & panels

One screen, four regions (ratatui `Layout` with `Constraint`s):

```
┌ Program ───────────┐┌ State ─────────────────────┐
│  h 0               ││ Generalized Tableau (2q):  │
│▶ cnot 0 1          ││  Destabilizers: [ ... ]    │   ← state_string()
│  measure 0         ││  Stabilizers:   [ ... ]    │
└────────────────────┘└────────────────────────────┘
┌ Measurement record ─────────────────────────────── ┐
│ 0 1                                                 │
└──────────────────────────────────────────────────── ┘
 ppvm> cnot 0 1                    Enter=step  :c  :q
```

- **Program (left), contextual.** When a program is loaded: the compiled
  module's instruction listing (each instruction via `Display`), with a `▶`
  marker at `current_pc()`, auto-scrolled to keep the pc visible. In a pure REPL
  session (no program): a scrollback of entered commands and their inline
  results (`h 0`, `measure 0  => 1`). Both are the same local `CodeView` widget
  with different content and cursor semantics.
- **State (right).** `machine.state_string()` verbatim. This already covers all
  backends and includes coefficients + per-qubit loss. Scrollable if it
  overflows (basic offset; fancy scrolling deferred).
- **Measurement record (bottom band).** `machine.measurement_record()` rendered
  with the established flat convention: `Zero→0`, `One→1`, `Lost→2`, grouped per
  measurement event.
- **Command line (footer).** `ppvm> ` + the current input buffer, a contextual
  hint, and the latest status/error message.

## 5. Command grammar & key handling

The command line is always live for text entry, so single-letter hotkeys would
collide with gate names (`s` = the S gate, etc.). The grammar is therefore
**prefix-disambiguated** and conflict-free:

- **Bare tokens = REPL gate ops.** Ported from the removed REPL's `gate_spec`
  map (name → `CircuitInstruction` + qubit/float arity):
  `device N`; `x y z h s sadj sqrtx sqrty sqrtxadj sqrtyadj t tadj reset measure <q>`;
  `cnot <c> <t>`; `cz <a> <b>`; `rx ry rz <q> <θ>`; `r <q> <axis> <θ>`;
  `rxx ryy rzz <a> <b> <θ>`; `u3 <q> <θ> <φ> <λ>`; `depolarize loss <q> <p>`;
  `depolarize2 <a> <b> <p>`; `paulierror <q> <px> <py> <pz>`;
  `correlatedloss <a> <b> <p0> <p1> <p2>`. Applied via
  `PPVM::apply_circuit_instruction` (already bounds-checks qubit indices).
  New measurement outcomes echo inline as `=> <bits>`.
- **`:`-prefixed = meta / debug commands.** `:load <path>`, `:continue` / `:c`,
  `:step` / `:s`, `:reset`, `:quit` / `:q`.
- **Empty line + Enter = step** when paused. Matches the removed `debug` loop's
  bare-Enter-steps ergonomics — fast repeated stepping without reaching for a
  prefix.
- **Footer hint** is contextual: `Enter=step  :c=continue  :q=quit` when a
  program is loaded and paused; `ppvm>  (type a gate, or :load <file>)`
  otherwise.

`handle_key` maps: printable chars → push to buffer; Backspace → pop; Enter →
`dispatch(buffer)` then clear; Ctrl-C / Esc on an empty buffer → quit (Ctrl-C on
a non-empty buffer clears it, shell-like). Up/Down reserved for scrollback
(basic history/scroll; full history optional).

Stepping semantics reuse the `debug_loop` logic: `step_once()` returns a
`StepOutcome`; `Breakpoint` sets `paused = true` and shows `-- breakpoint hit --`;
`Return`/`Halt` shows "Program finished." `:continue` steps in a tight loop
until the next `Breakpoint` or program end (the loop yields between engine steps
so the UI can repaint).

## 6. Engine integration

All public on `PPVM` today — **no engine changes for v1**:

- REPL device: `PPVM::with_qubits(n)` (Tableau-backed) for `device N`.
- Apply gate: `apply_circuit_instruction(inst, &qubits, &params)`.
- Step: `step_once() -> StepOutcome` (`Continue`/`Breakpoint`/`Return`/`Halt`,
  re-exported from `ppvm_vihaco::composite`).
- Position: `current_pc()`, `current_instruction()`.
- Records: `measurement_record() -> Vec<MeasurementResult>`,
  `trace_record() -> Vec<f64>`.
- State text: `state_string() -> String`.

**Program loading** builds the code listing without any private accessor: load
the module via the public `ppvm_vihaco::load_module_file(path)`, format its
public `code` field (each `PPVMInstruction` implements `Display`) into the
`CodeView`, then `machine.load(&module)` + `machine.init()`, leaving the machine
paused at pc 0. `:load` reuses this exact path. (If a program listing ever needs
richer data than `load_module_file` exposes, that is an additive engine change,
out of v1 scope.)

## 7. Error handling

Non-fatal, exactly like the removed REPL and the `debug` loop: a bad gate name,
wrong arity, out-of-range qubit, missing device (`device N` not yet run), or a
parse error on `:load` is caught and written to the status line; the loop
continues. `eyre::Result` throughout the library boundary. The terminal is
always restored via the RAII guard, including on panic, so a crash never leaves
the user's terminal in raw mode.

## 8. Testing

Mirror stellarscope's input tests and the old `repl_loop` / `debug_loop` tests —
**no terminal required**:

- **Command grammar** (`command.rs`): parsing a line into a `Command`, the
  `gate_spec` arity table, `:`-prefix routing, empty-line-steps.
- **State transitions** (`app.rs`): construct an `AppState`, feed `KeyEvent`s or
  call `dispatch`, assert engine effects — `device 1` then `x 0` then `measure 0`
  yields `=> 1`; `device 2; x 0; cnot 0 1; measure 0; measure 1` → `1,1`;
  loading `BREAKPOINT_PROGRAM` and stepping pauses at the breakpoint;
  out-of-range qubit surfaces an error and keeps looping; `:quit` sets
  `should_exit`.
- **Rendering** (optional, light): a smoke test with ratatui's `TestBackend`
  that a panelled frame renders without panic and contains expected substrings
  (e.g. the `▶` pc marker, a measurement bit).

Workspace gates unchanged: `cargo test --workspace`, `cargo fmt` before commit.

## 9. Composability contract (embedding into stellarscope later)

The whole point of the lib/bin split. Guarantees `ppvm-tui` upholds:

1. **Matching deps** — `ratatui 0.29` + `crossterm 0.29`, identical to
   stellarscope, so both can share one `Terminal` and one event stream.
2. **Terminal-agnostic components** — every view is `impl Widget for &SomeView`
   with the stock `render(self, area: Rect, buf: &mut Buffer)` signature; none
   own the terminal or run a loop.
3. **Poll-friendly state** — `AppState::handle_key(&mut self, KeyEvent) -> bool`
   and `dispatch(&str)` let a host app forward events and drive state without
   `ppvm-tui` blocking.

A future stellarscope integration then looks like: add `ppvm-tui` as a dep, hold
a ppvm `AppState`, `frame.render_widget(&StateView(&app), area)` into its layout,
and forward relevant key events to `app.handle_key`. Nothing in `ppvm-tui` needs
to change for that.

**Why not reuse stellarscope's `CodeView` directly?** It lives in
`stellarscope-inspect`, an unpublished crate. Depending on it would (a) make
ppvm depend on stellarscope by filesystem path (backwards, non-publishable), and
(b) drag in stellarscope's local `vihaco 0.1.0` path dep + `stellarscope-fpga`,
clashing with ppvm's registry `vihaco 0.1.1`. The widget is ~50 trivial lines,
so we reimplement a same-shaped `CodeView<T: Display>` locally. If real sharing
is wanted later, the clean move is extracting the dependency-light
`CodeView`/`Window` widgets into a **third** shared crate both repos depend on —
not either repo depending on the other.

## 10. Scope / deferred (v1)

Deferred by explicit request or to keep v1 tight:

- **PauliSum / LossyPauliSum in the REPL** — no `observable` / `trace` setup;
  `device N` is Tableau-only. (Loaded `.sst` programs may still declare any
  backend; `state_string()` renders them.)
- **Mid-program ad-hoc gate injection** — applying a REPL gate while stepping a
  loaded program is not specially designed (the pc/appended-code interaction of
  `execute_single_instruction` is a footgun). Gate ops are intended for REPL
  sessions in v1.
- **Richer structured tableau rendering** and any `GeneralizedTableau::Display`
  overhaul — use `state_string()` for now.
- **The actual stellarscope integration** — we ship only the embeddable surface
  (§9), not the cross-repo wiring.
- Mouse, theming/config, and full readline-style history/editing beyond basic
  buffer edit + scrollback.

## 11. Success criteria

- `ppvm` (no args) opens the TUI; `q`/`:q`/Ctrl-C exits cleanly with the
  terminal restored.
- REPL: `device 2`, apply gates, `measure` — outcomes echo inline and the
  measurement-record panel updates; the State panel reflects the tableau.
- Debug: `ppvm prog.sst` opens paused at pc 0; Enter steps (Program panel `▶`
  advances, records update); `:c` runs to the next breakpoint / end; an authored
  `breakpoint` pauses.
- `:load other.sst` swaps the program without relaunching.
- Existing `run`/`debug`/`parse`/`dump` behavior and tests are unchanged.
- `cargo test --workspace` passes; `ppvm-tui` has no dependency on a terminal in
  its unit tests.
