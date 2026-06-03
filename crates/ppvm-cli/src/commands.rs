// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use eyre::{Result, WrapErr};
use ppvm_vihaco::composite::{PPVM, StepOutcome};
use ppvm_vihaco::measurements::MeasurementResult;
use std::io::{BufRead, Write};
use std::path::Path;

/// Output format for `parse`.
#[derive(Clone, Debug, clap::ValueEnum)]
pub enum Format {
    Pretty,
    Debug,
    Json,
}

/// Output format for the measurement record from `run`.
#[derive(Clone, Debug, clap::ValueEnum)]
pub enum MeasurementFormat {
    /// One flat bit string per shot: `0`/`1`, lost qubit = `2`.
    Bits,
    /// Raw debug representation of all shots.
    Debug,
}

pub fn run(
    file: &str,
    shots: usize,
    seed: Option<u64>,
    output: Option<&str>,
    quiet: bool,
    format: MeasurementFormat,
) -> Result<()> {
    // Compile once, then run every shot against the shared module. The thread
    // pool is sized once in `main` via the top-level `--threads` flag.
    let module =
        ppvm_vihaco::load_module_file(file).wrap_err_with(|| format!("failed to load {file}"))?;
    let records = ppvm_vihaco::shots::run_shots(&module, shots, seed)
        .wrap_err_with(|| format!("failed to run {file}"))?;
    if quiet {
        return Ok(());
    }

    let text = match format {
        MeasurementFormat::Bits => format_shots(&records),
        MeasurementFormat::Debug => format!("{records:?}"),
    };

    match output {
        Some(path) => {
            std::fs::write(path, format!("{text}\n"))
                .wrap_err_with(|| format!("failed to write {path}"))?;
            eprintln!("Results written to {path}");
        }
        None => println!("{text}"),
    }
    Ok(())
}

/// Render one shot per line, each as a flat bit string.
fn format_shots(records: &[Vec<MeasurementResult>]) -> String {
    records
        .iter()
        .map(|shot| format_shot(shot))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Render a shot's full measurement record as one flat bit string, all events
/// and qubits concatenated: `Zero` → `0`, `One` → `1`, `Lost` → `2` (the
/// outcome's own enum value). An empty record renders as the empty string.
fn format_shot(record: &[MeasurementResult]) -> String {
    record
        .iter()
        .flatten()
        .map(|outcome| char::from(b'0' + *outcome as u8))
        .collect()
}

pub fn parse(file: &str, format: Format) -> Result<()> {
    let source =
        std::fs::read_to_string(file).wrap_err_with(|| format!("failed to read {file}"))?;
    let parsed = ppvm_vihaco::parse_program(&source)?;

    match format {
        Format::Json => {
            eprintln!("Warning: JSON format not yet supported for AST, using debug format");
            println!("{:#?}", parsed);
        }
        Format::Debug => {
            println!("{:#?}", parsed);
        }
        Format::Pretty => {
            println!("Module:");
            println!("  Headers: {}", parsed.headers.len());
            for (i, header) in parsed.headers.iter().enumerate() {
                println!("    [{}] {:?}", i, header);
            }
            println!("  Functions: {}", parsed.functions.len());
            for (i, func) in parsed.functions.iter().enumerate() {
                println!(
                    "    [{}] {}({} params, {} body items)",
                    i,
                    func.name,
                    func.params.len(),
                    func.body.len()
                );
            }
        }
    }

    Ok(())
}

pub fn dump(file: &str, output: Option<&str>, force: bool) -> Result<()> {
    let output_file = match output {
        Some(output_file_name) => output_file_name.to_string(),
        None => Path::new(file)
            .with_extension("ssb")
            .to_string_lossy()
            .into_owned(),
    };

    // Don't clobber an existing file unless asked to.
    if !force && Path::new(&output_file).exists() {
        return Err(eyre::eyre!(
            "{output_file} already exists; pass --force to overwrite"
        ));
    }

    ppvm_vihaco::dump_file(file, &output_file)
        .wrap_err_with(|| format!("failed to dump {file}"))?;
    eprintln!("Bytecode written to {output_file}");
    Ok(())
}

/// A command entered at the debugger prompt.
enum DebugCommand {
    Step,
    Continue,
    Quit,
}

/// Step through a program interactively, pausing at `breakpoint` instructions.
/// With `break_at_start`, also pauses before the first instruction so any
/// program can be stepped from the beginning.
pub fn debug(file: &str, break_at_start: bool) -> Result<()> {
    let stdin = std::io::stdin();
    let mut input = stdin.lock();
    let mut output = std::io::stdout();
    debug_loop(file, break_at_start, &mut input, &mut output)
}

/// Core debugger loop, generic over its IO so it can be driven by tests.
fn debug_loop(
    file: &str,
    break_at_start: bool,
    input: &mut impl BufRead,
    output: &mut impl Write,
) -> Result<()> {
    let mut machine = PPVM::default();
    machine
        .load_file(file)
        .wrap_err_with(|| format!("failed to load {file}"))?;
    machine.init()?;

    let mut paused = break_at_start;
    let mut ever_paused = paused;

    loop {
        // Safety net: stop if execution has run off the end of the code.
        if machine.current_instruction().is_none() {
            writeln!(output, "Program counter past end of code.")?;
            break;
        }

        if paused {
            print_location(&machine, output)?;
            match read_command(input, output)? {
                DebugCommand::Quit => {
                    writeln!(output, "Quit.")?;
                    return Ok(());
                }
                DebugCommand::Continue => paused = false,
                DebugCommand::Step => {}
            }
        }

        match machine.step_once()? {
            StepOutcome::Continue => {}
            StepOutcome::Breakpoint => {
                paused = true;
                ever_paused = true;
                writeln!(output, "-- breakpoint hit --")?;
            }
            StepOutcome::Return | StepOutcome::Halt => {
                writeln!(output, "Program finished.")?;
                break;
            }
        }
    }

    writeln!(
        output,
        "Measurements: {}",
        format_shot(&machine.measurement_record())
    )?;
    if !ever_paused {
        writeln!(
            output,
            "(no breakpoint was hit; pass --break-at-start to step from the beginning)"
        )?;
    }
    Ok(())
}

/// Print the program counter, the next instruction, and measurements so far.
fn print_location(machine: &PPVM, output: &mut impl Write) -> Result<()> {
    let pc = machine.current_pc();
    match machine.current_instruction() {
        Some(inst) => writeln!(output, "pc={pc}  next: {inst}")?,
        None => writeln!(output, "pc={pc}  (end of code)")?,
    }
    writeln!(
        output,
        "measurements: {}",
        format_shot(&machine.measurement_record())
    )?;
    Ok(())
}

/// Prompt for and read a debugger command. A bare Enter steps; EOF quits.
fn read_command(input: &mut impl BufRead, output: &mut impl Write) -> Result<DebugCommand> {
    loop {
        write!(output, "> s step | c continue | q quit: ")?;
        output.flush()?;

        let mut line = String::new();
        if input.read_line(&mut line)? == 0 {
            // EOF (e.g. stdin closed): treat as quit so we never spin.
            return Ok(DebugCommand::Quit);
        }
        match line.trim() {
            "" | "s" | "step" => return Ok(DebugCommand::Step),
            "c" | "continue" => return Ok(DebugCommand::Continue),
            "q" | "quit" => return Ok(DebugCommand::Quit),
            other => writeln!(output, "unknown command: {other:?}")?,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_vihaco::measurements::MeasurementOutcome;
    use std::fs;

    /// Minimal program that compiles and measures q0 in |0> (deterministic).
    const PROGRAM: &str =
        "device circuit.n_qubits 1;\nfn @main() { const.u64 0\n gate measure\n ret }\n";

    fn row(outcomes: &[MeasurementOutcome]) -> MeasurementResult {
        outcomes.iter().copied().collect()
    }

    /// Write `contents` to a uniquely-named temp file and return its path.
    fn temp_file(name: &str, contents: &str) -> String {
        let path = std::env::temp_dir().join(name);
        fs::write(&path, contents).unwrap();
        path.to_string_lossy().into_owned()
    }

    // ─── format_shot ─────────────────────────────────────────────────────

    #[test]
    fn format_shot_empty_record_is_empty() {
        assert_eq!(format_shot(&[]), "");
    }

    #[test]
    fn format_shot_concatenates_qubits_within_an_event() {
        let record = vec![row(&[
            MeasurementOutcome::One,
            MeasurementOutcome::Zero,
            MeasurementOutcome::One,
        ])];
        assert_eq!(format_shot(&record), "101");
    }

    #[test]
    fn format_shot_flattens_events_with_no_separator() {
        let record = vec![
            row(&[MeasurementOutcome::One]),
            row(&[MeasurementOutcome::Zero]),
        ];
        assert_eq!(format_shot(&record), "10");
    }

    #[test]
    fn format_shot_renders_lost_qubit_as_two() {
        let record = vec![
            row(&[MeasurementOutcome::One, MeasurementOutcome::Lost]),
            row(&[MeasurementOutcome::Zero]),
        ];
        assert_eq!(format_shot(&record), "120");
    }

    #[test]
    fn format_shots_joins_shots_with_newlines() {
        let shots = vec![
            vec![row(&[MeasurementOutcome::One])],
            vec![row(&[MeasurementOutcome::Zero])],
        ];
        assert_eq!(format_shots(&shots), "1\n0");
    }

    // ─── run ───────────────────────────────────────────────────────────

    #[test]
    fn run_succeeds_on_valid_file() {
        let src = temp_file("ppvm_cli_run_ok.sst", PROGRAM);
        let res = run(&src, 3, None, None, true, MeasurementFormat::Bits);
        let _ = fs::remove_file(&src);
        assert!(res.is_ok(), "got: {res:?}");
    }

    #[test]
    fn run_writes_one_line_per_shot_to_output_file() {
        let src = temp_file("ppvm_cli_run_output.sst", PROGRAM);
        let out = std::env::temp_dir().join("ppvm_cli_run_output.txt");
        let _ = fs::remove_file(&out);

        run(&src, 4, None, out.to_str(), false, MeasurementFormat::Bits).unwrap();
        let contents = fs::read_to_string(&out).unwrap();
        // Four deterministic shots of |0>, one per line.
        assert_eq!(contents, "0\n0\n0\n0\n");

        let _ = fs::remove_file(&src);
        let _ = fs::remove_file(&out);
    }

    #[test]
    fn run_errors_with_context_on_missing_file() {
        let err = run(
            "/no/such/file.sst",
            1,
            None,
            None,
            false,
            MeasurementFormat::Bits,
        )
        .unwrap_err();
        assert!(err.to_string().contains("failed to load"), "got: {err}");
    }

    // ─── parse ─────────────────────────────────────────────────────────

    #[test]
    fn parse_succeeds_on_valid_file() {
        let src = temp_file("ppvm_cli_parse_ok.sst", PROGRAM);
        let res = parse(&src, Format::Debug);
        let _ = fs::remove_file(&src);
        assert!(res.is_ok(), "got: {res:?}");
    }

    #[test]
    fn parse_errors_with_context_on_missing_file() {
        let err = parse("/no/such/file.sst", Format::Pretty).unwrap_err();
        assert!(err.to_string().contains("failed to read"), "got: {err}");
    }

    // ─── dump ──────────────────────────────────────────────────────────

    #[test]
    fn dump_writes_default_ssb_path_when_output_omitted() {
        let src = temp_file("ppvm_cli_dump_default.sst", PROGRAM);
        let expected = Path::new(&src).with_extension("ssb");
        let _ = fs::remove_file(&expected); // clear any leftover from a prior run

        dump(&src, None, false).unwrap();
        assert!(expected.exists(), "default .ssb should have been written");

        let _ = fs::remove_file(&src);
        let _ = fs::remove_file(&expected);
    }

    #[test]
    fn dump_writes_explicit_output_path() {
        let src = temp_file("ppvm_cli_dump_explicit.sst", PROGRAM);
        let out = std::env::temp_dir().join("ppvm_cli_dump_explicit_out.ssb");
        let _ = fs::remove_file(&out);

        dump(&src, Some(out.to_str().unwrap()), false).unwrap();
        assert!(out.exists());

        let _ = fs::remove_file(&src);
        let _ = fs::remove_file(&out);
    }

    #[test]
    fn dump_refuses_to_clobber_without_force() {
        let src = temp_file("ppvm_cli_dump_clobber.sst", PROGRAM);
        let out = std::env::temp_dir().join("ppvm_cli_dump_clobber_out.ssb");
        fs::write(&out, b"existing").unwrap();

        let err = dump(&src, Some(out.to_str().unwrap()), false).unwrap_err();
        assert!(err.to_string().contains("already exists"), "got: {err}");
        // The existing file must be left untouched.
        assert_eq!(fs::read(&out).unwrap(), b"existing");

        let _ = fs::remove_file(&src);
        let _ = fs::remove_file(&out);
    }

    #[test]
    fn dump_overwrites_existing_with_force() {
        let src = temp_file("ppvm_cli_dump_force.sst", PROGRAM);
        let out = std::env::temp_dir().join("ppvm_cli_dump_force_out.ssb");
        fs::write(&out, b"existing").unwrap();

        dump(&src, Some(out.to_str().unwrap()), true).unwrap();
        // Replaced with real bytecode, not the placeholder.
        assert_ne!(fs::read(&out).unwrap(), b"existing");

        let _ = fs::remove_file(&src);
        let _ = fs::remove_file(&out);
    }

    // ─── debug ─────────────────────────────────────────────────────────

    /// Program with a `breakpoint` before measuring q0 in |0> (deterministic).
    const BREAKPOINT_PROGRAM: &str = "device circuit.n_qubits 1;\nfn @main() { breakpoint\n const.u64 0\n gate measure\n ret }\n";

    /// Drive `debug_loop` with scripted input, returning the captured output.
    fn run_debug(program: &str, name: &str, break_at_start: bool, script: &str) -> String {
        let src = temp_file(name, program);
        let mut input = script.as_bytes();
        let mut output: Vec<u8> = Vec::new();
        debug_loop(&src, break_at_start, &mut input, &mut output).unwrap();
        let _ = fs::remove_file(&src);
        String::from_utf8(output).unwrap()
    }

    #[test]
    fn debug_break_at_start_steps_through_to_finish() {
        // PROGRAM is const.u64 0 / gate measure / ret = 3 steps.
        let out = run_debug(PROGRAM, "ppvm_cli_debug_step.sst", true, "s\ns\ns\n");
        assert!(
            out.contains("next: Measure"),
            "should display the gate: {out}"
        );
        assert!(out.contains("Program finished."), "{out}");
        assert!(
            out.contains("Measurements: 0"),
            "q0 in |0> measures 0: {out}"
        );
    }

    #[test]
    fn debug_continue_runs_to_end() {
        let out = run_debug(PROGRAM, "ppvm_cli_debug_continue.sst", true, "c\n");
        assert!(out.contains("Program finished."), "{out}");
        assert!(out.contains("Measurements: 0"), "{out}");
    }

    #[test]
    fn debug_honors_authored_breakpoint() {
        // Not breaking at start: must run until the `breakpoint` pauses it.
        let out = run_debug(BREAKPOINT_PROGRAM, "ppvm_cli_debug_bp.sst", false, "c\n");
        assert!(out.contains("-- breakpoint hit --"), "{out}");
        assert!(out.contains("Program finished."), "{out}");
        // A breakpoint was hit, so no "use --break-at-start" hint.
        assert!(!out.contains("no breakpoint was hit"), "{out}");
    }

    #[test]
    fn debug_quit_stops_before_finishing() {
        let out = run_debug(PROGRAM, "ppvm_cli_debug_quit.sst", true, "q\n");
        assert!(out.contains("Quit."), "{out}");
        assert!(!out.contains("Program finished."), "{out}");
        assert!(
            !out.contains("Measurements:"),
            "quit prints no record: {out}"
        );
    }

    #[test]
    fn debug_without_breakpoint_prints_hint() {
        // No breakpoint, no break-at-start, empty input: runs straight through
        // and tells the user how to step.
        let out = run_debug(PROGRAM, "ppvm_cli_debug_hint.sst", false, "");
        assert!(out.contains("Program finished."), "{out}");
        assert!(out.contains("no breakpoint was hit"), "{out}");
    }
}
