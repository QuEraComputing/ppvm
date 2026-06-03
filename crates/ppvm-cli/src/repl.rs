// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Interactive REPL — a quantum-circuit playground. Allocate a fixed-size
//! device with `device N`, then apply gates and measurements one line at a time
//! against a single persistent machine, seeing measurement outcomes inline.

use eyre::{Result, WrapErr, bail, eyre};
use ppvm_vihaco::CircuitInstruction;
use ppvm_vihaco::composite::PPVM;
use ppvm_vihaco::measurements::MeasurementResult;
use rustyline::DefaultEditor;
use rustyline::error::ReadlineError;
#[cfg(test)]
use std::io::BufRead;
use std::io::Write;
use std::path::PathBuf;

/// How a gate command lowers: the engine instruction plus how many qubit and
/// float operands it consumes (qubits first, then floats — the order
/// `apply_circuit_instruction` expects).
struct GateSpec {
    inst: CircuitInstruction,
    qubits: usize,
    floats: usize,
}

/// Resolve a gate command name to its spec, or `None` if it is not a gate.
/// `TwoQubitPauliError` is intentionally absent — its tableau arm is `todo!()`.
fn gate_spec(name: &str) -> Option<GateSpec> {
    use CircuitInstruction::*;
    let (inst, qubits, floats) = match name {
        "x" => (X, 1, 0),
        "y" => (Y, 1, 0),
        "z" => (Z, 1, 0),
        "h" => (H, 1, 0),
        "s" => (S, 1, 0),
        "sadj" => (SAdj, 1, 0),
        "sqrtx" => (SqrtX, 1, 0),
        "sqrty" => (SqrtY, 1, 0),
        "sqrtxadj" => (SqrtXAdj, 1, 0),
        "sqrtyadj" => (SqrtYAdj, 1, 0),
        "t" => (T, 1, 0),
        "tadj" => (TAdj, 1, 0),
        "measure" => (Measure, 1, 0),
        "reset" => (Reset, 1, 0),
        "cnot" => (CNOT, 2, 0),
        "cz" => (CZ, 2, 0),
        "rx" => (RX, 1, 1),
        "ry" => (RY, 1, 1),
        "rz" => (RZ, 1, 1),
        "r" => (R, 1, 2),
        "rxx" => (RXX, 2, 1),
        "ryy" => (RYY, 2, 1),
        "rzz" => (RZZ, 2, 1),
        "u3" => (U3, 1, 3),
        "depolarize" => (Depolarize, 1, 1),
        "depolarize2" => (Depolarize2, 2, 1),
        "loss" => (Loss, 1, 1),
        "paulierror" => (PauliError, 1, 3),
        "correlatedloss" => (CorrelatedLoss, 2, 3),
        _ => return None,
    };
    Some(GateSpec {
        inst,
        qubits,
        floats,
    })
}

/// Whether the loop should continue prompting or exit.
enum Outcome {
    Continue,
    Quit,
}

/// Launch the interactive REPL with line editing and command history.
/// History recall (up/down arrows), cursor movement, and Ctrl-R search come
/// from rustyline; the per-command logic is the same `dispatch` the scripted
/// tests drive through `repl_loop`.
pub fn repl() -> Result<()> {
    let mut rl = DefaultEditor::new()?;
    let history = history_path();
    if let Some(path) = &history {
        let _ = rl.load_history(path); // best-effort: a missing file is fine
    }

    let mut machine: Option<PPVM> = None;
    let mut output = std::io::stdout();
    loop {
        match rl.readline("ppvm> ") {
            Ok(line) => {
                let trimmed = line.trim();
                if trimmed.is_empty() {
                    continue;
                }
                let _ = rl.add_history_entry(trimmed);
                match dispatch(trimmed, &mut machine, &mut output) {
                    Ok(Outcome::Continue) => {}
                    Ok(Outcome::Quit) => break,
                    Err(e) => writeln!(output, "error: {e}")?,
                }
            }
            // Ctrl-C: abandon the current line, keep the session (shell-like).
            Err(ReadlineError::Interrupted) => continue,
            // Ctrl-D on an empty line: leave cleanly.
            Err(ReadlineError::Eof) => break,
            Err(e) => return Err(e).wrap_err("readline failed"),
        }
    }

    if let Some(path) = &history {
        let _ = rl.save_history(path); // best-effort: don't fail the session
    }
    Ok(())
}

/// Where to persist history across sessions: `$HOME/.ppvm_history`. `None`
/// (so history is session-only) when there's no `HOME` — e.g. on Windows,
/// where you'd reach for the `dirs` crate instead.
fn history_path() -> Option<PathBuf> {
    std::env::var_os("HOME").map(|home| PathBuf::from(home).join(".ppvm_history"))
}

/// Core REPL loop, generic over its IO so tests can drive it with scripted
/// input. Holds a single `Option<PPVM>` — `None` until `device N`. Command-level
/// errors are printed and the loop continues; only `quit`/`exit`/EOF exit.
///
/// Test-only: the interactive entry point is `repl`, which runs the same
/// `dispatch` under a rustyline editor for history and line editing.
#[cfg(test)]
fn repl_loop(input: &mut impl BufRead, output: &mut impl Write) -> Result<()> {
    let mut machine: Option<PPVM> = None;
    loop {
        write!(output, "ppvm> ")?;
        output.flush()?;

        let mut line = String::new();
        if input.read_line(&mut line)? == 0 {
            // EOF: leave cleanly, ending the dangling prompt line.
            writeln!(output)?;
            break;
        }
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        match dispatch(line, &mut machine, output) {
            Ok(Outcome::Continue) => {}
            Ok(Outcome::Quit) => break,
            Err(e) => writeln!(output, "error: {e}")?,
        }
    }
    Ok(())
}

/// Parse and run one command line.
fn dispatch(line: &str, machine: &mut Option<PPVM>, output: &mut impl Write) -> Result<Outcome> {
    let mut tokens = line.split_whitespace();
    let cmd = tokens.next().expect("line is non-empty after trim");
    let args: Vec<&str> = tokens.collect();

    match cmd {
        "quit" | "exit" => return Ok(Outcome::Quit),
        "help" => print_help(output)?,
        "device" => cmd_device(&args, machine, output)?,
        "show" => cmd_show(machine, output)?,
        _ => cmd_gate(cmd, &args, machine, output)?,
    }
    Ok(Outcome::Continue)
}

/// `device N` — (re)create a fresh N-qubit device, discarding any prior state.
fn cmd_device(args: &[&str], machine: &mut Option<PPVM>, output: &mut impl Write) -> Result<()> {
    let [n] = args else {
        bail!("usage: device N");
    };
    let n: usize = n
        .parse()
        .wrap_err_with(|| format!("invalid qubit count {n:?}"))?;
    if n == 0 {
        bail!("device must have at least 1 qubit");
    }

    let existed = machine.is_some();
    *machine = Some(PPVM::with_qubits(n)?);
    if existed {
        writeln!(
            output,
            "ok: fresh {n}-qubit device (previous state discarded)"
        )?;
    } else {
        writeln!(output, "ok: fresh {n}-qubit device")?;
    }
    Ok(())
}

/// `show` — print the current tableau / Pauli state.
fn cmd_show(machine: &mut Option<PPVM>, output: &mut impl Write) -> Result<()> {
    let machine = require_device(machine)?;
    writeln!(output, "{}", machine.state_string())?;
    Ok(())
}

/// A gate command: `<name> <qubit…> [param…]`. Applies the gate and, if it
/// produced any measurement outcomes, prints them as `=> <bits>`.
fn cmd_gate(
    name: &str,
    args: &[&str],
    machine: &mut Option<PPVM>,
    output: &mut impl Write,
) -> Result<()> {
    let spec = gate_spec(name).ok_or_else(|| eyre!("unknown command {name:?}; try \"help\""))?;
    let machine = require_device(machine)?;

    let expected = spec.qubits + spec.floats;
    if args.len() != expected {
        bail!(
            "{name} takes {} qubit(s) and {} param(s), got {} argument(s)",
            spec.qubits,
            spec.floats,
            args.len()
        );
    }

    let (qubit_args, param_args) = args.split_at(spec.qubits);
    let qubits = qubit_args
        .iter()
        .map(|t| {
            t.parse::<usize>()
                .wrap_err_with(|| format!("invalid qubit index {t:?}"))
        })
        .collect::<Result<Vec<_>>>()?;
    let params = param_args
        .iter()
        .map(|t| {
            t.parse::<f64>()
                .wrap_err_with(|| format!("invalid parameter {t:?}"))
        })
        .collect::<Result<Vec<_>>>()?;

    // New entries in the measurement record are this command's outcomes.
    let before = machine.measurement_record().len();
    machine.apply_circuit_instruction(spec.inst, &qubits, &params)?;
    let record = machine.measurement_record();
    if record.len() > before {
        writeln!(output, "=> {}", format_outcomes(&record[before..]))?;
    }
    Ok(())
}

/// Borrow the device, or error if none has been allocated yet.
fn require_device(machine: &mut Option<PPVM>) -> Result<&mut PPVM> {
    machine
        .as_mut()
        .ok_or_else(|| eyre!("no device; run \"device N\" first"))
}

/// Render measurement outcomes as a flat bit string: `0`/`1`, lost qubit = `2`
/// (the same convention as the `run` command's output).
fn format_outcomes(records: &[MeasurementResult]) -> String {
    records
        .iter()
        .flatten()
        .map(|outcome| char::from(b'0' + *outcome as u8))
        .collect()
}

fn print_help(output: &mut impl Write) -> Result<()> {
    writeln!(output, "meta:")?;
    writeln!(
        output,
        "  device N            (re)create a fresh N-qubit device"
    )?;
    writeln!(output, "  show                print the current state")?;
    writeln!(output, "  help                show this help")?;
    writeln!(output, "  quit | exit | EOF   leave the REPL")?;
    writeln!(
        output,
        "gates (<q> = qubit index, angles/probs are floats):"
    )?;
    writeln!(
        output,
        "  x y z h s sadj sqrtx sqrty sqrtxadj sqrtyadj t tadj reset measure <q>"
    )?;
    writeln!(output, "  cnot <c> <t> | cz <q0> <q1>")?;
    writeln!(output, "  rx ry rz <q> <angle> | r <q> <axis> <angle>")?;
    writeln!(output, "  u3 <q> <theta> <phi> <lam>")?;
    writeln!(output, "  rxx ryy rzz <q0> <q1> <angle>")?;
    writeln!(
        output,
        "  depolarize loss <q> <p> | depolarize2 <q0> <q1> <p>"
    )?;
    writeln!(
        output,
        "  paulierror <q> <px> <py> <pz> | correlatedloss <q0> <q1> <p0> <p1> <p2>"
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Drive the loop with a scripted session, returning everything it wrote.
    fn session(script: &str) -> String {
        let mut input = std::io::Cursor::new(script.as_bytes().to_vec());
        let mut output: Vec<u8> = Vec::new();
        repl_loop(&mut input, &mut output).unwrap();
        String::from_utf8(output).unwrap()
    }

    #[test]
    fn device_then_x_then_measure_is_one() {
        let out = session("device 1\nx 0\nmeasure 0\nquit\n");
        assert!(out.contains("ok: fresh 1-qubit device"), "{out}");
        assert!(out.contains("=> 1"), "X|0> measured should be 1:\n{out}");
    }

    #[test]
    fn fresh_measure_is_zero() {
        let out = session("device 1\nmeasure 0\nquit\n");
        assert!(out.contains("=> 0"), "|0> measured should be 0:\n{out}");
    }

    #[test]
    fn gate_before_device_errors_and_continues() {
        // The error is printed and the loop keeps going: the later device+measure
        // still works.
        let out = session("x 0\ndevice 1\nmeasure 0\nquit\n");
        assert!(out.contains("no device"), "{out}");
        assert!(
            out.contains("=> 0"),
            "loop should continue after the error:\n{out}"
        );
    }

    #[test]
    fn device_twice_reports_discarded_state() {
        let out = session("device 1\ndevice 2\nquit\n");
        assert!(out.contains("previous state discarded"), "{out}");
    }

    #[test]
    fn show_renders_the_state() {
        let out = session("device 1\nshow\nquit\n");
        let expected = PPVM::with_qubits(1).unwrap().state_string();
        assert!(out.contains(expected.trim()), "show output missing:\n{out}");
    }

    #[test]
    fn unknown_command_errors_and_continues() {
        let out = session("device 1\nbogus 0\nmeasure 0\nquit\n");
        assert!(out.contains("unknown command"), "{out}");
        assert!(out.contains("=> 0"), "loop should continue:\n{out}");
    }

    #[test]
    fn eof_exits_cleanly() {
        // No quit line: the session ends at EOF without hanging or erroring.
        let out = session("device 1\n");
        assert!(out.contains("ok: fresh 1-qubit device"), "{out}");
    }

    #[test]
    fn cnot_respects_control_and_target_order() {
        // x 0 -> |10>; cnot 0 1 (control q0=1) flips q1 -> |11>. Both measure 1.
        // If control/target were swapped, q1 would stay 0.
        let out = session("device 2\nx 0\ncnot 0 1\nmeasure 0\nmeasure 1\nquit\n");
        // Each scripted line follows a "ppvm> " prompt with no echo, so a result
        // line reads "ppvm> => 1"; take the text after "=> ".
        let measurements: Vec<&str> = out.lines().filter_map(|l| l.split("=> ").nth(1)).collect();
        assert_eq!(measurements, vec!["1", "1"], "{out}");
    }

    #[test]
    fn two_qubit_float_gate_runs() {
        // rxx (2 qubits + 1 float) should parse and apply without error.
        let out = session("device 2\nrxx 0 1 0.5\nquit\n");
        assert!(!out.contains("error"), "rxx should run cleanly:\n{out}");
    }

    #[test]
    fn bad_arity_errors() {
        let out = session("device 1\nx\nquit\n");
        assert!(out.contains("takes"), "arity error expected:\n{out}");
    }

    #[test]
    fn out_of_range_qubit_errors_not_panics() {
        let out = session("device 1\nx 3\nquit\n");
        assert!(out.contains("out of range"), "{out}");
    }
}
