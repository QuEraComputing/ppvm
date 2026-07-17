// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The TUI command grammar. Bare tokens are gate ops; `:`-prefixed tokens are
//! meta/debug commands; an empty line means "step". The gate table is ported
//! from the removed rustyline REPL.

use eyre::{Result, WrapErr, bail, eyre};
use ppvm_vihaco::CircuitInstruction;

/// How a gate command lowers: the engine instruction plus how many qubit and
/// float operands it consumes (qubits first, then floats).
pub struct GateSpec {
    pub inst: CircuitInstruction,
    pub qubits: usize,
    pub floats: usize,
}

/// Resolve a gate name to its spec, or `None` if it is not a gate.
/// `TwoQubitPauliError` is intentionally absent — its tableau arm is `todo!()`.
pub fn gate_spec(name: &str) -> Option<GateSpec> {
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

/// One parsed command-line entry.
#[derive(Debug, PartialEq)]
pub enum Command {
    /// `device N` — (re)create a fresh N-qubit tableau device.
    Device(usize),
    /// A gate op, e.g. `cnot 0 1` or `rx 0 0.5`.
    Gate {
        inst: CircuitInstruction,
        qubits: Vec<usize>,
        params: Vec<f64>,
    },
    /// Advance one instruction (also the meaning of an empty line).
    Step,
    /// Run to the next breakpoint or program end.
    Continue,
    /// Reset the loaded program / device to its initial state.
    Reset,
    /// Load a `.sst`/`.ssb` file.
    Load(String),
    /// Toggle the help overlay.
    Help,
    /// Toggle detailed state rendering in the State panel.
    ToggleState,
    /// Leave the TUI.
    Quit,
}

/// Parse one command-line string.
pub fn parse_command(line: &str) -> Result<Command> {
    let line = line.trim();
    if line.is_empty() {
        return Ok(Command::Step);
    }

    // `:`-prefixed meta / debug commands.
    if let Some(rest) = line.strip_prefix(':') {
        let mut it = rest.split_whitespace();
        let cmd = it.next().unwrap_or("");
        return match cmd {
            "q" | "quit" => Ok(Command::Quit),
            "c" | "continue" => Ok(Command::Continue),
            "s" | "step" => Ok(Command::Step),
            "reset" => Ok(Command::Reset),
            "help" | "h" | "?" => Ok(Command::Help),
            "state" => Ok(Command::ToggleState),
            "load" => {
                let path = it.next().ok_or_else(|| eyre!(":load needs a file path"))?;
                if it.next().is_some() {
                    bail!(":load takes a single file path");
                }
                Ok(Command::Load(path.to_string()))
            }
            other => bail!("unknown command :{other}"),
        };
    }

    // Bare tokens: `device N` or a gate op.
    let mut it = line.split_whitespace();
    let head = it.next().unwrap();
    let args: Vec<&str> = it.collect();

    if head == "device" {
        if args.len() != 1 {
            bail!("device takes a single qubit count, got {}", args.len());
        }
        let n = args[0].parse::<usize>().wrap_err("invalid qubit count")?;
        return Ok(Command::Device(n));
    }

    let spec =
        gate_spec(head).ok_or_else(|| eyre!("unknown command {head:?}; try :load or device N"))?;
    let expected = spec.qubits + spec.floats;
    if args.len() != expected {
        bail!(
            "{head} takes {} qubit(s) and {} param(s), got {}",
            spec.qubits,
            spec.floats,
            args.len()
        );
    }
    let (qs, ps) = args.split_at(spec.qubits);
    let qubits = qs
        .iter()
        .map(|t| {
            t.parse::<usize>()
                .wrap_err_with(|| format!("invalid qubit index {t:?}"))
        })
        .collect::<Result<Vec<_>>>()?;
    let params = ps
        .iter()
        .map(|t| {
            t.parse::<f64>()
                .wrap_err_with(|| format!("invalid parameter {t:?}"))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(Command::Gate {
        inst: spec.inst,
        qubits,
        params,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_line_is_step() {
        assert_eq!(parse_command("   ").unwrap(), Command::Step);
    }

    #[test]
    fn device_parses_count() {
        assert_eq!(parse_command("device 3").unwrap(), Command::Device(3));
    }

    #[test]
    fn single_qubit_gate() {
        assert_eq!(
            parse_command("h 0").unwrap(),
            Command::Gate {
                inst: CircuitInstruction::H,
                qubits: vec![0],
                params: vec![],
            }
        );
    }

    #[test]
    fn two_qubit_gate_keeps_operand_order() {
        assert_eq!(
            parse_command("cnot 0 1").unwrap(),
            Command::Gate {
                inst: CircuitInstruction::CNOT,
                qubits: vec![0, 1],
                params: vec![],
            }
        );
    }

    #[test]
    fn rotation_parses_float_param() {
        assert_eq!(
            parse_command("rx 0 0.5").unwrap(),
            Command::Gate {
                inst: CircuitInstruction::RX,
                qubits: vec![0],
                params: vec![0.5],
            }
        );
    }

    #[test]
    fn meta_commands() {
        assert_eq!(parse_command(":q").unwrap(), Command::Quit);
        assert_eq!(parse_command(":continue").unwrap(), Command::Continue);
        assert_eq!(parse_command(":s").unwrap(), Command::Step);
        assert_eq!(parse_command(":reset").unwrap(), Command::Reset);
        assert_eq!(parse_command(":help").unwrap(), Command::Help);
        assert_eq!(parse_command(":h").unwrap(), Command::Help);
        assert_eq!(parse_command(":state").unwrap(), Command::ToggleState);
        assert_eq!(
            parse_command(":load foo.sst").unwrap(),
            Command::Load("foo.sst".to_string())
        );
    }

    #[test]
    fn unknown_gate_errors() {
        assert!(parse_command("bogus 0").is_err());
    }

    #[test]
    fn wrong_arity_errors() {
        assert!(parse_command("x").is_err());
        assert!(parse_command("cnot 0").is_err());
    }

    #[test]
    fn load_rejects_trailing_tokens() {
        assert!(parse_command(":load a.sst").is_ok());
        assert!(parse_command(":load a.sst junk").is_err());
    }

    #[test]
    fn tree_scroll_command_is_not_supported() {
        assert!(parse_command(":tree older").is_err());
    }

    #[test]
    fn device_rejects_trailing_tokens() {
        assert!(parse_command("device 2").is_ok());
        assert!(parse_command("device 2 extra").is_err());
    }
}
