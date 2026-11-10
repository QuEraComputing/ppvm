// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The TUI command grammar. Bare tokens are gate ops; `:`-prefixed tokens are
//! meta/debug commands; an empty line means "step". Gate names are the ISA
//! mnemonics, i.e. the same spelling as `circuit::circuit.<name>` in `.sst`.

use eyre::{Result, WrapErr, bail, eyre};
use ppvm_vihaco::circuit;
use ppvm_vihaco::component::{CircuitOp, TwoQubitPauliErrorMsg};
use ppvm_vihaco::composite::{PPVM, Record};
use vihaco::PopFromStack;

/// One parsed command-line entry.
#[derive(Debug, PartialEq)]
pub enum Command {
    /// `device N` — (re)create a fresh N-qubit tableau device.
    Device(usize),
    /// A gate op, e.g. `cnot 0 1` or `rx 0 0.5`; run with [`send_gate`].
    Gate { name: String, args: Vec<String> },
    /// Advance one instruction (also the meaning of an empty line).
    Step,
    /// Run to the next cpu::cpu.breakpoint or program end.
    Continue,
    /// Reset the loaded program / device to its initial state.
    Reset,
    /// Load a `.sst`/`.ssb` file.
    Load(String),
    /// Toggle the help overlay.
    Help,
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
    let args: Vec<String> = it.map(str::to_string).collect();

    if head == "device" {
        if args.len() != 1 {
            bail!("device takes a single qubit count, got {}", args.len());
        }
        let n = args[0].parse::<usize>().wrap_err("invalid qubit count")?;
        return Ok(Command::Device(n));
    }

    Ok(Command::Gate {
        name: head.to_string(),
        args,
    })
}

/// Send gate `name` to the circuit, building its message from `args`.
///
/// Operands are parsed and range-checked before anything executes, so a
/// malformed command leaves the device untouched.
pub fn send_gate(m: &mut PPVM, name: &str, args: &[&str]) -> Result<()> {
    use circuit::*;
    match name {
        "x" => send(m, name, args, &X),
        "y" => send(m, name, args, &Y),
        "z" => send(m, name, args, &Z),
        "h" => send(m, name, args, &H),
        "s" => send(m, name, args, &S),
        "s_adj" => send(m, name, args, &SAdj),
        "sqrt_x" => send(m, name, args, &SqrtX),
        "sqrt_y" => send(m, name, args, &SqrtY),
        "sqrt_x_adj" => send(m, name, args, &SqrtXAdj),
        "sqrt_y_adj" => send(m, name, args, &SqrtYAdj),
        "t" => send(m, name, args, &T),
        "t_adj" => send(m, name, args, &TAdj),
        "measure" => send(m, name, args, &Measure),
        "reset" => send(m, name, args, &Reset),
        "cnot" => send(m, name, args, &CNOT),
        "cz" => send(m, name, args, &CZ),
        "rx" => send(m, name, args, &RX),
        "ry" => send(m, name, args, &RY),
        "rz" => send(m, name, args, &RZ),
        "r" => send(m, name, args, &R),
        "rxx" => send(m, name, args, &RXX),
        "ryy" => send(m, name, args, &RYY),
        "rzz" => send(m, name, args, &RZZ),
        "u3" => send(m, name, args, &U3),
        "depolarize" => send(m, name, args, &Depolarize),
        "depolarize2" => send(m, name, args, &Depolarize2),
        "loss" => send(m, name, args, &Loss),
        "pauli_error" => send(m, name, args, &PauliError),
        "two_qubit_pauli_error" => send(m, name, args, &TwoQubitPauliError),
        "correlated_loss" => send(m, name, args, &CorrelatedLoss),
        // Push order is q, p0, p1 (probabilities of leaking into pinned |0⟩ / |1⟩).
        "leakage" => send(m, name, args, &Leakage),
        _ => bail!("unknown command {name:?}; try :load or device N"),
    }
}

fn send<I>(m: &mut PPVM, name: &str, args: &[&str], instruction: &I) -> Result<()>
where
    I: CircuitOp,
    I::Message: FromTokens,
    PPVM: Record<I::Effect>,
{
    let message = I::Message::from_tokens(args, m.n_qubits()).map_err(|e| eyre!("{name}: {e}"))?;
    m.send(instruction, message)
}

/// A circuit message built from command-line operand tokens.
trait FromTokens: Sized {
    fn from_tokens(args: &[&str], n_qubits: usize) -> Result<Self>;
}

impl FromTokens for PopFromStack<(u64,)> {
    fn from_tokens(args: &[&str], n: usize) -> Result<Self> {
        let [q] = exactly(args)?;
        Ok(PopFromStack((qubit(q, n)?,)))
    }
}

impl FromTokens for PopFromStack<(u64, u64)> {
    fn from_tokens(args: &[&str], n: usize) -> Result<Self> {
        let [q0, q1] = exactly(args)?;
        Ok(PopFromStack((qubit(q0, n)?, qubit(q1, n)?)))
    }
}

impl FromTokens for PopFromStack<(u64, f64)> {
    fn from_tokens(args: &[&str], n: usize) -> Result<Self> {
        let [q, p] = exactly(args)?;
        Ok(PopFromStack((qubit(q, n)?, param(p)?)))
    }
}

impl FromTokens for PopFromStack<(u64, u64, f64)> {
    fn from_tokens(args: &[&str], n: usize) -> Result<Self> {
        let [q0, q1, p] = exactly(args)?;
        Ok(PopFromStack((qubit(q0, n)?, qubit(q1, n)?, param(p)?)))
    }
}

impl FromTokens for PopFromStack<(u64, f64, f64)> {
    fn from_tokens(args: &[&str], n: usize) -> Result<Self> {
        let [q, p0, p1] = exactly(args)?;
        Ok(PopFromStack((qubit(q, n)?, param(p0)?, param(p1)?)))
    }
}

impl FromTokens for PopFromStack<(u64, f64, f64, f64)> {
    fn from_tokens(args: &[&str], n: usize) -> Result<Self> {
        let [q, p0, p1, p2] = exactly(args)?;
        Ok(PopFromStack((
            qubit(q, n)?,
            param(p0)?,
            param(p1)?,
            param(p2)?,
        )))
    }
}

impl FromTokens for PopFromStack<(u64, u64, f64, f64, f64)> {
    fn from_tokens(args: &[&str], n: usize) -> Result<Self> {
        let [q0, q1, p0, p1, p2] = exactly(args)?;
        Ok(PopFromStack((
            qubit(q0, n)?,
            qubit(q1, n)?,
            param(p0)?,
            param(p1)?,
            param(p2)?,
        )))
    }
}

impl FromTokens for TwoQubitPauliErrorMsg {
    fn from_tokens(args: &[&str], n: usize) -> Result<Self> {
        let [q0, q1, ps @ ..] = exactly::<17>(args)?;
        let mut probs = [0.0; 15];
        for (p, tok) in probs.iter_mut().zip(ps) {
            *p = param(tok)?;
        }
        Ok(TwoQubitPauliErrorMsg {
            q0: qubit(q0, n)? as usize,
            q1: qubit(q1, n)? as usize,
            ps: probs,
        })
    }
}

fn exactly<'a, const N: usize>(args: &[&'a str]) -> Result<[&'a str; N]> {
    <[&str; N]>::try_from(args).map_err(|_| eyre!("expected {N} operands, got {}", args.len()))
}

fn qubit(tok: &str, n_qubits: usize) -> Result<u64> {
    let q = tok
        .parse::<usize>()
        .wrap_err_with(|| format!("invalid qubit index {tok:?}"))?;
    if q >= n_qubits {
        bail!("qubit {q} out of range for {n_qubits}-qubit device");
    }
    Ok(q as u64)
}

fn param(tok: &str) -> Result<f64> {
    tok.parse::<f64>()
        .wrap_err_with(|| format!("invalid parameter {tok:?}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_vihaco::composite::PPVMDeviceInfo;
    use ppvm_vihaco::measurements::MeasurementOutcome::{self, One, Zero};

    fn device(n_qubits: usize) -> PPVM {
        let info = PPVMDeviceInfo {
            n_qubits,
            ..Default::default()
        };
        PPVM::new(&info, Some(0)).unwrap()
    }

    /// Run gate lines in order, then return every recorded outcome.
    fn run(n_qubits: usize, lines: &[&str]) -> Result<Vec<MeasurementOutcome>> {
        let mut m = device(n_qubits);
        for line in lines {
            let mut toks = line.split_whitespace();
            let name = toks.next().unwrap();
            send_gate(&mut m, name, &toks.collect::<Vec<_>>())?;
        }
        Ok(m.measurement_record().iter().flatten().copied().collect())
    }

    #[test]
    fn empty_line_is_step() {
        assert_eq!(parse_command("   ").unwrap(), Command::Step);
    }

    #[test]
    fn device_parses_count() {
        assert_eq!(parse_command("device 3").unwrap(), Command::Device(3));
    }

    #[test]
    fn gate_keeps_operand_order() {
        assert_eq!(
            parse_command("rx 0 0.5").unwrap(),
            Command::Gate {
                name: "rx".to_string(),
                args: vec!["0".to_string(), "0.5".to_string()],
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
        assert_eq!(
            parse_command(":load foo.sst").unwrap(),
            Command::Load("foo.sst".to_string())
        );
    }

    #[test]
    fn load_rejects_trailing_tokens() {
        assert!(parse_command(":load a.sst").is_ok());
        assert!(parse_command(":load a.sst junk").is_err());
    }

    #[test]
    fn device_rejects_trailing_tokens() {
        assert!(parse_command("device 2").is_ok());
        assert!(parse_command("device 2 extra").is_err());
    }

    #[test]
    fn x_flips_and_state_persists_across_commands() -> Result<()> {
        assert_eq!(run(1, &["measure 0"])?, [Zero]);
        assert_eq!(run(1, &["x 0", "measure 0"])?, [One]);
        assert_eq!(run(1, &["x 0", "x 0", "measure 0"])?, [Zero]);
        Ok(())
    }

    #[test]
    fn cnot_control_is_first_operand() -> Result<()> {
        assert_eq!(
            run(2, &["x 0", "cnot 0 1", "measure 0", "measure 1"])?,
            [One, One]
        );
        assert_eq!(
            run(2, &["x 1", "cnot 0 1", "measure 0", "measure 1"])?,
            [Zero, One]
        );
        Ok(())
    }

    #[test]
    fn adjoint_spellings_invert_their_gate() -> Result<()> {
        // In the X basis, S·S† = 1 and T·T† = 1; the H sandwich makes a missing
        // or wrong inverse show up as a |1⟩ outcome.
        for (g, g_adj) in [("s", "s_adj"), ("t", "t_adj")] {
            let lines = [
                "h 0",
                &format!("{g} 0"),
                &format!("{g_adj} 0"),
                "h 0",
                "measure 0",
            ];
            assert_eq!(run(1, &lines)?, [Zero], "{g} / {g_adj}");
        }
        for (g, g_adj) in [("sqrt_x", "sqrt_x_adj"), ("sqrt_y", "sqrt_y_adj")] {
            let lines = [
                &format!("{g} 0") as &str,
                &format!("{g_adj} 0"),
                "measure 0",
            ];
            assert_eq!(run(1, &lines)?, [Zero], "{g} / {g_adj}");
        }
        // √X · √X = X.
        assert_eq!(run(1, &["sqrt_x 0", "sqrt_x 0", "measure 0"])?, [One]);
        Ok(())
    }

    #[test]
    fn rx_pi_flips() -> Result<()> {
        let pi = std::f64::consts::PI.to_string();
        assert_eq!(run(1, &[&format!("rx 0 {pi}"), "measure 0"])?, [One]);
        Ok(())
    }

    #[test]
    fn bad_operands_error_and_leave_state_untouched() {
        let mut m = device(2);
        for (name, args) in [
            ("x", &[][..]),
            ("cnot", &["0"][..]),
            ("rx", &["0"][..]),
            ("x", &["2"][..]),
            ("x", &["q"][..]),
            ("rx", &["0", "half"][..]),
            ("bogus", &["0"][..]),
        ] {
            assert!(send_gate(&mut m, name, args).is_err(), "{name} {args:?}");
        }
        send_gate(&mut m, "measure", &["0"]).unwrap();
        send_gate(&mut m, "measure", &["1"]).unwrap();
        let outcomes: Vec<_> = m.measurement_record().iter().flatten().copied().collect();
        assert_eq!(outcomes, [Zero, Zero]);
    }
}
