// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use stim_parser::ast::{GateName, MeasureName, NoiseName, Target};
use stim_parser::extended::{ExtendedInstruction, ExtendedProgram, RawPassthrough};

#[derive(Debug, thiserror::Error, PartialEq)]
pub enum ExecError {
    #[error("unsupported instruction '{name}' at line {line}")]
    Unsupported { name: String, line: usize },
    #[error("invalid probability {value} for '{name}' at line {line}; expected value in [0, 1]")]
    InvalidProbability {
        name: String,
        line: usize,
        value: f64,
    },
    #[error("invalid measurement-record control for '{name}' at line {line}: {message}")]
    InvalidRecordControl {
        name: String,
        line: usize,
        message: String,
    },
}

pub fn prepare(program: &ExtendedProgram) -> Result<(), ExecError> {
    validate_slice(&program.instructions)
}

fn validate_slice(instructions: &[ExtendedInstruction]) -> Result<(), ExecError> {
    for instr in instructions {
        match instr {
            ExtendedInstruction::Raw(RawPassthrough::Gate {
                name,
                targets,
                line,
                ..
            }) => {
                check_gate_supported(*name, *line)?;
                check_record_controls(*name, targets, *line)?;
            }
            ExtendedInstruction::Raw(RawPassthrough::Noise { name, line, .. }) => {
                check_noise_supported(*name, *line)?;
            }
            ExtendedInstruction::Raw(RawPassthrough::Measure {
                name, args, line, ..
            }) => {
                check_measure_supported(*name, *line)?;
                if let Some(&p) = args.first() {
                    check_probability(p, name.canonical_name(), *line)?;
                }
            }
            ExtendedInstruction::Repeat { body, .. } => validate_slice(body)?,
            ExtendedInstruction::Raw(RawPassthrough::Annotation { .. })
            | ExtendedInstruction::T { .. }
            | ExtendedInstruction::TDag { .. }
            | ExtendedInstruction::Rotation { .. }
            | ExtendedInstruction::U3 { .. }
            | ExtendedInstruction::Loss { .. }
            | ExtendedInstruction::CorrelatedLoss { .. } => {}
            ExtendedInstruction::MPad { prob, line, .. } => {
                if let Some(p) = prob {
                    check_probability(*p, "MPAD", *line)?;
                }
            }
        }
    }
    Ok(())
}

fn check_gate_supported(name: GateName, line: usize) -> Result<(), ExecError> {
    use GateName::*;
    match name {
        Reset | ResetZ | ResetX | ResetY | X | Y | Z | H | HXZ | S | SqrtZ | SDag | SqrtZDag
        | SqrtX | SqrtXDag | SqrtY | SqrtYDag | Identity | CX | ZCX | CNot | CY | ZCY | CZ
        | ZCZ => Ok(()),
        Swap | ISwap | ISwapDag | SqrtXX | SqrtYY | SqrtZZ | CXSwap | SwapCX | XCX | XCY | XCZ
        | YCX | YCY | YCZ | CXYZ | CZYX | HXY | HYZ => Err(ExecError::Unsupported {
            name: name.canonical_name().to_string(),
            line,
        }),
    }
}

/// A measurement-record control (`rec[-k]`) is only meaningful as the *control*
/// of a controlled Pauli — `CX`/`CNOT`, `CY`, `CZ`. Reject it anywhere else, and
/// reject it in the *target* slot (Stim: "measurement record editing is not
/// supported"). Mirrors `TableauSimulator::single_cx`'s target check.
fn check_record_controls(name: GateName, targets: &[Target], line: usize) -> Result<(), ExecError> {
    use GateName::*;
    if !targets.iter().any(|t| matches!(t, Target::Rec(_))) {
        return Ok(());
    }
    let controlled = matches!(name, CX | ZCX | CNot | CY | ZCY | CZ | ZCZ);
    if !controlled {
        return Err(ExecError::InvalidRecordControl {
            name: name.canonical_name().to_string(),
            line,
            message: "measurement-record controls (rec[-k]) are only valid on CX/CNOT, CY and CZ"
                .to_string(),
        });
    }
    // Targets come in (control, target) pairs; the target may never be a record.
    for pair in targets.chunks_exact(2) {
        if matches!(pair[1], Target::Rec(_)) {
            return Err(ExecError::InvalidRecordControl {
                name: name.canonical_name().to_string(),
                line,
                message: "measurement record cannot be a gate target".to_string(),
            });
        }
    }
    Ok(())
}

fn check_noise_supported(name: NoiseName, line: usize) -> Result<(), ExecError> {
    use NoiseName::*;
    match name {
        Depolarize1 | Depolarize2 | PauliChannel1 | PauliChannel2 | XError | YError | ZError => {
            Ok(())
        }
        IError | HeraldedErase | HeraldedPauliChannel1 | CorrelatedError | ElseCorrelatedError => {
            Err(ExecError::Unsupported {
                name: name.canonical_name().to_string(),
                line,
            })
        }
    }
}

fn check_measure_supported(name: MeasureName, line: usize) -> Result<(), ExecError> {
    use MeasureName::*;
    match name {
        M | MZ | MR | MX | MY | MRX | MRY => Ok(()),
        MXX | MYY | MZZ | MPP => Err(ExecError::Unsupported {
            name: name.canonical_name().to_string(),
            line,
        }),
    }
}

fn check_probability(p: f64, name: &str, line: usize) -> Result<(), ExecError> {
    if p.is_finite() && (0.0..=1.0).contains(&p) {
        Ok(())
    } else {
        Err(ExecError::InvalidProbability {
            name: name.to_string(),
            line,
            value: p,
        })
    }
}
