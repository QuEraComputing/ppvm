// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use stim_parser::ast::{GateName, MeasureName, NoiseName, PauliFactor, Target};
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
    #[error("invalid 'MPP' Pauli product at line {line}: {message}")]
    InvalidPauliProduct { line: usize, message: String },
}

pub fn validate(program: &ExtendedProgram) -> Result<(), ExecError> {
    let mut measurements = 0usize;
    validate_slice(&program.instructions, &mut measurements)
}

/// Validate a slice of instructions, threading `measurements` — the number of
/// recorded bits produced so far — so `rec[-k]` controls can be range-checked
/// against the live measurement record (see [`check_record_controls`]).
fn validate_slice(
    instructions: &[ExtendedInstruction],
    measurements: &mut usize,
) -> Result<(), ExecError> {
    for instr in instructions {
        match instr {
            ExtendedInstruction::Raw(RawPassthrough::Gate {
                name,
                targets,
                line,
                ..
            }) => {
                check_gate_supported(*name, *line)?;
                check_record_controls(*name, targets, *measurements, *line)?;
            }
            ExtendedInstruction::Raw(RawPassthrough::Noise { name, line, .. }) => {
                check_noise_supported(*name, *line)?;
            }
            ExtendedInstruction::Raw(RawPassthrough::Measure {
                name,
                args,
                targets,
                line,
                ..
            }) => {
                check_measure_supported(*name, *line)?;
                if let Some(&p) = args.first() {
                    check_probability(p, name.canonical_name(), *line)?;
                }
                *measurements = measurements.saturating_add(targets.len());
            }
            ExtendedInstruction::Mpp {
                args,
                products,
                line,
                ..
            } => {
                if let Some(&p) = args.first() {
                    check_probability(p, "MPP", *line)?;
                }
                for product in products {
                    check_mpp_distinct_qubits(product, *line)?;
                }
                *measurements = measurements.saturating_add(products.len());
            }
            ExtendedInstruction::Repeat { count, body, .. } => {
                // Validate one iteration starting from the current record length.
                // That is the *first* iteration — the shortest record any body
                // instruction ever sees — so a `rec[-k]` that is in range here is
                // in range on every later iteration too. Then advance the count by
                // all `count` iterations' worth of measurements.
                let before = *measurements;
                validate_slice(body, measurements)?;
                let per_iter = measurements.saturating_sub(before);
                let count = usize::try_from(*count).unwrap_or(usize::MAX);
                *measurements = before.saturating_add(per_iter.saturating_mul(count));
            }
            ExtendedInstruction::Raw(RawPassthrough::Annotation { .. })
            | ExtendedInstruction::T { .. }
            | ExtendedInstruction::TDag { .. }
            | ExtendedInstruction::Rotation { .. }
            | ExtendedInstruction::U3 { .. }
            | ExtendedInstruction::Loss { .. }
            | ExtendedInstruction::CorrelatedLoss { .. } => {}
            ExtendedInstruction::MPad {
                prob, bits, line, ..
            } => {
                if let Some(p) = prob {
                    check_probability(*p, "MPAD", *line)?;
                }
                *measurements = measurements.saturating_add(bits.len());
            }
        }
    }
    Ok(())
}

fn check_gate_supported(name: GateName, line: usize) -> Result<(), ExecError> {
    use GateName::*;
    match name {
        Reset | ResetZ | ResetX | ResetY | X | Y | Z | H | HXZ | S | SqrtZ | SDag | SqrtZDag
        | SqrtX | SqrtXDag | SqrtY | SqrtYDag | T | TDag | Identity | CX | ZCX | CNot | CY
        | ZCY | CZ | ZCZ => Ok(()),
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
///
/// `measurements` is the number of recorded bits produced before this gate; a
/// `rec[-k]` with `k > measurements` looks back before the start of the record,
/// which Stim raises as an `IndexError`. Rejecting it here keeps such circuits
/// from silently no-op'ing in the executor.
fn check_record_controls(
    name: GateName,
    targets: &[Target],
    measurements: usize,
    line: usize,
) -> Result<(), ExecError> {
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
    // Targets come in (control, target) pairs; the target may never be a record,
    // and a control may never look back past the start of the record.
    for pair in targets.chunks_exact(2) {
        if matches!(pair[1], Target::Rec(_)) {
            return Err(ExecError::InvalidRecordControl {
                name: name.canonical_name().to_string(),
                line,
                message: "measurement record cannot be a gate target".to_string(),
            });
        }
        if let Target::Rec(k) = pair[0]
            && k > measurements
        {
            return Err(ExecError::InvalidRecordControl {
                name: name.canonical_name().to_string(),
                line,
                message: format!(
                    "rec[-{k}] looks back before the start of the measurement record \
                     ({measurements} measurement(s) recorded so far)"
                ),
            });
        }
    }
    Ok(())
}

/// Reject an `MPP` product that names the same qubit more than once, e.g.
/// `MPP X0*X0` or the anti-Hermitian `MPP Z0*X0`. Stim folds repeated factors
/// via Pauli multiplication — collapsing `X0*X0` to identity (always 0) and
/// rejecting `Z0*X0` (`= iY0`) as anti-Hermitian — but the non-destructive
/// CX-ladder gadget in the executor assumes one factor per qubit, so a repeat
/// would silently misbehave. We reject all repeats with a single clear error
/// rather than emulate the folding.
fn check_mpp_distinct_qubits(product: &[PauliFactor], line: usize) -> Result<(), ExecError> {
    for (i, factor) in product.iter().enumerate() {
        if product[..i].iter().any(|prev| prev.qubit == factor.qubit) {
            return Err(ExecError::InvalidPauliProduct {
                line,
                message: format!(
                    "qubit {} appears more than once; MPP products must act on distinct qubits. \
                     If you have a use case for repeated Paulis in a product, please open an \
                     issue at https://github.com/QuEraComputing/ppvm/issues",
                    factor.qubit
                ),
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
