use stim_parser::ast::{GateName, MeasureName, NoiseName, RawInstruction};
use stim_parser::extended::{ExtendedInstruction, ExtendedProgram};

#[derive(Debug, thiserror::Error, PartialEq)]
pub enum ExecError {
    #[error("unsupported instruction '{name}' at line {line}")]
    Unsupported { name: String, line: usize },
    /// Raised for `ExtendedInstruction::Raw(_)` values that the extended
    /// interpreter would have lowered to typed variants (`MPad`, `Repeat`).
    /// `parse_extended` never produces these; only a caller hand-constructing
    /// an `ExtendedProgram` can. Reported as a recoverable error rather than
    /// a panic.
    #[error(
        "malformed ExtendedProgram: Raw({kind}) at line {line} should have been lowered to ExtendedInstruction::{kind} by the interpreter"
    )]
    Malformed { kind: &'static str, line: usize },
}

pub fn prepare(program: &ExtendedProgram) -> Result<(), ExecError> {
    validate_slice(&program.instructions)
}

fn validate_slice(instructions: &[ExtendedInstruction]) -> Result<(), ExecError> {
    for instr in instructions {
        match instr {
            ExtendedInstruction::Raw(RawInstruction::Gate { name, line, .. }) => {
                check_gate_supported(*name, *line)?;
            }
            ExtendedInstruction::Raw(RawInstruction::Noise { name, line, .. }) => {
                check_noise_supported(*name, *line)?;
            }
            ExtendedInstruction::Raw(RawInstruction::Measure { name, line, .. }) => {
                check_measure_supported(*name, *line)?;
            }
            ExtendedInstruction::Repeat { body, .. } => validate_slice(body)?,
            ExtendedInstruction::Raw(RawInstruction::Annotation { .. })
            | ExtendedInstruction::MPad { .. }
            | ExtendedInstruction::T { .. }
            | ExtendedInstruction::TDag { .. }
            | ExtendedInstruction::Rotation { .. }
            | ExtendedInstruction::U3 { .. }
            | ExtendedInstruction::Loss { .. }
            | ExtendedInstruction::CorrelatedLoss { .. } => {}
            ExtendedInstruction::Raw(RawInstruction::MPad { line, .. }) => {
                return Err(ExecError::Malformed {
                    kind: "MPad",
                    line: *line,
                });
            }
            ExtendedInstruction::Raw(RawInstruction::Repeat { line, .. }) => {
                return Err(ExecError::Malformed {
                    kind: "Repeat",
                    line: *line,
                });
            }
        }
    }
    Ok(())
}

fn check_gate_supported(name: GateName, line: usize) -> Result<(), ExecError> {
    use GateName::*;
    match name {
        Reset | ResetZ | X | Y | Z | H | HXZ | S | SqrtZ | SDag | SqrtZDag | SqrtX | SqrtXDag
        | SqrtY | SqrtYDag | Identity | CX | ZCX | CNot | CY | ZCY | CZ | ZCZ => Ok(()),
        Swap | ISwap | ISwapDag | SqrtXX | SqrtYY | SqrtZZ | CXSwap | SwapCX | XCX | XCY | XCZ
        | YCX | YCY | YCZ | CXYZ | CZYX | HXY | HYZ => Err(ExecError::Unsupported {
            name: name.canonical_name().to_string(),
            line,
        }),
    }
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
        M | MZ | MR => Ok(()),
        MX | MY | MRX | MRY | MXX | MYY | MZZ | MPP => Err(ExecError::Unsupported {
            name: name.canonical_name().to_string(),
            line,
        }),
    }
}
