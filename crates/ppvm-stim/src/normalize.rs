use crate::parser::ast::{
    AnnotationKind, GateName, MeasureName, NoiseName, Program, RawInstruction,
    Tag, TagParam,
};
use crate::tableau_program::{
    GateKind, Instruction, MeasureKind, NoiseKind, TableauProgram,
};

#[derive(Debug, thiserror::Error, PartialEq)]
#[non_exhaustive]
pub enum NormalizeError {
    #[error("unsupported instruction '{name}' at line {line} (phase 1)")]
    Unsupported { name: String, line: usize },

    #[error("invalid tag '{tag}' on '{instruction}' at line {line}: {message}")]
    InvalidTag {
        tag: String,
        instruction: String,
        line: usize,
        message: String,
    },
}

pub fn to_tableau(program: &Program) -> Result<TableauProgram, NormalizeError> {
    let mut out = Vec::with_capacity(program.instructions.len());
    let mut count = 0usize;
    normalize_slice(&program.instructions, &mut out, &mut count, 1)?;
    Ok(TableauProgram {
        instructions: out,
        expected_measurement_count: count,
    })
}

fn normalize_slice(
    src: &[RawInstruction],
    out: &mut Vec<Instruction>,
    measure_count: &mut usize,
    enclosing_repeat_factor: u64,
) -> Result<(), NormalizeError> {
    for raw in src {
        match raw {
            RawInstruction::Gate { name, tags, args, targets, line } => {
                let kind = gate_to_kind(*name, tags, args, *line)?;
                out.push(Instruction::Gate { kind, targets: targets.clone(), line: *line });
            }
            RawInstruction::Noise { name, tags, args, targets, line } => {
                let (kind, normalized_args) = noise_to_kind(*name, tags, args, *line)?;
                out.push(Instruction::Noise {
                    kind,
                    targets: targets.clone(),
                    args: normalized_args,
                    line: *line,
                });
            }
            RawInstruction::Measure { name, args, targets, line, .. } => {
                let kind = measure_to_kind(*name, *line)?;
                let noise = args.first().copied().unwrap_or(0.0);
                *measure_count = measure_count.saturating_add(
                    targets.len().saturating_mul(enclosing_repeat_factor as usize),
                );
                out.push(Instruction::Measure {
                    kind,
                    targets: targets.clone(),
                    noise,
                    line: *line,
                });
            }
            RawInstruction::Annotation { line, .. } => {
                out.push(Instruction::Annotation { line: *line });
            }
            RawInstruction::Repeat { count, body, line } => {
                let mut inner = Vec::with_capacity(body.len());
                normalize_slice(
                    body,
                    &mut inner,
                    measure_count,
                    enclosing_repeat_factor.saturating_mul(*count),
                )?;
                out.push(Instruction::Repeat { count: *count, body: inner, line: *line });
            }
        }
    }
    Ok(())
}

fn gate_to_kind(
    name: GateName,
    tags: &[Tag],
    _args: &[f64],
    line: usize,
) -> Result<GateKind, NormalizeError> {
    use GateName::*;
    Ok(match name {
        Reset | ResetZ => GateKind::Reset,
        X => GateKind::X,
        Y => GateKind::Y,
        Z => GateKind::Z,
        H | HXZ => GateKind::H,
        SqrtZ => GateKind::S,
        SqrtZDag => GateKind::SDag,
        SqrtX => GateKind::SqrtX,
        SqrtXDag => GateKind::SqrtXDag,
        SqrtY => GateKind::SqrtY,
        SqrtYDag => GateKind::SqrtYDag,
        S => match find_tag(tags, "T") {
            Some(t) => {
                require_no_params(t, "S", line)?;
                GateKind::T
            }
            None => GateKind::S,
        },
        SDag => match find_tag(tags, "T") {
            Some(t) => {
                require_no_params(t, "S_DAG", line)?;
                GateKind::TDag
            }
            None => GateKind::SDag,
        },
        Identity => identity_to_kind(tags, line)?,
        CX | ZCX | CNot => GateKind::CX,
        CY | ZCY => GateKind::CY,
        CZ | ZCZ => GateKind::CZ,
        Swap | ISwap | ISwapDag | SqrtXX | SqrtYY | SqrtZZ
        | CXSwap | SwapCX | XCX | XCY | XCZ | YCX | YCY | YCZ
        | CXYZ | CZYX | HXY | HYZ => {
            return Err(NormalizeError::Unsupported {
                name: name.canonical_name().to_string(),
                line,
            });
        }
    })
}

fn identity_to_kind(tags: &[Tag], line: usize) -> Result<GateKind, NormalizeError> {
    let tag = match tags {
        [t] => t,
        [] => return Err(NormalizeError::InvalidTag {
            tag: String::new(),
            instruction: "I".to_string(),
            line,
            message: "expected a single dialect tag like [R_X(theta=…)]".into(),
        }),
        _ => return Err(NormalizeError::InvalidTag {
            tag: tags[0].name.clone(),
            instruction: "I".to_string(),
            line,
            message: "expected exactly one tag".into(),
        }),
    };

    let lookup_named = |key: &str| -> Result<f64, NormalizeError> {
        tag.params
            .iter()
            .find_map(|p| match p {
                TagParam::Named { key: k, value } if k == key => Some(*value),
                _ => None,
            })
            .ok_or(NormalizeError::InvalidTag {
                tag: tag.name.clone(),
                instruction: "I".to_string(),
                line,
                message: format!("missing required named parameter '{key}'"),
            })
    };

    Ok(match tag.name.as_str() {
        "R_X" => GateKind::RX { theta: lookup_named("theta")? },
        "R_Y" => GateKind::RY { theta: lookup_named("theta")? },
        "R_Z" => GateKind::RZ { theta: lookup_named("theta")? },
        "U3" => GateKind::U3 {
            theta: lookup_named("theta")?,
            phi: lookup_named("phi")?,
            lambda: lookup_named("lambda")?,
        },
        other => return Err(NormalizeError::InvalidTag {
            tag: other.to_string(),
            instruction: "I".to_string(),
            line,
            message: "unrecognized identity-tag name".into(),
        }),
    })
}

fn noise_to_kind(
    name: NoiseName,
    tags: &[Tag],
    args: &[f64],
    line: usize,
) -> Result<(NoiseKind, Vec<f64>), NormalizeError> {
    use NoiseName::*;
    Ok(match name {
        Depolarize1 => (NoiseKind::Depolarize1, args.to_vec()),
        Depolarize2 => (NoiseKind::Depolarize2, args.to_vec()),
        PauliChannel1 => (NoiseKind::PauliChannel1, args.to_vec()),
        PauliChannel2 => (NoiseKind::PauliChannel2, args.to_vec()),
        XError => (NoiseKind::XError, args.to_vec()),
        YError => (NoiseKind::YError, args.to_vec()),
        ZError => (NoiseKind::ZError, args.to_vec()),
        IError => match tags {
            [t] if t.name == "loss" => {
                if args.len() != 1 {
                    return Err(NormalizeError::InvalidTag {
                        tag: "loss".into(),
                        instruction: "I_ERROR".into(),
                        line,
                        message: format!("[loss] expects 1 arg, got {}", args.len()),
                    });
                }
                (NoiseKind::Loss, args.to_vec())
            }
            [t] if t.name == "correlated_loss" => {
                let normalized = match args.len() {
                    1 => vec![args[0], 0.0, 0.0],
                    3 => args.to_vec(),
                    n => return Err(NormalizeError::InvalidTag {
                        tag: "correlated_loss".into(),
                        instruction: "I_ERROR".into(),
                        line,
                        message: format!("[correlated_loss] expects 1 or 3 args, got {n}"),
                    }),
                };
                (NoiseKind::CorrelatedLoss, normalized)
            }
            [] => return Err(NormalizeError::InvalidTag {
                tag: String::new(),
                instruction: "I_ERROR".into(),
                line,
                message: "I_ERROR requires a [loss] or [correlated_loss] tag".into(),
            }),
            [t] => return Err(NormalizeError::InvalidTag {
                tag: t.name.clone(),
                instruction: "I_ERROR".into(),
                line,
                message: "expected [loss] or [correlated_loss]".into(),
            }),
            _ => return Err(NormalizeError::InvalidTag {
                tag: tags[0].name.clone(),
                instruction: "I_ERROR".into(),
                line,
                message: "expected exactly one tag".into(),
            }),
        },
        HeraldedErase | HeraldedPauliChannel1 | CorrelatedError | ElseCorrelatedError => {
            return Err(NormalizeError::Unsupported {
                name: name.canonical_name().to_string(),
                line,
            });
        }
    })
}

fn measure_to_kind(name: MeasureName, line: usize) -> Result<MeasureKind, NormalizeError> {
    use MeasureName::*;
    Ok(match name {
        M | MZ => MeasureKind::M,
        MR => MeasureKind::MR,
        other => return Err(NormalizeError::Unsupported {
            name: other.canonical_name().to_string(),
            line,
        }),
    })
}

fn find_tag<'a>(tags: &'a [Tag], name: &str) -> Option<&'a Tag> {
    tags.iter().find(|t| t.name == name)
}

fn require_no_params(tag: &Tag, instruction: &str, line: usize) -> Result<(), NormalizeError> {
    if !tag.params.is_empty() {
        return Err(NormalizeError::InvalidTag {
            tag: tag.name.clone(),
            instruction: instruction.to_string(),
            line,
            message: "tag must have no parameters".into(),
        });
    }
    Ok(())
}

#[allow(dead_code)]
const _: AnnotationKind = AnnotationKind::Detector;
