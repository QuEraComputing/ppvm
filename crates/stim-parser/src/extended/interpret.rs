//! Post-pass interpretation from the vanilla Stim AST to the extended AST.

use crate::ast::{GateName, NoiseName, Program, RawInstruction, Tag, TagParam};
use crate::extended::ast::{Axis, ExtendedInstruction, ExtendedProgram};
use crate::extended::parser::ExtendedParseError;

pub(crate) fn interpret(prog: Program) -> Result<ExtendedProgram, ExtendedParseError> {
    let mut out = Vec::with_capacity(prog.instructions.len());
    interpret_slice(prog.instructions, &mut out)?;
    Ok(ExtendedProgram { instructions: out })
}

fn interpret_slice(
    src: Vec<RawInstruction>,
    out: &mut Vec<ExtendedInstruction>,
) -> Result<(), ExtendedParseError> {
    for raw in src {
        out.push(interpret_one(raw)?);
    }
    Ok(())
}

fn interpret_one(raw: RawInstruction) -> Result<ExtendedInstruction, ExtendedParseError> {
    match raw {
        RawInstruction::Gate {
            name,
            tags,
            args,
            targets,
            line,
        } => interpret_gate(name, tags, args, targets, line),
        RawInstruction::Noise {
            name,
            tags,
            args,
            targets,
            line,
        } => interpret_noise(name, tags, args, targets, line),
        m @ RawInstruction::Measure { .. } => Ok(ExtendedInstruction::Raw(m)),
        a @ RawInstruction::Annotation { .. } => Ok(ExtendedInstruction::Raw(a)),
        RawInstruction::MPad {
            tags,
            prob,
            bits,
            line,
        } => Ok(ExtendedInstruction::MPad {
            tags,
            prob,
            bits: convert_mpad_bits(&bits, line)?,
            line,
        }),
        RawInstruction::Repeat { count, body, line } => {
            let mut inner = Vec::with_capacity(body.len());
            interpret_slice(body, &mut inner)?;
            Ok(ExtendedInstruction::Repeat {
                count,
                body: inner,
                line,
            })
        }
    }
}

fn interpret_gate(
    name: GateName,
    tags: Vec<Tag>,
    args: Vec<f64>,
    targets: Vec<usize>,
    line: usize,
) -> Result<ExtendedInstruction, ExtendedParseError> {
    use GateName::*;

    match name {
        S | SDag => match tags.as_slice() {
            [] => Ok(ExtendedInstruction::Raw(RawInstruction::Gate {
                name,
                tags,
                args,
                targets,
                line,
            })),
            [t] if t.name == "T" => {
                require_no_params(t, name.canonical_name(), line)?;
                Ok(if matches!(name, S) {
                    ExtendedInstruction::T { targets, line }
                } else {
                    ExtendedInstruction::TDag { targets, line }
                })
            }
            [t] => Err(invalid_tag(
                t.name.clone(),
                name.canonical_name(),
                line,
                "expected [T]",
            )),
            _ => Err(invalid_tag(
                tags[0].name.clone(),
                name.canonical_name(),
                line,
                "expected exactly one tag",
            )),
        },
        Identity => match tags.as_slice() {
            [] => Ok(ExtendedInstruction::Raw(RawInstruction::Gate {
                name,
                tags,
                args,
                targets,
                line,
            })),
            [t] => interpret_identity_tag(t, targets, line),
            _ => Err(invalid_tag(
                tags[0].name.clone(),
                "I",
                line,
                "expected exactly one tag",
            )),
        },
        _ => Ok(ExtendedInstruction::Raw(RawInstruction::Gate {
            name,
            tags,
            args,
            targets,
            line,
        })),
    }
}

fn require_no_params(tag: &Tag, instruction: &str, line: usize) -> Result<(), ExtendedParseError> {
    if !tag.params.is_empty() {
        return Err(invalid_tag(
            tag.name.clone(),
            instruction,
            line,
            "tag must have no parameters",
        ));
    }
    Ok(())
}

fn interpret_identity_tag(
    tag: &Tag,
    targets: Vec<usize>,
    line: usize,
) -> Result<ExtendedInstruction, ExtendedParseError> {
    let axis = match tag.name.as_str() {
        "R_X" => Some(Axis::X),
        "R_Y" => Some(Axis::Y),
        "R_Z" => Some(Axis::Z),
        _ => None,
    };

    if let Some(axis) = axis {
        let [theta] = exact_named_params(tag, ["theta"], "I", line)?;
        return Ok(ExtendedInstruction::Rotation {
            axis,
            theta,
            targets,
            line,
        });
    }

    if tag.name == "U3" {
        let [theta, phi, lambda] = exact_named_params(tag, ["theta", "phi", "lambda"], "I", line)?;
        return Ok(ExtendedInstruction::U3 {
            theta,
            phi,
            lambda,
            targets,
            line,
        });
    }

    Err(invalid_tag(
        tag.name.clone(),
        "I",
        line,
        "unrecognized tag (expected R_X / R_Y / R_Z / U3)",
    ))
}

fn interpret_noise(
    name: NoiseName,
    tags: Vec<Tag>,
    args: Vec<f64>,
    targets: Vec<usize>,
    line: usize,
) -> Result<ExtendedInstruction, ExtendedParseError> {
    use NoiseName::*;

    match (name, tags.as_slice()) {
        (IError, [t]) if t.name == "loss" => {
            require_no_params(t, "I_ERROR", line)?;
            if args.len() != 1 {
                return Err(invalid_tag(
                    "loss",
                    "I_ERROR",
                    line,
                    format!("[loss] expects 1 arg, got {}", args.len()),
                ));
            }
            Ok(ExtendedInstruction::Loss {
                p: args[0],
                targets,
                line,
            })
        }
        (IError, [t]) if t.name == "correlated_loss" => {
            require_no_params(t, "I_ERROR", line)?;
            if targets.is_empty() || !targets.len().is_multiple_of(2) {
                return Err(invalid_tag(
                    "correlated_loss",
                    "I_ERROR",
                    line,
                    format!(
                        "[correlated_loss] expects a nonzero even target count, got {}",
                        targets.len()
                    ),
                ));
            }
            let ps = match args.len() {
                1 => [args[0], 0.0, 0.0],
                3 => [args[0], args[1], args[2]],
                n => {
                    return Err(invalid_tag(
                        "correlated_loss",
                        "I_ERROR",
                        line,
                        format!("[correlated_loss] expects 1 or 3 args, got {n}"),
                    ));
                }
            };
            Ok(ExtendedInstruction::CorrelatedLoss {
                ps,
                targets: pair_targets(&targets),
                line,
            })
        }
        (IError, []) => Err(invalid_tag(
            "<missing>",
            "I_ERROR",
            line,
            "I_ERROR requires a [loss] or [correlated_loss] tag",
        )),
        (IError, [t]) => Err(invalid_tag(
            t.name.clone(),
            "I_ERROR",
            line,
            "expected [loss] or [correlated_loss]",
        )),
        (IError, _) => Err(invalid_tag(
            tags[0].name.clone(),
            "I_ERROR",
            line,
            "expected exactly one tag",
        )),
        _ => Ok(ExtendedInstruction::Raw(RawInstruction::Noise {
            name,
            tags,
            args,
            targets,
            line,
        })),
    }
}

fn invalid_tag(
    tag_name: impl Into<String>,
    instruction: &str,
    line: usize,
    message: impl Into<String>,
) -> ExtendedParseError {
    ExtendedParseError::InvalidTag {
        tag: tag_name.into(),
        instruction: instruction.to_string(),
        line,
        message: message.into(),
    }
}

fn convert_mpad_bits(bits: &[usize], line: usize) -> Result<Vec<bool>, ExtendedParseError> {
    let mut out = Vec::with_capacity(bits.len());
    for (index, value) in bits.iter().copied().enumerate() {
        match value {
            0 => out.push(false),
            1 => out.push(true),
            _ => return Err(ExtendedParseError::InvalidMPadBit { line, index, value }),
        }
    }
    Ok(out)
}

fn pair_targets(targets: &[usize]) -> Vec<(usize, usize)> {
    targets
        .chunks_exact(2)
        .map(|pair| (pair[0], pair[1]))
        .collect()
}

fn exact_named_params<const N: usize>(
    tag: &Tag,
    required: [&str; N],
    instruction: &str,
    line: usize,
) -> Result<[f64; N], ExtendedParseError> {
    let mut values = [0.0; N];
    let mut seen = [false; N];

    for param in &tag.params {
        match param {
            TagParam::Positional(_) => {
                return Err(invalid_tag(
                    tag.name.clone(),
                    instruction,
                    line,
                    "tag parameters must be named",
                ));
            }
            TagParam::Named { key, value } => {
                let Some(index) = required.iter().position(|required_key| key == required_key)
                else {
                    return Err(invalid_tag(
                        tag.name.clone(),
                        instruction,
                        line,
                        format!("unexpected named parameter '{key}'"),
                    ));
                };
                if seen[index] {
                    return Err(invalid_tag(
                        tag.name.clone(),
                        instruction,
                        line,
                        format!("duplicate named parameter '{key}'"),
                    ));
                }
                seen[index] = true;
                values[index] = *value;
            }
        }
    }

    for (index, key) in required.iter().enumerate() {
        if !seen[index] {
            return Err(invalid_tag(
                tag.name.clone(),
                instruction,
                line,
                format!("missing required named parameter '{key}'"),
            ));
        }
    }

    Ok(values)
}
