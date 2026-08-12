// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use stim_parser::prelude::{ExtendedInstruction, ExtendedProgram};

use crate::validate::{ExecError, validate};

use super::{run_serial, shots};
use crate::executor::BackendTableau;
use crate::executor::adapter::{SelectedBackend, TableauType};

/// Execute shots serially with entropy-seeded per-shot streams.
pub fn sample_serial<C, I, S, F>(
    program: &ExtendedProgram,
    num_shots: usize,
    make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
where
    SelectedBackend: TableauType<C, I, S>,
    F: Fn(usize) -> BackendTableau<C, I, S>,
{
    sample_serial_with_rng(program, num_shots, make_tableau, |_| {
        rand::make_rng::<rand::rngs::SmallRng>()
    })
}

/// Execute shots serially with streams made in ascending shot order.
pub fn sample_serial_with_rng<C, I, S, F, G, R>(
    program: &ExtendedProgram,
    num_shots: usize,
    make_tableau: F,
    make_rng: G,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
where
    SelectedBackend: TableauType<C, I, S>,
    F: Fn(usize) -> BackendTableau<C, I, S>,
    G: FnMut(usize) -> R,
    R: rand::Rng,
{
    validate(program)?;
    Ok(sample_serial_validated_with_rng(
        &program.instructions,
        program.measurement_count(),
        num_shots,
        make_tableau,
        make_rng,
    ))
}

/// Serial sampling for an already validated instruction stream.
pub fn sample_serial_validated<C, I, S, F>(
    instructions: &[ExtendedInstruction],
    measurement_count: usize,
    num_shots: usize,
    make_tableau: F,
) -> Vec<Vec<Option<bool>>>
where
    SelectedBackend: TableauType<C, I, S>,
    F: Fn(usize) -> BackendTableau<C, I, S>,
{
    sample_serial_validated_with_rng(
        instructions,
        measurement_count,
        num_shots,
        make_tableau,
        |_| rand::make_rng::<rand::rngs::SmallRng>(),
    )
}

/// Serial validated sampling with explicitly constructed per-shot streams.
pub fn sample_serial_validated_with_rng<C, I, S, F, G, R>(
    instructions: &[ExtendedInstruction],
    measurement_count: usize,
    num_shots: usize,
    make_tableau: F,
    make_rng: G,
) -> Vec<Vec<Option<bool>>>
where
    SelectedBackend: TableauType<C, I, S>,
    F: Fn(usize) -> BackendTableau<C, I, S>,
    G: FnMut(usize) -> R,
    R: rand::Rng,
{
    run_serial(
        instructions,
        measurement_count,
        shots(num_shots, make_tableau, make_rng),
    )
}
