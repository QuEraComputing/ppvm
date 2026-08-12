// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use stim_parser::prelude::{ExtendedInstruction, ExtendedProgram};

use crate::validate::{ExecError, validate};

use super::{Shot, shots};
use crate::executor::adapter::{SelectedBackend, TableauType};
use crate::executor::{BackendTableau, execute_validated_with_rng};

pub(super) fn run<C, I, S, R>(
    instructions: &[ExtendedInstruction],
    measurement_count: usize,
    shots: Vec<Shot<C, I, S, R>>,
) -> Vec<Vec<Option<bool>>>
where
    SelectedBackend: TableauType<C, I, S>,
    BackendTableau<C, I, S>: Send,
    R: rand::Rng + Send,
{
    use rayon::prelude::*;
    shots
        .into_par_iter()
        .map(|(mut tab, mut rng)| {
            let mut results = Vec::with_capacity(measurement_count);
            execute_validated_with_rng(instructions, &mut tab, &mut results, &mut rng);
            results
        })
        .collect()
}

/// Execute every shot in Rayon's global thread pool.
pub fn sample_parallel_with_rng<C, I, S, F, G, R>(
    program: &ExtendedProgram,
    num_shots: usize,
    make_tableau: F,
    make_rng: G,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
where
    SelectedBackend: TableauType<C, I, S>,
    BackendTableau<C, I, S>: Send,
    F: Fn(usize) -> BackendTableau<C, I, S>,
    G: FnMut(usize) -> R,
    R: rand::Rng + Send,
{
    validate(program)?;
    Ok(sample_parallel_validated_with_rng(
        &program.instructions,
        program.measurement_count(),
        num_shots,
        make_tableau,
        make_rng,
    ))
}

/// Parallel validated sampling with serially prederived streams.
pub fn sample_parallel_validated_with_rng<C, I, S, F, G, R>(
    instructions: &[ExtendedInstruction],
    measurement_count: usize,
    num_shots: usize,
    make_tableau: F,
    make_rng: G,
) -> Vec<Vec<Option<bool>>>
where
    SelectedBackend: TableauType<C, I, S>,
    BackendTableau<C, I, S>: Send,
    F: Fn(usize) -> BackendTableau<C, I, S>,
    G: FnMut(usize) -> R,
    R: rand::Rng + Send,
{
    run(
        instructions,
        measurement_count,
        shots(num_shots, make_tableau, make_rng),
    )
}

/// Execute every shot in Rayon's global thread pool with entropy-seeded streams.
pub fn sample_parallel<C, I, S, F>(
    program: &ExtendedProgram,
    num_shots: usize,
    make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
where
    SelectedBackend: TableauType<C, I, S>,
    BackendTableau<C, I, S>: Send,
    F: Fn(usize) -> BackendTableau<C, I, S>,
{
    sample_parallel_with_rng(program, num_shots, make_tableau, |_| {
        rand::make_rng::<rand::rngs::SmallRng>()
    })
}

/// Parallel sampling for an already validated instruction stream.
pub fn sample_parallel_validated<C, I, S, F>(
    instructions: &[ExtendedInstruction],
    measurement_count: usize,
    num_shots: usize,
    make_tableau: F,
) -> Vec<Vec<Option<bool>>>
where
    SelectedBackend: TableauType<C, I, S>,
    BackendTableau<C, I, S>: Send,
    F: Fn(usize) -> BackendTableau<C, I, S>,
{
    sample_parallel_validated_with_rng(
        instructions,
        measurement_count,
        num_shots,
        make_tableau,
        |_| rand::make_rng::<rand::rngs::SmallRng>(),
    )
}
