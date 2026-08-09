// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

#[cfg(feature = "rayon")]
mod parallel;
mod serial;

#[cfg(feature = "rayon")]
pub use parallel::{
    sample_parallel, sample_parallel_validated, sample_parallel_validated_with_rng,
    sample_parallel_with_rng,
};
pub use serial::{
    sample_serial, sample_serial_validated, sample_serial_validated_with_rng,
    sample_serial_with_rng,
};

use stim_parser::prelude::{ExtendedInstruction, ExtendedProgram};

use crate::validate::{ExecError, validate};

use super::adapter::{SelectedBackend, TableauType};
use super::{BackendTableau, execute_validated_with_rng};

type Shot<C, I, S, R> = (BackendTableau<C, I, S>, R);

fn shots<C, I, S, F, G, R>(
    num_shots: usize,
    make_tableau: F,
    mut make_rng: G,
) -> Vec<Shot<C, I, S, R>>
where
    SelectedBackend: TableauType<C, I, S>,
    F: Fn(usize) -> BackendTableau<C, I, S>,
    G: FnMut(usize) -> R,
{
    (0..num_shots)
        .map(|i| (make_tableau(i), make_rng(i)))
        .collect()
}

fn run_serial<C, I, S, R>(
    instructions: &[ExtendedInstruction],
    measurement_count: usize,
    shots: Vec<Shot<C, I, S, R>>,
) -> Vec<Vec<Option<bool>>>
where
    SelectedBackend: TableauType<C, I, S>,
    R: rand::Rng,
{
    shots
        .into_iter()
        .map(|(mut tab, mut rng)| {
            let mut results = Vec::with_capacity(measurement_count);
            execute_validated_with_rng(instructions, &mut tab, &mut results, &mut rng);
            results
        })
        .collect()
}

/// Execute shots, selecting serial or Rayon execution by batch size.
pub fn sample<C, I, S, F>(
    program: &ExtendedProgram,
    num_shots: usize,
    make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
where
    SelectedBackend: TableauType<C, I, S>,
    BackendTableau<C, I, S>: Send,
    F: Fn(usize) -> BackendTableau<C, I, S>,
{
    sample_with_rng(program, num_shots, make_tableau, |_| {
        rand::make_rng::<rand::rngs::SmallRng>()
    })
}

/// Execute shots with per-shot RNGs constructed serially before any Rayon work.
pub fn sample_with_rng<C, I, S, F, G, R>(
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
    Ok(sample_validated_with_rng(
        &program.instructions,
        program.measurement_count(),
        num_shots,
        make_tableau,
        make_rng,
    ))
}

/// Sampling for an already validated instruction stream.
pub fn sample_validated<C, I, S, F>(
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
    sample_validated_with_rng(
        instructions,
        measurement_count,
        num_shots,
        make_tableau,
        |_| rand::make_rng::<rand::rngs::SmallRng>(),
    )
}

/// Validated sampling with streams prederived in ascending shot order.
pub fn sample_validated_with_rng<C, I, S, F, G, R>(
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
    let shots = shots(num_shots, make_tableau, make_rng);
    #[cfg(feature = "rayon")]
    if rayon::current_num_threads() > 1 && num_shots >= 4 * rayon::current_num_threads() {
        return parallel::run(instructions, measurement_count, shots);
    }
    run_serial(instructions, measurement_count, shots)
}
