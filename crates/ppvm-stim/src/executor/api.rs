// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use stim_parser::prelude::{ExtendedInstruction, ExtendedProgram};

use crate::validate::{ExecError, validate};

use super::adapter::{SelectedBackend, TableauType};
use super::execute_validated;

#[doc(hidden)]
pub type BackendTableau<C, I, S> = <SelectedBackend as TableauType<C, I, S>>::Tableau;

/// Validate and execute a parsed program against a tableau.
pub fn execute<C, I, S>(
    program: &ExtendedProgram,
    tab: &mut BackendTableau<C, I, S>,
) -> Result<Vec<Option<bool>>, ExecError>
where
    SelectedBackend: TableauType<C, I, S>,
{
    validate(program)?;
    let mut results = Vec::with_capacity(program.measurement_count());
    execute_validated(&program.instructions, tab, &mut results);
    Ok(results)
}

/// Execute shots serially, constructing a fresh tableau for each shot index.
pub fn sample_serial<C, I, S, F>(
    program: &ExtendedProgram,
    num_shots: usize,
    make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
where
    SelectedBackend: TableauType<C, I, S>,
    F: Fn(usize) -> BackendTableau<C, I, S>,
{
    validate(program)?;
    Ok(sample_serial_validated(
        &program.instructions,
        program.measurement_count(),
        num_shots,
        make_tableau,
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
    (0..num_shots)
        .map(|i| {
            let mut tab = make_tableau(i);
            let mut results = Vec::with_capacity(measurement_count);
            execute_validated(instructions, &mut tab, &mut results);
            results
        })
        .collect()
}

/// Execute shots, using rayon for sufficiently large batches when enabled.
#[cfg(not(feature = "rayon"))]
pub fn sample<C, I, S, F>(
    program: &ExtendedProgram,
    num_shots: usize,
    make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
where
    SelectedBackend: TableauType<C, I, S>,
    F: Fn(usize) -> BackendTableau<C, I, S>,
{
    validate(program)?;
    Ok(sample_serial_validated(
        &program.instructions,
        program.measurement_count(),
        num_shots,
        make_tableau,
    ))
}

/// Execute shots, using rayon for sufficiently large batches when enabled.
#[cfg(feature = "rayon")]
pub fn sample<C, I, S, F>(
    program: &ExtendedProgram,
    num_shots: usize,
    make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
where
    SelectedBackend: TableauType<C, I, S>,
    BackendTableau<C, I, S>: Send,
    F: Fn(usize) -> BackendTableau<C, I, S> + Sync,
{
    validate(program)?;
    Ok(sample_validated(
        &program.instructions,
        program.measurement_count(),
        num_shots,
        make_tableau,
    ))
}

/// Sampling for an already validated instruction stream.
#[cfg(not(feature = "rayon"))]
pub fn sample_validated<C, I, S, F>(
    instructions: &[ExtendedInstruction],
    measurement_count: usize,
    num_shots: usize,
    make_tableau: F,
) -> Vec<Vec<Option<bool>>>
where
    SelectedBackend: TableauType<C, I, S>,
    F: Fn(usize) -> BackendTableau<C, I, S>,
{
    sample_serial_validated(instructions, measurement_count, num_shots, make_tableau)
}

/// Sampling for an already validated instruction stream.
#[cfg(feature = "rayon")]
pub fn sample_validated<C, I, S, F>(
    instructions: &[ExtendedInstruction],
    measurement_count: usize,
    num_shots: usize,
    make_tableau: F,
) -> Vec<Vec<Option<bool>>>
where
    SelectedBackend: TableauType<C, I, S>,
    BackendTableau<C, I, S>: Send,
    F: Fn(usize) -> BackendTableau<C, I, S> + Sync,
{
    let n_threads = rayon::current_num_threads();
    if n_threads > 1 && num_shots >= 4 * n_threads {
        return sample_parallel_validated(instructions, measurement_count, num_shots, make_tableau);
    }
    sample_serial_validated(instructions, measurement_count, num_shots, make_tableau)
}

/// Execute shots in parallel across rayon's global thread pool.
#[cfg(feature = "rayon")]
pub fn sample_parallel<C, I, S, F>(
    program: &ExtendedProgram,
    num_shots: usize,
    make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
where
    SelectedBackend: TableauType<C, I, S>,
    BackendTableau<C, I, S>: Send,
    F: Fn(usize) -> BackendTableau<C, I, S> + Sync,
{
    validate(program)?;
    Ok(sample_parallel_validated(
        &program.instructions,
        program.measurement_count(),
        num_shots,
        make_tableau,
    ))
}

/// Parallel sampling for an already validated instruction stream.
#[cfg(feature = "rayon")]
pub fn sample_parallel_validated<C, I, S, F>(
    instructions: &[ExtendedInstruction],
    measurement_count: usize,
    num_shots: usize,
    make_tableau: F,
) -> Vec<Vec<Option<bool>>>
where
    SelectedBackend: TableauType<C, I, S>,
    BackendTableau<C, I, S>: Send,
    F: Fn(usize) -> BackendTableau<C, I, S> + Sync,
{
    use rayon::prelude::*;
    (0..num_shots)
        .into_par_iter()
        .map(|i| {
            let mut tab = make_tableau(i);
            let mut results = Vec::with_capacity(measurement_count);
            execute_validated(instructions, &mut tab, &mut results);
            results
        })
        .collect()
}
