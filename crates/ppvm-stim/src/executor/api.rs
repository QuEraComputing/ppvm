// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use stim_parser::prelude::ExtendedProgram;

use crate::validate::{ExecError, validate};

use super::adapter::{SelectedBackend, TableauType};
use super::{execute_validated, execute_validated_with_rng};

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

/// Validate and execute a parsed program with caller-supplied randomness.
pub fn execute_with_rng<C, I, S, R: rand::Rng + ?Sized>(
    program: &ExtendedProgram,
    tab: &mut BackendTableau<C, I, S>,
    rng: &mut R,
) -> Result<Vec<Option<bool>>, ExecError>
where
    SelectedBackend: TableauType<C, I, S>,
{
    validate(program)?;
    let mut results = Vec::with_capacity(program.measurement_count());
    execute_validated_with_rng(&program.instructions, tab, &mut results, rng);
    Ok(results)
}
