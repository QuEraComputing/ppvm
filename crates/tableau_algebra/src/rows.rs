// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_tableau_2::Tableau;
use ppvm_traits_2::{Pauli, StabilizerFrame};

// Row indices here refer to the stabilizer half.
pub(crate) fn site(tableau: &Tableau, row: usize, column: usize) -> Pauli {
    tableau.row_site(tableau.n_qubits() + row, column)
}

pub(crate) fn has_x(site: Pauli) -> bool {
    matches!(site, Pauli::X | Pauli::Y)
}

pub(crate) fn has_z(site: Pauli) -> bool {
    matches!(site, Pauli::Z | Pauli::Y)
}

pub(crate) fn multiply(tableau: &mut Tableau, src: usize, dst: usize) {
    let n = tableau.n_qubits();
    tableau.row_multiply(n + src, n + dst);
    // The inverse basis change on the dual rows preserves their pairing.
    tableau.row_multiply(dst, src);
}

pub(crate) fn swap(tableau: &mut Tableau, a: usize, b: usize) {
    if a != b {
        // Commuting Hermitian rows square to identity: three multiplies swap.
        multiply(tableau, a, b);
        multiply(tableau, b, a);
        multiply(tableau, a, b);
    }
}

pub(crate) fn eliminate(
    tableau: &mut Tableau,
    pivot: usize,
    column: usize,
    component: fn(Pauli) -> bool,
) {
    for row in 0..tableau.n_qubits() {
        if row != pivot && component(site(tableau, row, column)) {
            multiply(tableau, pivot, row);
        }
    }
}

// Reduce one component, retaining the pivot columns for Gottesman form.
pub(crate) fn forward_pass(
    tableau: &mut Tableau,
    start: usize,
    columns: std::ops::Range<usize>,
    component: fn(Pauli) -> bool,
) -> Vec<usize> {
    let mut pivots = Vec::new();
    for column in columns {
        let pivot = start + pivots.len();
        if let Some(row) =
            (pivot..tableau.n_qubits()).find(|&row| component(site(tableau, row, column)))
        {
            swap(tableau, row, pivot);
            eliminate(tableau, pivot, column, component);
            pivots.push(column);
        }
    }
    pivots
}
