// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_tableau_2::Tableau;
use ppvm_traits_2::Pauli;

use crate::rows::{multiply, site, swap};

/// Put a pure stabilizer state in the clipped gauge, in place.
///
/// > **Warning:** Requires a valid full-rank frame: `n` destabilizers followed by
/// > their `n` paired stabilizers on `n` qubits (`2n × n` Pauli entries).
/// > Stabilizer-only and rank-deficient tableaux are unsupported.
///
/// This is QuantumClifford's `canonicalize_clip!` form (Nahum et al., 2017;
/// Li et al., 2019, Appendix A). A left-to-right pass fixes the left endpoints;
/// a right-to-left pass shortens the rows while preserving those endpoints.
/// Each qubit has two endpoints in total, counting both ends of a single-site
/// row. Rows sharing a left endpoint, or sharing a right endpoint, have
/// different nonidentity Paulis there.
///
/// Endpoint pairs describe the state's entanglement across cuts in the qubit
/// ordering.
///
/// ```
/// use ppvm_tableau_2::Tableau;
/// use ppvm_traits_2::Clifford;
/// use tableau_algebra::canonicalize_clip;
///
/// let mut ghz = Tableau::new(3);
/// ghz.h(0);
/// ghz.cnot(0, 1);
/// ghz.cnot(1, 2);
/// canonicalize_clip(&mut ghz);
/// // +XXX, +ZZI, +IZZ have endpoint pairs (0,2), (0,1), (1,2).
/// ```
pub fn canonicalize_clip(tableau: &mut Tableau) {
    let n = tableau.n_qubits();
    let mut next = 0;
    for column in 0..n {
        let Some(first) = (next..n).find(|&row| site(tableau, row, column) != Pauli::I) else {
            continue;
        };
        let second = (first + 1..n).find(|&row| {
            let p = site(tableau, row, column);
            p != Pauli::I && p != site(tableau, first, column)
        });
        swap(tableau, first, next);
        let second = second.map(|row| {
            swap(tableau, row, next + 1);
            next + 1
        });
        let used = 1 + usize::from(second.is_some());
        for row in next + used..n {
            clear_site(tableau, row, column, next, second);
        }
        next += used;
    }

    let mut active: Vec<_> = (0..n).rev().collect();
    for column in (0..n).rev() {
        let Some(first) = active
            .iter()
            .position(|&row| site(tableau, row, column) != Pauli::I)
        else {
            continue;
        };
        let first_row = active[first];
        let second = (first + 1..active.len()).find(|&index| {
            let p = site(tableau, active[index], column);
            p != Pauli::I && p != site(tableau, first_row, column)
        });
        let second_row = second.map(|index| active[index]);
        for &row in &active[first + 1..] {
            if Some(row) != second_row {
                clear_site(tableau, row, column, first_row, second_row);
            }
        }
        if let Some(index) = second {
            active.remove(index);
        }
        active.remove(first);
    }
}

fn clear_site(
    tableau: &mut Tableau,
    row: usize,
    column: usize,
    first: usize,
    second: Option<usize>,
) {
    let p = site(tableau, row, column);
    if p == Pauli::I {
        return;
    }
    if p == site(tableau, first, column) {
        multiply(tableau, first, row);
    } else if let Some(second) = second {
        if p != site(tableau, second, column) {
            // The third nonidentity Pauli needs both pivots to cancel.
            multiply(tableau, first, row);
        }
        multiply(tableau, second, row);
    }
}
