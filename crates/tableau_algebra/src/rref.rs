// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_tableau_2::Tableau;

use crate::rows::{eliminate, has_x, has_z, site, swap};

/// Reduce selected qubit columns, alternating X and Z elimination.
///
/// > **Warning:** Requires a valid full-rank frame: `n` destabilizers followed by
/// > their `n` paired stabilizers on `n` qubits (`2n × n` Pauli entries).
/// > Stabilizer-only and rank-deficient tableaux are unsupported.
///
/// This is QuantumClifford's `canonicalize_rref!` form (Audenaert and Plenio,
/// 2005). Pivots go at the bottom. Returns the number of leading stabilizer
/// rows that are identity on every selected qubit. Those rows generate the
/// reduced state after tracing out the selected qubits.
///
/// Columns are zero-based and processed in the supplied order. Empty and
/// repeated selections are allowed. For a full reduction, supply `0..n`.
///
/// # Panics
///
/// Panics before changing the tableau if any column is out of bounds.
///
/// ```
/// use ppvm_tableau_2::Tableau;
/// use ppvm_traits_2::Clifford;
/// use tableau_algebra::canonicalize_rref;
///
/// let mut ghz = Tableau::new(3);
/// ghz.h(0);
/// ghz.cnot(0, 1);
/// ghz.cnot(1, 2);
/// assert_eq!(canonicalize_rref(&mut ghz, &[0]), 1);
/// // The leading stabilizer +IZZ survives tracing out qubit 0.
/// ```
pub fn canonicalize_rref(tableau: &mut Tableau, columns: &[usize]) -> usize {
    let n = tableau.n_qubits();
    assert!(columns.iter().all(|&q| q < n), "column out of bounds");
    let mut remaining = n;
    for &column in columns {
        for component in [has_x, has_z] {
            if let Some(row) = (0..remaining).find(|&row| component(site(tableau, row, column))) {
                remaining -= 1;
                swap(tableau, row, remaining);
                eliminate(tableau, remaining, column, component);
            }
        }
    }
    remaining
}
