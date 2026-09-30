// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_tableau_2::Tableau;

use crate::rows::{forward_pass, has_x, has_z};

/// Reduce all X columns, then all Z columns, in place.
///
/// > **Warning:** Requires a valid full-rank frame: `n` destabilizers followed by
/// > their `n` paired stabilizers on `n` qubits (`2n × n` Pauli entries).
/// > Stabilizer-only and rank-deficient tableaux are unsupported.
///
/// This is QuantumClifford's `canonicalize!` form, used there for stabilizer
/// inner products (Garcia et al., 2012). Returns `(x_rank, total_rank)`, the
/// numbers of pivots after the X pass and after both passes. For a valid
/// pure-state frame, `total_rank` equals the number of qubits.
///
/// ```
/// use ppvm_tableau_2::Tableau;
/// use ppvm_traits_2::Clifford;
/// use tableau_algebra::canonicalize;
///
/// let mut ghz = Tableau::new(3);
/// ghz.h(0);
/// ghz.cnot(0, 1);
/// ghz.cnot(1, 2);
/// assert_eq!(canonicalize(&mut ghz), (1, 3));
/// // Stabilizers are now +XXX, +ZIZ, +IZZ.
/// ```
pub fn canonicalize(tableau: &mut Tableau) -> (usize, usize) {
    let n = tableau.n_qubits();
    let x_rank = forward_pass(tableau, 0, 0..n, has_x).len();
    let z_rank = forward_pass(tableau, x_rank, 0..n, has_z).len();
    (x_rank, x_rank + z_rank)
}
