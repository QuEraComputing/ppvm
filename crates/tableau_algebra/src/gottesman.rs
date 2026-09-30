// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_tableau_2::Tableau;
use ppvm_traits_2::Clifford;

use crate::rows::{forward_pass, has_x, has_z};

/// Ranks and qubit permutations from [`canonicalize_gott`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GottesmanForm {
    /// Number of X pivots.
    pub x_rank: usize,
    /// Number of Z pivots in the second pass (not the total rank).
    pub z_rank: usize,
    /// New position to preceding position, after the X pass.
    pub x_permutation: Vec<usize>,
    /// New position to preceding position, after the Z pass.
    /// The final qubit `q` came from `x_permutation[z_permutation[q]]`.
    pub z_permutation: Vec<usize>,
}

/// Put the stabilizers in Gottesman's standard form, reordering qubits.
///
/// > **Warning:** Requires a valid full-rank frame: `n` destabilizers followed by
/// > their `n` paired stabilizers on `n` qubits (`2n × n` Pauli entries).
/// > Stabilizer-only and rank-deficient tableaux are unsupported.
///
/// Like QuantumClifford's `canonicalize_gott!`, first reduce X and move its
/// pivot columns to the left. Then reduce Z on the remaining columns and
/// move those pivots next. The pivot blocks are identity matrices. This form
/// is used in stabilizer-code algebra, including QuantumClifford's logical
/// operator construction (Gottesman, 1997).
///
/// Both halves of the frame follow the qubit permutations. These are physical
/// qubit reorderings, so use the returned permutations to track qubit labels.
///
/// ```
/// use ppvm_tableau_2::Tableau;
/// use ppvm_traits_2::Clifford;
/// use tableau_algebra::canonicalize_gott;
///
/// let mut tableau = Tableau::new(3);
/// tableau.h(2);
/// let form = canonicalize_gott(&mut tableau);
/// assert_eq!((form.x_rank, form.z_rank), (1, 2));
/// assert_eq!(form.x_permutation, [2, 0, 1]);
/// assert_eq!(form.z_permutation, [0, 1, 2]);
/// // Stabilizers are now +XII, +IZI, +IIZ.
/// ```
pub fn canonicalize_gott(tableau: &mut Tableau) -> GottesmanForm {
    let n = tableau.n_qubits();
    let x_pivots = forward_pass(tableau, 0, 0..n, has_x);
    let x_rank = x_pivots.len();
    let x_permutation = pivot_permutation(n, 0, &x_pivots);
    permute_qubits(tableau, &x_permutation);

    let z_pivots = forward_pass(tableau, x_rank, x_rank..n, has_z);
    let z_permutation = pivot_permutation(n, x_rank, &z_pivots);
    permute_qubits(tableau, &z_permutation);
    GottesmanForm {
        x_rank,
        z_rank: z_pivots.len(),
        x_permutation,
        z_permutation,
    }
}

fn pivot_permutation(n: usize, start: usize, pivots: &[usize]) -> Vec<usize> {
    (0..start)
        .chain(pivots.iter().copied())
        .chain((start..n).filter(|q| !pivots.contains(q)))
        .collect()
}

fn permute_qubits(tableau: &mut Tableau, permutation: &[usize]) {
    let mut current: Vec<_> = (0..tableau.n_qubits()).collect();
    for (new, &old) in permutation.iter().enumerate() {
        let from = current.iter().position(|&q| q == old).unwrap();
        if new != from {
            // SWAP via the public Clifford interface, including both halves.
            tableau.cnot(new, from);
            tableau.cnot(from, new);
            tableau.cnot(new, from);
            current.swap(new, from);
        }
    }
}
