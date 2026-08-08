// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::complex::Complex64;
use ppvm_traits_2::{Clifford, Pauli};

use crate::{Bitstring, GeneralizedTableau, RowStorage};

fn amplitudes_equal<I: Bitstring>(
    left: &GeneralizedTableau<impl RowStorage, I, impl Sized>,
    right: &GeneralizedTableau<impl RowStorage, I, impl Sized>,
) -> bool {
    if left.coefficients.len() != right.coefficients.len() {
        return false;
    }
    let cutoff_sq = left.coefficient_threshold * left.coefficient_threshold;
    left.coefficients.iter().all(|(value, index)| {
        let delta: Complex64 = *value - right.coefficients.get(index);
        delta.norm_sqr() < cutoff_sq
    })
}

pub(crate) fn structurally_equal<A: RowStorage, I: Bitstring, H>(
    left: &GeneralizedTableau<A, I, H>,
    right: &GeneralizedTableau<A, I, H>,
) -> bool {
    left.is_lost == right.is_lost
        && left.tableau.rows().eq(right.tableau.rows())
        && amplitudes_equal(left, right)
}

#[derive(Clone, Copy)]
pub(crate) enum Mutation {
    Pauli {
        pauli: Pauli,
        qubit: usize,
    },
    Pauli2 {
        first: Pauli,
        second: Pauli,
        qubit0: usize,
        qubit1: usize,
    },
    Loss {
        qubit: usize,
    },
    Loss2 {
        qubit0: usize,
        qubit1: usize,
    },
}

pub(crate) fn apply_mutation<A: RowStorage, I: Bitstring, H>(
    tab: &mut GeneralizedTableau<A, I, H>,
    mutation: Mutation,
) {
    match mutation {
        Mutation::Pauli { pauli, qubit } => match pauli {
            Pauli::X => tab.x(qubit),
            Pauli::Y => tab.y(qubit),
            Pauli::Z => tab.z(qubit),
            Pauli::I => {}
        },
        Mutation::Loss { qubit } => tab.is_lost[qubit] = true,
        Mutation::Pauli2 {
            first,
            second,
            qubit0,
            qubit1,
        } => {
            apply_mutation(
                tab,
                Mutation::Pauli {
                    pauli: first,
                    qubit: qubit0,
                },
            );
            apply_mutation(
                tab,
                Mutation::Pauli {
                    pauli: second,
                    qubit: qubit1,
                },
            );
        }
        Mutation::Loss2 { qubit0, qubit1 } => {
            tab.is_lost[qubit0] = true;
            tab.is_lost[qubit1] = true;
        }
    }
}

pub(crate) fn structurally_equal_mutated<A: RowStorage, I: Bitstring, H>(
    existing: &GeneralizedTableau<A, I, H>,
    parent: &GeneralizedTableau<A, I, H>,
    mutation: Mutation,
) -> bool {
    let loss_equal = match mutation {
        Mutation::Loss { qubit } => existing
            .is_lost
            .iter()
            .zip(&parent.is_lost)
            .enumerate()
            .all(|(i, (&actual, &old))| actual == (old || i == qubit)),
        Mutation::Loss2 { qubit0, qubit1 } => existing
            .is_lost
            .iter()
            .zip(&parent.is_lost)
            .enumerate()
            .all(|(i, (&actual, &old))| actual == (old || i == qubit0 || i == qubit1)),
        Mutation::Pauli { .. } | Mutation::Pauli2 { .. } => existing.is_lost == parent.is_lost,
    };
    if !loss_equal || !amplitudes_equal(existing, parent) {
        return false;
    }
    existing
        .tableau
        .rows()
        .zip(parent.tableau.rows())
        .enumerate()
        .all(|(row, ((ex, ez, ep), (px, pz, pp)))| {
            if ex != px || ez != pz {
                return false;
            }
            let flip = match mutation {
                Mutation::Loss { .. } | Mutation::Loss2 { .. } => false,
                Mutation::Pauli { pauli, qubit } => {
                    pauli_flip(pauli, parent.tableau.row_site(row, qubit))
                }
                Mutation::Pauli2 {
                    first,
                    second,
                    qubit0,
                    qubit1,
                } => {
                    pauli_flip(first, parent.tableau.row_site(row, qubit0))
                        ^ pauli_flip(second, parent.tableau.row_site(row, qubit1))
                }
            };
            ep == pp ^ ((flip as u8) << 1)
        })
}

fn pauli_flip(applied: Pauli, row_site: Pauli) -> bool {
    matches!(
        (applied, row_site),
        (Pauli::X, Pauli::Y | Pauli::Z)
            | (Pauli::Y, Pauli::X | Pauli::Z)
            | (Pauli::Z, Pauli::X | Pauli::Y)
    )
}
