// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Compile-and-run guards for the wide Python-facing bitstring tiers.

use ppvm_tableau_2::prelude::*;

fn exercise<A: RowStorage, I: Bitstring>(mut tableau: GeneralizedTableau<A, I>) {
    tableau.h(0);
    tableau.t(0);
    tableau.rz(0, 0.17);
    assert!(!tableau.coefficients.is_empty());
}

#[test]
fn bnum_indices_drive_generalized_tableaux() {
    exercise(GeneralizedTableau::<[usize; 4], U256>::new(200, 1e-12));
    exercise(GeneralizedTableau::<[usize; 8], U512>::new(400, 1e-12));
    exercise(GeneralizedTableau::<[usize; 16], U1024>::new(800, 1e-12));
    exercise(GeneralizedTableau::<[usize; 32], U2048>::new(1600, 1e-12));
}
