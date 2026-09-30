// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod common;

use common::{assert_frame, ghz, stabilizers};
use ppvm_tableau_2::Tableau;
use ppvm_traits_2::{Clifford, Pauli};
use tableau_algebra::{canonicalize, canonicalize_gott};

#[test]
fn moves_x_pivots_and_preserves_signs() {
    let mut tableau = Tableau::new(3);
    tableau.h(2);
    tableau.z(2);
    tableau.x(0);
    let form = canonicalize_gott(&mut tableau);
    assert_eq!((form.x_rank, form.z_rank), (1, 2));
    assert_eq!(form.x_permutation, [2, 0, 1]);
    assert_eq!(form.z_permutation, [0, 1, 2]);
    assert_eq!(stabilizers(&tableau), ["-XII", "-IZI", "+IIZ"]);
    assert_frame(&tableau);
}

#[test]
fn ghz_has_identity_pivot_blocks() {
    for n in [0, 1, 31, 32, 33, 63, 64, 65, 127, 128, 129] {
        let mut tableau = ghz(n);
        let mut expected = tableau.clone();
        canonicalize(&mut expected);
        let form = canonicalize_gott(&mut tableau);
        assert_eq!(form.x_rank + form.z_rank, n);
        assert_eq!(form.x_permutation, (0..n).collect::<Vec<_>>());
        assert_eq!(form.z_permutation, (0..n).collect::<Vec<_>>());
        for row in form.x_rank..n {
            for q in form.x_rank..n {
                assert_eq!(
                    tableau.row_site(n + row, q),
                    if row == q { Pauli::Z } else { Pauli::I }
                );
            }
        }
        assert_frame(&tableau);
        canonicalize(&mut tableau);
        assert_eq!(stabilizers(&tableau), stabilizers(&expected));
    }
}

#[test]
fn permutation_moves_a_pivot_across_word_boundaries() {
    for n in [33, 65, 129] {
        let mut tableau = Tableau::new(n);
        tableau.h(n - 1);
        let form = canonicalize_gott(&mut tableau);
        assert_eq!(form.x_rank, 1);
        assert_eq!(form.x_permutation[0], n - 1);
        assert_eq!(form.x_permutation[1..], (0..n - 1).collect::<Vec<_>>());
        for row in 0..n {
            for q in 0..n {
                let expected = if row != q {
                    Pauli::I
                } else if row == 0 {
                    Pauli::X
                } else {
                    Pauli::Z
                };
                assert_eq!(tableau.row_site(n + row, q), expected);
            }
        }
        assert_frame(&tableau);
    }
}
