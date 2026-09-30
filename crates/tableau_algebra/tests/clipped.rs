// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod common;

use common::{assert_frame, ghz, stabilizers};
use ppvm_traits_2::{Clifford, Pauli};
use tableau_algebra::{canonicalize, canonicalize_clip};

#[test]
fn ghz_is_already_clipped() {
    let mut tableau = ghz(3);
    canonicalize_clip(&mut tableau);
    assert_eq!(stabilizers(&tableau), ["+XXX", "+ZZI", "+IZZ"]);
    assert_frame(&tableau);
}

#[test]
fn endpoints_and_state_are_preserved_across_word_boundaries() {
    for n in [0, 1, 2, 7, 31, 32, 33, 63, 64, 65, 127, 128, 129] {
        let mut tableau = ghz(n);
        for q in 0..n {
            tableau.s(q);
            if q % 3 == 0 {
                tableau.h(q);
            }
            if n > 1 {
                tableau.cnot(q, (q + n - 1) % n);
            }
        }
        let mut expected = tableau.clone();
        canonicalize(&mut expected);
        canonicalize_clip(&mut tableau);

        let mut left = vec![Vec::new(); n];
        let mut right = vec![Vec::new(); n];
        for row in n..2 * n {
            let first = (0..n)
                .find(|&q| tableau.row_site(row, q) != Pauli::I)
                .unwrap();
            let last = (0..n)
                .rfind(|&q| tableau.row_site(row, q) != Pauli::I)
                .unwrap();
            left[first].push(tableau.row_site(row, first));
            right[last].push(tableau.row_site(row, last));
        }
        for q in 0..n {
            assert_eq!(left[q].len() + right[q].len(), 2);
            for endpoints in [&left[q], &right[q]] {
                if endpoints.len() == 2 {
                    assert_ne!(endpoints[0], endpoints[1]);
                }
            }
        }
        assert_frame(&tableau);
        let clipped = tableau.clone();
        canonicalize_clip(&mut tableau);
        assert_eq!(tableau, clipped);

        canonicalize(&mut tableau);
        assert_eq!(stabilizers(&tableau), stabilizers(&expected));
    }
}
