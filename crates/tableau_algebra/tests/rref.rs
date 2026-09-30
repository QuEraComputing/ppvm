// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod common;

use common::{assert_frame, ghz, stabilizers};
use ppvm_tableau_2::Tableau;
use ppvm_traits_2::Pauli;
use tableau_algebra::canonicalize_rref;

#[test]
fn ghz_partial_trace_keeps_the_untouched_correlation() {
    let mut tableau = ghz(3);
    assert_eq!(canonicalize_rref(&mut tableau, &[0]), 1);
    assert_eq!(stabilizers(&tableau), ["+IZZ", "+ZZI", "+XXX"]);
    assert_frame(&tableau);
}

#[test]
fn selections_across_word_boundaries() {
    for n in [1, 31, 32, 33, 63, 64, 65, 127, 128, 129] {
        let mut tableau = ghz(n);
        let columns = [n - 1, 0, n / 2, n - 1];
        let remaining = canonicalize_rref(&mut tableau, &columns);
        for row in 0..remaining {
            for q in columns {
                assert_eq!(tableau.row_site(n + row, q), Pauli::I);
            }
        }
        assert_frame(&tableau);

        let all = (0..n).rev().collect::<Vec<_>>();
        assert_eq!(canonicalize_rref(&mut tableau, &all), 0);
        assert_frame(&tableau);
    }
}

#[test]
fn empty_selection_leaves_the_frame_unchanged() {
    for n in [0, 3] {
        let mut tableau = ghz(n);
        let before = tableau.clone();
        assert_eq!(canonicalize_rref(&mut tableau, &[]), n);
        assert_eq!(tableau, before);
    }
}

#[test]
fn invalid_column_panics_before_mutation() {
    let mut tableau = Tableau::new(3);
    let before = tableau.clone();
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        canonicalize_rref(&mut tableau, &[0, 3]);
    }));
    assert!(result.is_err());
    assert_eq!(tableau, before);
}
