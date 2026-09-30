// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod common;

use common::{assert_frame, ghz, stabilizers};
use ppvm_tableau_2::Tableau;
use ppvm_traits_2::Clifford;
use tableau_algebra::canonicalize;

#[test]
fn ghz_matches_quantumclifford_example() {
    let mut tableau = ghz(4);
    assert_eq!(canonicalize(&mut tableau), (1, 4));
    assert_eq!(stabilizers(&tableau), ["+XXXX", "+ZIIZ", "+IZIZ", "+IIZZ"]);
    assert_frame(&tableau);
}

#[test]
fn signed_y_rows_preserve_the_state() {
    let mut tableau = ghz(2);
    tableau.s(0);
    tableau.x(1);
    assert_eq!(stabilizers(&tableau), ["+YX", "-ZZ"]);
    canonicalize(&mut tableau);
    assert_eq!(stabilizers(&tableau), ["-XY", "-ZZ"]);
    assert_frame(&tableau);
}

#[test]
fn boundary_sizes_and_idempotence() {
    for n in [0, 1, 31, 32, 33, 63, 64, 65, 127, 128, 129] {
        let mut tableau = ghz(n);
        assert_eq!(canonicalize(&mut tableau), (usize::from(n > 0), n));
        let canonical = tableau.clone();
        canonicalize(&mut tableau);
        assert_eq!(tableau, canonical);
        assert_frame(&tableau);

        let mut zero = Tableau::new(n);
        assert_eq!(canonicalize(&mut zero), (0, n));
        assert_eq!(zero, Tableau::new(n));
    }
}
