// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_tableau_2::Tableau;
use ppvm_traits_2::{Clifford, Pauli};

pub fn ghz(n: usize) -> Tableau {
    let mut tableau = Tableau::new(n);
    if n > 0 {
        tableau.h(0);
        for q in 1..n {
            tableau.cnot(q - 1, q);
        }
    }
    tableau
}

pub fn stabilizers(tableau: &Tableau) -> Vec<String> {
    let n = tableau.n_qubits();
    (n..2 * n)
        .map(|row| {
            let mut result = match tableau.row_phase(row) {
                0 => "+".to_owned(),
                2 => "-".to_owned(),
                phase => panic!("nonreal phase {phase}"),
            };
            for q in 0..n {
                result.push(match tableau.row_site(row, q) {
                    Pauli::I => 'I',
                    Pauli::X => 'X',
                    Pauli::Y => 'Y',
                    Pauli::Z => 'Z',
                });
            }
            result
        })
        .collect()
}

// A canonicalized tableau must still be a valid simulator frame.
pub fn assert_frame(tableau: &Tableau) {
    let n = tableau.n_qubits();
    for a in 0..2 * n {
        assert_eq!(tableau.row_phase(a) % 2, 0);
        for b in a + 1..2 * n {
            let anticommutes = (0..n).fold(false, |parity, q| {
                let p = tableau.row_site(a, q);
                let r = tableau.row_site(b, q);
                parity ^ (p != Pauli::I && r != Pauli::I && p != r)
            });
            assert_eq!(anticommutes, a < n && b == a + n);
        }
    }
}
