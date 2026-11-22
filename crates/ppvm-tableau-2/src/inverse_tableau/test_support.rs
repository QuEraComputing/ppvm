// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Helpers shared by the tableau's unit tests.

use ppvm_pauli_kernels::Bits;
use ppvm_traits_2::Clifford;

use super::{InverseTableau, RowRef};

pub(super) const SIZES: [usize; 7] = [1, 2, 7, 63, 64, 65, 130];

/// Whether two rows anticommute: odd symplectic product.
fn anticommute(a: RowRef, b: RowRef) -> bool {
    let mut parity = 0;
    for i in 0..a.x.len() {
        parity ^= ((a.x[i] & b.z[i]) ^ (a.z[i] & b.x[i])).popcount() & 1;
    }
    parity == 1
}

/// Conjugation preserves commutation: `ix_a`, `iz_b` anticommute iff `a == b`.
pub(super) fn assert_canonical_commutation(t: &InverseTableau) {
    let n = t.num_qubits();
    for a in 0..n {
        for b in 0..n {
            assert_eq!(anticommute(t.ix(a), t.iz(b)), a == b);
            assert!(!anticommute(t.ix(a), t.ix(b)));
            assert!(!anticommute(t.iz(a), t.iz(b)));
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) enum Gate {
    One(usize, usize),
    Two(usize, usize, usize),
}

/// `gate` as a physical gate (or on any `Clifford`).
pub(super) fn apply<T: Clifford>(t: &mut T, gate: Gate) {
    match gate {
        Gate::One(g, q) => match g {
            0 => t.x(q),
            1 => t.y(q),
            2 => t.z(q),
            3 => t.h(q),
            4 => t.s(q),
            5 => t.s_dag(q),
            6 => t.sqrt_x(q),
            7 => t.sqrt_x_dag(q),
            8 => t.sqrt_y(q),
            _ => t.sqrt_y_dag(q),
        },
        Gate::Two(g, a, b) => match g {
            0 => t.cnot(a, b),
            1 => t.cz(a, b),
            _ => t.cy(a, b),
        },
    }
}

type Rule1 = fn(bool, bool) -> ([bool; 2], bool);
type Rule2 = fn([bool; 4]) -> ([bool; 4], bool);

/// `gate` as a change of basis on virtual qubits, `C ← C·V`.
pub(super) fn change_basis(t: &mut InverseTableau, gate: Gate) {
    use ppvm_pauli_kernels::heisenberg as rules;
    let one: [Rule1; 10] = [
        rules::x,
        rules::y,
        rules::z,
        rules::h,
        rules::s,
        rules::s_dag,
        rules::sqrt_x,
        rules::sqrt_x_dag,
        rules::sqrt_y,
        rules::sqrt_y_dag,
    ];
    let two: [Rule2; 3] = [rules::cnot, rules::cz, rules::cy];
    match gate {
        Gate::One(g, v) => t.change_basis1(v, one[g]),
        Gate::Two(g, a, b) => t.change_basis2(a, b, two[g]),
    }
}

/// A random circuit from a xorshift seed.
pub(super) fn random_circuit(n: usize, len: usize, mut seed: u64) -> Vec<Gate> {
    let mut next = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    };
    (0..len)
        .map(|_| {
            if n < 2 || next() % 3 == 0 {
                Gate::One(next() % 10, next() % n)
            } else {
                let a = next() % n;
                let b = (a + 1 + next() % (n - 1)) % n;
                Gate::Two(next() % 3, a, b)
            }
        })
        .collect()
}
