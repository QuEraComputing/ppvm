// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::complex::Complex64;
use ppvm_pauli_sum::config::fxhash::ByteF64;
use ppvm_tableau::prelude::*;

type TestTableau = GeneralizedTableau<ByteF64<1>>;

fn ccz_via_t(tab: &mut TestTableau, a: usize, b: usize, c: usize) {
    tab.t(a);
    tab.t(b);
    tab.t(c);
    tab.cnot(a, b);
    tab.t_dag(b);
    tab.cnot(a, b);
    tab.cnot(a, c);
    tab.t_dag(c);
    tab.cnot(b, c);
    tab.t(c);
    tab.cnot(a, c);
    tab.t_dag(c);
    tab.cnot(b, c);
}

#[test]
fn ccz_matches_seven_t_gates() {
    let mut native = TestTableau::new(3, 1e-12);
    for q in 0..3 {
        native.h(q);
    }
    native.rx(0, 0.31);
    native.ry(1, -0.27);
    native.rz(2, 0.61);
    native.cnot(0, 2);
    native.s(1);
    let mut naive = native.clone();

    native.ccz(2, 0, 1);
    ccz_via_t(&mut naive, 2, 0, 1);

    // The decomposition restores the same tableau, so compare all amplitudes.
    assert_eq!(native.tableau.data, naive.tableau.data);
    for idx in 0..8 {
        let actual = SparseVector::get(&native.coefficients, &idx);
        let expected = SparseVector::get(&naive.coefficients, &idx);
        assert!((actual - expected).norm() < 1e-12);
    }
}

#[test]
fn ccz_truncates_after_merging_branches() {
    let mut tab = TestTableau::new(3, 1e-12);
    for q in 0..3 {
        tab.h(q);
    }
    ccz_via_t(&mut tab, 0, 1, 2);

    // The inverse must recover |+++>, including contributions below the cutoff.
    tab.coefficient_threshold = 0.3;
    tab.ccz(0, 1, 2);
    assert_eq!(tab.coefficients.len(), 1);
    assert!((SparseVector::get(&tab.coefficients, &0) - Complex64::new(1.0, 0.0)).norm() < 1e-12);
}

#[test]
fn ccz_skips_lost_or_leaked_operands() {
    for status in [QubitStatus::Lost, QubitStatus::Leaked] {
        for inactive in 0..3 {
            let mut tab = TestTableau::new(3, 1e-12);
            for q in 0..3 {
                tab.h(q);
            }
            tab.qubit_status[inactive] = status;
            let before = tab.clone();
            tab.ccz(0, 1, 2);
            assert_eq!(tab.tableau.data, before.tableau.data);
            assert_eq!(tab.coefficients, before.coefficients);
            assert_eq!(tab.qubit_status, before.qubit_status);
        }
    }
}
