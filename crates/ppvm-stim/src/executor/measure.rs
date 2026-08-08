// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use stim_parser::prelude::{MeasureName, MeasureOp, MppOp, PauliAxis};

use super::StimTableau;
use super::helpers::measure_reset_z;

pub(super) fn execute<T: StimTableau>(
    op: &MeasureOp,
    tab: &mut T,
    results: &mut Vec<Option<bool>>,
) {
    let MeasureOp {
        name,
        args,
        targets,
        ..
    } = op;
    let noise = args.first().copied().unwrap_or(0.0);
    match name {
        MeasureName::M | MeasureName::MZ => {
            if noise > 0.0 {
                for &q in targets {
                    results.push(tab.measure_noisy(q, noise));
                }
            } else {
                results.extend(tab.measure_many(targets));
            }
        }
        MeasureName::MR => {
            for &q in targets {
                results.push(measure_reset_z(tab, q, noise));
            }
        }
        MeasureName::MX => {
            for &q in targets {
                tab.h(q);
                results.push(tab.measure_noisy(q, noise));
                tab.h(q);
            }
        }
        MeasureName::MY => {
            for &q in targets {
                tab.s_dag(q);
                tab.h(q);
                results.push(tab.measure_noisy(q, noise));
                tab.h(q);
                tab.s(q);
            }
        }
        MeasureName::MRX => {
            for &q in targets {
                tab.h(q);
                results.push(measure_reset_z(tab, q, noise));
                tab.h(q);
            }
        }
        MeasureName::MRY => {
            for &q in targets {
                tab.s_dag(q);
                tab.h(q);
                results.push(measure_reset_z(tab, q, noise));
                tab.h(q);
                tab.s(q);
            }
        }
        MeasureName::MXX | MeasureName::MYY | MeasureName::MZZ | MeasureName::MPP => {
            unreachable!("unsupported measure {name:?} should have been rejected by validate")
        }
    }
}

pub(super) fn execute_mpp<T: StimTableau>(
    op: &MppOp,
    tab: &mut T,
    results: &mut Vec<Option<bool>>,
) {
    let noise = op.args.first().copied().unwrap_or(0.0);
    for product in &op.products {
        for factor in product {
            basis_to_z(tab, factor.axis, factor.qubit);
        }
        let q0 = product[0].qubit;
        for factor in &product[1..] {
            tab.cnot(factor.qubit, q0);
        }
        results.push(tab.measure_noisy(q0, noise));
        for factor in product[1..].iter().rev() {
            tab.cnot(factor.qubit, q0);
        }
        for factor in product {
            basis_from_z(tab, factor.axis, factor.qubit);
        }
    }
}

fn basis_to_z<T: StimTableau>(tab: &mut T, axis: PauliAxis, q: usize) {
    match axis {
        PauliAxis::X => tab.h(q),
        PauliAxis::Y => {
            tab.s_dag(q);
            tab.h(q);
        }
        PauliAxis::Z => {}
    }
}

fn basis_from_z<T: StimTableau>(tab: &mut T, axis: PauliAxis, q: usize) {
    match axis {
        PauliAxis::X => tab.h(q),
        PauliAxis::Y => {
            tab.h(q);
            tab.s(q);
        }
        PauliAxis::Z => {}
    }
}
