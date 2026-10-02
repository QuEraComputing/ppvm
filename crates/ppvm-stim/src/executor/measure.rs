// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use stim_parser::prelude::{MeasureName, MeasureOp, MppOp, PauliAxis};

use super::StimTableau;
use super::helpers::has_repeats;

pub(super) fn execute<T: StimTableau, R: rand::Rng + ?Sized>(
    op: &MeasureOp,
    tab: &mut T,
    results: &mut Vec<Option<bool>>,
    rng: &mut R,
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
                tab.measure_noisy_many(targets, noise, rng, results);
            } else {
                results.extend(tab.measure_many(targets, rng));
            }
        }
        MeasureName::MR => tab.measure_reset_many(targets, noise, rng, results),
        MeasureName::MX | MeasureName::MY => {
            let axis = if *name == MeasureName::MX {
                PauliAxis::X
            } else {
                PauliAxis::Y
            };
            measure_in_basis(tab, targets, axis, results, |tab, q, results| {
                tab.measure_noisy_many(q, noise, rng, results)
            });
        }
        MeasureName::MRX | MeasureName::MRY => {
            let axis = if *name == MeasureName::MRX {
                PauliAxis::X
            } else {
                PauliAxis::Y
            };
            measure_in_basis(tab, targets, axis, results, |tab, q, results| {
                tab.measure_reset_many(q, noise, rng, results)
            });
        }
        MeasureName::MXX | MeasureName::MYY | MeasureName::MZZ | MeasureName::MPP => {
            unreachable!("unsupported measure {name:?} should have been rejected by validate")
        }
    }
}

pub(super) fn execute_mpp<T: StimTableau, R: rand::Rng + ?Sized>(
    op: &MppOp,
    tab: &mut T,
    results: &mut Vec<Option<bool>>,
    rng: &mut R,
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
        results.push(tab.measure_noisy(q0, noise, rng));
        for factor in product[1..].iter().rev() {
            tab.cnot(factor.qubit, q0);
        }
        for factor in product {
            basis_from_z(tab, factor.axis, factor.qubit);
        }
    }
}

/// Measure `targets` along `axis` with `measure`, a Z-basis batch. Distinct
/// targets rotate onto Z all at once around one batch; a repeated target keeps
/// the per-target order, where the rotations between its measurements matter.
fn measure_in_basis<T: StimTableau>(
    tab: &mut T,
    targets: &[usize],
    axis: PauliAxis,
    results: &mut Vec<Option<bool>>,
    mut measure: impl FnMut(&mut T, &[usize], &mut Vec<Option<bool>>),
) {
    if has_repeats(targets) {
        for &q in targets {
            basis_to_z(tab, axis, q);
            measure(tab, &[q], results);
            basis_from_z(tab, axis, q);
        }
        return;
    }
    targets.iter().for_each(|&q| basis_to_z(tab, axis, q));
    measure(tab, targets, results);
    targets.iter().for_each(|&q| basis_from_z(tab, axis, q));
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
