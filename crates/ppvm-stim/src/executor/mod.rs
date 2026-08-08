// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod adapter;
mod api;
mod gates;
mod helpers;
mod measure;
mod noise;

pub use adapter::StimTableau;
#[doc(hidden)]
pub use adapter::{SelectedBackend, TableauType};
pub use api::{
    BackendTableau, execute, sample, sample_serial, sample_serial_validated, sample_validated,
};
#[cfg(feature = "rayon")]
pub use api::{sample_parallel, sample_parallel_validated};

use stim_parser::prelude::{Axis, ExtendedInstruction};

/// Dispatch validated instructions, appending measurement bits to `results`.
pub fn execute_validated<C, I, S>(
    instructions: &[ExtendedInstruction],
    tab: &mut BackendTableau<C, I, S>,
    results: &mut Vec<Option<bool>>,
) where
    SelectedBackend: TableauType<C, I, S>,
{
    dispatch(instructions, tab, results);
}

fn dispatch<T: StimTableau>(
    instructions: &[ExtendedInstruction],
    tab: &mut T,
    results: &mut Vec<Option<bool>>,
) {
    for instruction in instructions {
        match instruction {
            ExtendedInstruction::Gate(op) => gates::execute(op, tab),
            ExtendedInstruction::T { targets, .. } => {
                targets.iter().for_each(|&q| tab.t(q));
            }
            ExtendedInstruction::TDag { targets, .. } => {
                targets.iter().for_each(|&q| tab.t_dag(q));
            }
            ExtendedInstruction::Rotation {
                axis,
                theta,
                targets,
                ..
            } => {
                for &q in targets {
                    match axis {
                        Axis::X => tab.rx(q, *theta),
                        Axis::Y => tab.ry(q, *theta),
                        Axis::Z => tab.rz(q, *theta),
                    }
                }
            }
            ExtendedInstruction::U3 {
                theta,
                phi,
                lambda,
                targets,
                ..
            } => {
                for &q in targets {
                    tab.u3(q, *theta, *phi, *lambda);
                }
            }
            ExtendedInstruction::Noise(op) => noise::execute(op, tab),
            ExtendedInstruction::Loss { p, targets, .. } => {
                for &q in targets {
                    tab.loss_channel(q, *p);
                }
            }
            ExtendedInstruction::CorrelatedLoss { ps, targets, .. } => {
                for &(a, b) in targets {
                    tab.correlated_loss_channel(a, b, *ps);
                }
            }
            ExtendedInstruction::Measure(op) => measure::execute(op, tab, results),
            ExtendedInstruction::MPad { bits, prob, .. } => {
                let noise = prob.unwrap_or(0.0);
                for &bit in bits {
                    let recorded = Some(tab.flip_with_prob(bit, noise));
                    tab.append_measurement_record(recorded);
                    results.push(recorded);
                }
            }
            ExtendedInstruction::Mpp(op) => measure::execute_mpp(op, tab, results),
            ExtendedInstruction::Annotation(_) => {}
            ExtendedInstruction::Repeat { count, body, .. } => {
                for _ in 0..*count {
                    dispatch(body, tab, results);
                }
            }
        }
    }
}
