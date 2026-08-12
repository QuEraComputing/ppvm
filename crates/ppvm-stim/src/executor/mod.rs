// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod adapter;
mod api;
mod gates;
mod helpers;
mod measure;
mod noise;
mod sample;

pub use adapter::StimTableau;
#[doc(hidden)]
pub use adapter::{SelectedBackend, TableauType};
pub use api::{BackendTableau, execute, execute_with_rng};
pub use sample::{
    sample, sample_serial, sample_serial_validated, sample_serial_validated_with_rng,
    sample_serial_with_rng, sample_validated, sample_validated_with_rng, sample_with_rng,
};
#[cfg(feature = "rayon")]
pub use sample::{
    sample_parallel, sample_parallel_validated, sample_parallel_validated_with_rng,
    sample_parallel_with_rng,
};

use stim_parser::prelude::{Axis, ExtendedInstruction};

/// Dispatch validated instructions, appending measurement bits to `results`.
pub fn execute_validated<C, I, S>(
    instructions: &[ExtendedInstruction],
    tab: &mut BackendTableau<C, I, S>,
    results: &mut Vec<Option<bool>>,
) where
    SelectedBackend: TableauType<C, I, S>,
{
    let mut rng = rand::make_rng::<rand::rngs::SmallRng>();
    execute_validated_with_rng(instructions, tab, results, &mut rng);
}

/// Execute an already validated instruction stream with caller-supplied randomness.
pub fn execute_validated_with_rng<C, I, S, R: rand::Rng + ?Sized>(
    instructions: &[ExtendedInstruction],
    tab: &mut BackendTableau<C, I, S>,
    results: &mut Vec<Option<bool>>,
    rng: &mut R,
) where
    SelectedBackend: TableauType<C, I, S>,
{
    dispatch(instructions, tab, results, rng);
}

fn dispatch<T: StimTableau, R: rand::Rng + ?Sized>(
    instructions: &[ExtendedInstruction],
    tab: &mut T,
    results: &mut Vec<Option<bool>>,
    rng: &mut R,
) {
    for instruction in instructions {
        match instruction {
            ExtendedInstruction::Gate(op) => gates::execute(op, tab, rng),
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
            ExtendedInstruction::Noise(op) => noise::execute(op, tab, rng),
            ExtendedInstruction::Loss { p, targets, .. } => {
                for &q in targets {
                    tab.loss_channel(q, *p, rng);
                }
            }
            ExtendedInstruction::CorrelatedLoss { ps, targets, .. } => {
                for &(a, b) in targets {
                    tab.correlated_loss_channel(a, b, *ps, rng);
                }
            }
            ExtendedInstruction::Measure(op) => measure::execute(op, tab, results, rng),
            ExtendedInstruction::MPad { bits, prob, .. } => {
                let noise = prob.unwrap_or(0.0);
                for &bit in bits {
                    let recorded = Some(tab.flip_with_prob(bit, noise, rng));
                    tab.append_measurement_record(recorded);
                    results.push(recorded);
                }
            }
            ExtendedInstruction::Mpp(op) => measure::execute_mpp(op, tab, results, rng),
            ExtendedInstruction::Annotation(_) => {}
            ExtendedInstruction::Repeat { count, body, .. } => {
                for _ in 0..*count {
                    dispatch(body, tab, results, rng);
                }
            }
        }
    }
}
