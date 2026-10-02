// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use itertools::Itertools;
use num::Integer;
use stim_parser::prelude::{NoiseName, NoiseOp};

use super::StimTableau;

pub(super) fn execute<T: StimTableau, R: rand::Rng + ?Sized>(
    op: &NoiseOp,
    tab: &mut T,
    rng: &mut R,
) {
    let NoiseOp {
        name,
        targets,
        args,
        ..
    } = op;
    match name {
        NoiseName::Depolarize1 => tab.depolarize1_many(targets, args[0], rng),
        NoiseName::Depolarize2 => tab.depolarize2_many(targets, args[0], rng),
        NoiseName::PauliChannel1 => tab.pauli_error_many(targets, [args[0], args[1], args[2]], rng),
        NoiseName::PauliChannel2 => {
            debug_assert!(targets.len().is_even());
            let p = std::array::from_fn(|i| args[i]);
            for (a, b) in targets.iter().copied().tuples() {
                tab.two_qubit_pauli_error(a, b, p, rng);
            }
        }
        NoiseName::XError | NoiseName::YError | NoiseName::ZError => {
            let zero = 0.0;
            let p = match name {
                NoiseName::XError => [args[0], zero, zero],
                NoiseName::YError => [zero, args[0], zero],
                NoiseName::ZError => [zero, zero, args[0]],
                _ => unreachable!(),
            };
            tab.pauli_error_many(targets, p, rng);
        }
        NoiseName::IError
        | NoiseName::HeraldedErase
        | NoiseName::HeraldedPauliChannel1
        | NoiseName::CorrelatedError
        | NoiseName::ElseCorrelatedError => {
            unreachable!("unsupported noise {name:?} should have been rejected by validate")
        }
    }
}
