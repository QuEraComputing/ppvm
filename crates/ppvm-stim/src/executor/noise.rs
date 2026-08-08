// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use itertools::Itertools;
use num::Integer;
use stim_parser::prelude::{NoiseName, NoiseOp};

use super::StimTableau;

pub(super) fn execute<T: StimTableau>(op: &NoiseOp, tab: &mut T) {
    let NoiseOp {
        name,
        targets,
        args,
        ..
    } = op;
    match name {
        NoiseName::Depolarize1 => {
            for &q in targets {
                tab.depolarize1(q, args[0]);
            }
        }
        NoiseName::Depolarize2 => {
            for (a, b) in targets.iter().copied().tuples() {
                tab.depolarize2(a, b, args[0]);
            }
        }
        NoiseName::PauliChannel1 => {
            let p = [args[0], args[1], args[2]];
            for &q in targets {
                tab.pauli_error(q, p);
            }
        }
        NoiseName::PauliChannel2 => {
            debug_assert!(targets.len().is_even());
            let p = std::array::from_fn(|i| args[i]);
            for (a, b) in targets.iter().copied().tuples() {
                tab.two_qubit_pauli_error(a, b, p);
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
            for &q in targets {
                tab.pauli_error(q, p);
            }
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
