// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use stim_parser::prelude::{GateName, GateOp, Target};

use super::StimTableau;
use super::helpers::{has_record_control, qubit, qubit_pairs, qubits, record_bit};

pub(super) fn execute<T: StimTableau, R: rand::Rng + ?Sized>(
    op: &GateOp,
    tab: &mut T,
    rng: &mut R,
) {
    let GateOp { name, targets, .. } = op;
    match name {
        GateName::Reset | GateName::ResetZ => {
            for &target in targets {
                tab.reset(qubit(target), rng);
            }
        }
        GateName::ResetX => {
            for &target in targets {
                let q = qubit(target);
                tab.reset(q, rng);
                tab.h(q);
            }
        }
        GateName::ResetY => {
            for &target in targets {
                let q = qubit(target);
                tab.reset(q, rng);
                tab.h(q);
                tab.s(q);
            }
        }
        GateName::X => tab.x_many(&qubits(targets)),
        GateName::Y => tab.y_many(&qubits(targets)),
        GateName::Z => tab.z_many(&qubits(targets)),
        GateName::H | GateName::HXZ => tab.h_many(&qubits(targets)),
        GateName::S | GateName::SqrtZ => tab.s_many(&qubits(targets)),
        GateName::SDag | GateName::SqrtZDag => tab.s_dag_many(&qubits(targets)),
        GateName::SqrtX => tab.sqrt_x_many(&qubits(targets)),
        GateName::SqrtXDag => tab.sqrt_x_dag_many(&qubits(targets)),
        GateName::SqrtY => tab.sqrt_y_many(&qubits(targets)),
        GateName::SqrtYDag => tab.sqrt_y_dag_many(&qubits(targets)),
        GateName::Identity => {}
        GateName::CX | GateName::ZCX | GateName::CNot => {
            controlled(targets, tab, Controlled::X);
        }
        GateName::CY | GateName::ZCY => {
            controlled(targets, tab, Controlled::Y);
        }
        GateName::CZ | GateName::ZCZ => {
            controlled(targets, tab, Controlled::Z);
        }
        GateName::Swap
        | GateName::ISwap
        | GateName::ISwapDag
        | GateName::SqrtXX
        | GateName::SqrtYY
        | GateName::SqrtZZ
        | GateName::CXSwap
        | GateName::SwapCX
        | GateName::XCX
        | GateName::XCY
        | GateName::XCZ
        | GateName::YCX
        | GateName::YCY
        | GateName::YCZ
        | GateName::CXYZ
        | GateName::CZYX
        | GateName::HXY
        | GateName::HYZ => {
            unreachable!("unsupported gate {name:?} should have been rejected by validate")
        }
        GateName::T | GateName::TDag => {
            unreachable!("T/T_DAG are lowered by interpret")
        }
    }
}

#[derive(Clone, Copy)]
enum Controlled {
    X,
    Y,
    Z,
}

fn controlled<T: StimTableau>(targets: &[Target], tab: &mut T, gate: Controlled) {
    if has_record_control(targets) {
        for pair in targets.chunks_exact(2) {
            match pair[0] {
                Target::Qubit(control) => match gate {
                    Controlled::X => tab.cnot(control, qubit(pair[1])),
                    Controlled::Y => tab.cy(control, qubit(pair[1])),
                    Controlled::Z => tab.cz(control, qubit(pair[1])),
                },
                Target::Rec(k) => {
                    if record_bit(tab.measurement_record(), k) {
                        match gate {
                            Controlled::X => tab.x(qubit(pair[1])),
                            Controlled::Y => tab.y(qubit(pair[1])),
                            Controlled::Z => tab.z(qubit(pair[1])),
                        }
                    }
                }
            }
        }
    } else {
        let pairs = qubit_pairs(targets);
        // Retain each backend's fused batch path.
        match gate {
            Controlled::X => tab.cnot_many(&pairs),
            Controlled::Y => tab.cy_many(&pairs),
            Controlled::Z => tab.cz_many(&pairs),
        }
    }
}
