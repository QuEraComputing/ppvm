use crate::instruction::CircuitInstruction;
use crate::message::CircuitMessage;
use eyre::{Result, eyre};
use ppvm_runtime::config::fxhash::ByteF64;
use ppvm_tableau::prelude::*;
use vihaco::{Event, ExecContext, component};

pub struct Circuit<const NBytes: usize, I: TableauIndex> {
    pub tab: GeneralizedTableau<ByteF64<NBytes>, I>,
}

#[derive(Debug, Clone, Event)]
pub struct MeasurementResult {
    qubit: usize,

    /// None if lost, else 0 or 1 according to outcome
    outcome: Option<bool>,
}

#[component(instruction = CircuitInstruction, message = CircuitMessage)]
impl<const NBytes: usize, I> Circuit<NBytes, I>
where
    I: TableauIndex + Send + Sync + std::fmt::Debug,
{
    fn execute(
        &mut self,
        inst: CircuitInstruction,
        msg: CircuitMessage,
        ctx: &mut ExecContext,
    ) -> Result<()> {
        match (inst, msg) {
            // Single-qubit Clifford
            (CircuitInstruction::X, CircuitMessage::Qubit(addr)) => self.tab.x(addr),
            (CircuitInstruction::Y, CircuitMessage::Qubit(addr)) => self.tab.y(addr),
            (CircuitInstruction::Z, CircuitMessage::Qubit(addr)) => self.tab.z(addr),
            (CircuitInstruction::H, CircuitMessage::Qubit(addr)) => self.tab.h(addr),
            (CircuitInstruction::S, CircuitMessage::Qubit(addr)) => self.tab.s(addr),
            (CircuitInstruction::SAdj, CircuitMessage::Qubit(addr)) => self.tab.s_adj(addr),
            (CircuitInstruction::SqrtX, CircuitMessage::Qubit(addr)) => self.tab.sqrt_x(addr),
            (CircuitInstruction::SqrtY, CircuitMessage::Qubit(addr)) => self.tab.sqrt_y(addr),
            (CircuitInstruction::SqrtXAdj, CircuitMessage::Qubit(addr)) => {
                self.tab.sqrt_x_adj(addr)
            }
            (CircuitInstruction::SqrtYAdj, CircuitMessage::Qubit(addr)) => {
                self.tab.sqrt_y_adj(addr)
            }

            // Controlled gates
            (CircuitInstruction::CNOT, CircuitMessage::TwoQubit(addr0, addr1)) => {
                self.tab.cnot(addr0, addr1)
            }
            (CircuitInstruction::CZ, CircuitMessage::TwoQubit(addr0, addr1)) => {
                self.tab.cz(addr0, addr1)
            }

            // T gate
            (CircuitInstruction::T, CircuitMessage::Qubit(addr)) => self.tab.t(addr),
            (CircuitInstruction::TAdj, CircuitMessage::Qubit(addr)) => self.tab.t_adj(addr),

            // Single-qubit rotations
            (CircuitInstruction::RX, CircuitMessage::QubitAndFloat(addr, angle)) => {
                self.tab.rx(addr, angle)
            }
            (CircuitInstruction::RY, CircuitMessage::QubitAndFloat(addr, angle)) => {
                self.tab.ry(addr, angle)
            }
            (CircuitInstruction::RZ, CircuitMessage::QubitAndFloat(addr, angle)) => {
                self.tab.rz(addr, angle)
            }

            // Two-qubit rotations
            (CircuitInstruction::RXX, CircuitMessage::TwoQubitAndFloat(addr0, addr1, angle)) => {
                self.tab.rxx(addr0, addr1, angle)
            }
            (CircuitInstruction::RYY, CircuitMessage::TwoQubitAndFloat(addr0, addr1, angle)) => {
                self.tab.ryy(addr0, addr1, angle)
            }
            (CircuitInstruction::RZZ, CircuitMessage::TwoQubitAndFloat(addr0, addr1, angle)) => {
                self.tab.rzz(addr0, addr1, angle)
            }

            // U3
            (CircuitInstruction::U3, CircuitMessage::QubitU3(addr, theta, phi, lam)) => {
                self.tab.u3(addr, theta, phi, lam)
            }

            // Measure & Reset
            (CircuitInstruction::Measure, CircuitMessage::Qubit(addr)) => {
                let outcome = self.tab.measure(addr);
                ctx.emit(MeasurementResult {
                    qubit: addr,
                    outcome: outcome,
                });
            }
            (CircuitInstruction::Reset, CircuitMessage::Qubit(addr)) => self.tab.reset(addr),

            // Noise
            (CircuitInstruction::Depolarize, CircuitMessage::QubitAndFloat(addr, p)) => {
                self.tab.depolarize(addr, p)
            }
            (
                CircuitInstruction::Depolarize2,
                CircuitMessage::TwoQubitAndFloat(addr0, addr1, p),
            ) => self.tab.depolarize2(addr0, addr1, p),
            (CircuitInstruction::PauliError, CircuitMessage::QubitAndFloatArr3(addr0, ps)) => {
                self.tab.pauli_error(addr0, ps)
            }
            (
                CircuitInstruction::TwoQubitPauliError,
                CircuitMessage::TwoQubitAndFloatArr15(addr0, addr1, ps),
            ) => self.tab.two_qubit_pauli_error(addr0, addr1, ps),

            // Loss
            (CircuitInstruction::Loss, CircuitMessage::QubitAndFloat(addr, p)) => {
                self.tab.loss_channel(addr, p)
            }
            (
                CircuitInstruction::CorrelatedLoss,
                CircuitMessage::TwoQubitAndFloatArr3(addr0, addr1, ps),
            ) => self.tab.correlated_loss_channel(addr0, addr1, ps),

            // Fallback
            _ => {
                return Err(eyre!(
                    "Invalid gate arguments {:?} for gate {:?}",
                    msg,
                    inst
                ));
            }
        };
        Ok(())
    }
}
