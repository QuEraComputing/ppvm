use crate::instruction::CircuitInstruction;
use crate::message::CircuitMessage;
use bitvec::view::BitView;
use bnum::types::{U256, U512, U1024, U2048};
use eyre::{Result, eyre};
use num::PrimInt;
use num::complex::Complex64;
use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;
use vihaco::{Effects, component, observe};

macro_rules! batch_for {
    ($tab:expr, $method:ident, $addrs:expr) => {
        for addr in $addrs { $tab.$method(*addr); }
    };
    ($tab:expr, $method:ident, $addrs:expr, $($arg:expr),+) => {
        for addr in $addrs { $tab.$method(*addr, $($arg),+); }
    };
}

pub struct CircuitExecutor<T: Config<Coeff = f64>, I: TableauIndex, C: SparseVector<Complex64, I>> {
    pub tab: GeneralizedTableau<T, I, C>,
}

pub type MeasurementResult = Vec<Option<bool>>;

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct MeasurementEffect {
    measurement_results: MeasurementResult,
}

#[component(instruction = CircuitInstruction, message = CircuitMessage, effect = MeasurementEffect)]
impl<T, I, C> CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    fn execute(
        &mut self,
        inst: CircuitInstruction,
        msg: CircuitMessage,
    ) -> Result<Effects<MeasurementEffect>> {
        self.execute_instruction(&inst, &msg)
    }

    fn execute_instruction(
        &mut self,
        inst: &CircuitInstruction,
        msg: &CircuitMessage,
    ) -> Result<Effects<MeasurementEffect>> {
        use CircuitInstruction::*;
        use CircuitMessage::*;

        match (inst, msg) {
            // Single-qubit Clifford
            (X, &Qubit(addr)) => self.tab.x(addr),
            (Y, &Qubit(addr)) => self.tab.y(addr),
            (Z, &Qubit(addr)) => self.tab.z(addr),
            (H, &Qubit(addr)) => self.tab.h(addr),
            (S, &Qubit(addr)) => self.tab.s(addr),
            (SAdj, &Qubit(addr)) => self.tab.s_adj(addr),
            (SqrtX, &Qubit(addr)) => self.tab.sqrt_x(addr),
            (SqrtY, &Qubit(addr)) => self.tab.sqrt_y(addr),
            (SqrtXAdj, &Qubit(addr)) => self.tab.sqrt_x_adj(addr),
            (SqrtYAdj, &Qubit(addr)) => self.tab.sqrt_y_adj(addr),

            // Controlled gates
            (CNOT, &TwoQubit(addr0, addr1)) => self.tab.cnot(addr0, addr1),
            (CZ, &TwoQubit(addr0, addr1)) => self.tab.cz(addr0, addr1),

            // T gate
            (T, &Qubit(addr)) => self.tab.t(addr),
            (TAdj, &Qubit(addr)) => self.tab.t_adj(addr),

            // Single-qubit rotations
            (RX, &QubitAndFloat(addr, angle)) => self.tab.rx(addr, angle),
            (RY, &QubitAndFloat(addr, angle)) => self.tab.ry(addr, angle),
            (RZ, &QubitAndFloat(addr, angle)) => self.tab.rz(addr, angle),

            // Two-qubit rotations
            (RXX, &TwoQubitAndFloat(addr0, addr1, angle)) => self.tab.rxx(addr0, addr1, angle),
            (RYY, &TwoQubitAndFloat(addr0, addr1, angle)) => self.tab.ryy(addr0, addr1, angle),
            (RZZ, &TwoQubitAndFloat(addr0, addr1, angle)) => self.tab.rzz(addr0, addr1, angle),

            // U3
            (U3, &QubitU3(addr, theta, phi, lam)) => self.tab.u3(addr, theta, phi, lam),

            // Measure & Reset
            (Measure, &Qubit(addr)) => {
                let outcome = self.tab.measure(addr);
                return Ok(Effects::one(MeasurementEffect {
                    measurement_results: vec![outcome],
                }));
            }
            (Reset, &Qubit(addr)) => self.tab.reset(addr),

            // Noise
            (Depolarize, &QubitAndFloat(addr, p)) => self.tab.depolarize(addr, p),
            (Depolarize2, &TwoQubitAndFloat(addr0, addr1, p)) => {
                self.tab.depolarize2(addr0, addr1, p)
            }
            (PauliError, QubitAndFloatArr3(addr0, ps)) => self.tab.pauli_error(*addr0, *ps),
            (TwoQubitPauliError, TwoQubitAndFloatArr15(addr0, addr1, ps)) => {
                self.tab.two_qubit_pauli_error(*addr0, *addr1, *ps)
            }

            // Loss
            (Loss, &QubitAndFloat(addr, p)) => self.tab.loss_channel(addr, p),
            (CorrelatedLoss, TwoQubitAndFloatArr3(addr0, addr1, ps)) => {
                self.tab.correlated_loss_channel(*addr0, *addr1, *ps)
            }

            /* BATCH OPERATIONS START HERE */
            // Batch: dedicated batch methods
            (SqrtX, QubitBatch(addrs)) => self.tab.sqrt_x_batch(addrs),
            (SqrtY, QubitBatch(addrs)) => self.tab.sqrt_y_batch(addrs),
            (SqrtXAdj, QubitBatch(addrs)) => self.tab.sqrt_x_adj_batch(addrs),
            (SqrtYAdj, QubitBatch(addrs)) => self.tab.sqrt_y_adj_batch(addrs),
            (H, QubitBatch(addrs)) => self.tab.h_batch(addrs),
            (CZ, TwoQubitBatch(pairs)) => self.tab.cz_batch(pairs),

            // TODO: replace things below by actual batched methods once they are available
            // Batch: single-qubit for loops
            (X, QubitBatch(addrs)) => batch_for!(self.tab, x, addrs),
            (Y, QubitBatch(addrs)) => batch_for!(self.tab, y, addrs),
            (Z, QubitBatch(addrs)) => batch_for!(self.tab, z, addrs),
            (S, QubitBatch(addrs)) => batch_for!(self.tab, s, addrs),
            (SAdj, QubitBatch(addrs)) => batch_for!(self.tab, s_adj, addrs),
            (T, QubitBatch(addrs)) => batch_for!(self.tab, t, addrs),
            (TAdj, QubitBatch(addrs)) => batch_for!(self.tab, t_adj, addrs),
            (Reset, QubitBatch(addrs)) => batch_for!(self.tab, reset, addrs),
            (RX, QubitBatchAndFloat(addrs, angle)) => batch_for!(self.tab, rx, addrs, *angle),
            (RY, QubitBatchAndFloat(addrs, angle)) => batch_for!(self.tab, ry, addrs, *angle),
            (RZ, QubitBatchAndFloat(addrs, angle)) => batch_for!(self.tab, rz, addrs, *angle),
            (Depolarize, QubitBatchAndFloat(addrs, p)) => {
                batch_for!(self.tab, depolarize, addrs, *p)
            }
            (Loss, QubitBatchAndFloat(addrs, p)) => batch_for!(self.tab, loss_channel, addrs, *p),
            (PauliError, QubitBatchAndFloatArr3(addrs, ps)) => {
                batch_for!(self.tab, pauli_error, addrs, *ps)
            }
            (U3, QubitBatchU3(addrs, theta, phi, lam)) => {
                batch_for!(self.tab, u3, addrs, *theta, *phi, *lam)
            }

            // Batch: two-qubit for loops
            (CNOT, TwoQubitBatch(pairs)) => {
                for &(a, b) in pairs {
                    self.tab.cnot(a, b);
                }
            }
            (RXX, TwoQubitBatchAndFloat(pairs, angle)) => {
                for &(a, b) in pairs {
                    self.tab.rxx(a, b, *angle);
                }
            }
            (RYY, TwoQubitBatchAndFloat(pairs, angle)) => {
                for &(a, b) in pairs {
                    self.tab.ryy(a, b, *angle);
                }
            }
            (RZZ, TwoQubitBatchAndFloat(pairs, angle)) => {
                for &(a, b) in pairs {
                    self.tab.rzz(a, b, *angle);
                }
            }
            (Depolarize2, TwoQubitBatchAndFloat(pairs, p)) => {
                for &(a, b) in pairs {
                    self.tab.depolarize2(a, b, *p);
                }
            }
            (TwoQubitPauliError, TwoQubitBatchAndFloatArr15(pairs, ps)) => {
                for &(a, b) in pairs {
                    self.tab.two_qubit_pauli_error(a, b, *ps);
                }
            }
            (CorrelatedLoss, TwoQubitBatchAndFloatArr3(pairs, ps)) => {
                for &(a, b) in pairs {
                    self.tab.correlated_loss_channel(a, b, *ps);
                }
            }

            // Batch: measure (emits per qubit)
            (Measure, QubitBatch(addrs)) => {
                let outcomes = addrs.iter().map(|&addr| self.tab.measure(addr));
                return Ok(Effects::one(MeasurementEffect {
                    measurement_results: outcomes.collect(),
                }));
            }

            // Fallback
            (inst, msg) => {
                return Err(eyre!(
                    "Invalid gate arguments {:?} for gate {:?}",
                    msg,
                    inst
                ));
            }
        };

        Ok(Effects::None)
    }
}

impl<T, I, C> vihaco::Reset for CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    fn reset(&mut self) {
        self.tab.reset_all();
    }
}

pub enum Circuit {
    Bits64(CircuitExecutor<Byte8F64<1>, usize, Vec<(Complex64, usize)>>),
    Bits128(CircuitExecutor<Byte8F64<2>, u128, Vec<(Complex64, u128)>>),
    Bits256(CircuitExecutor<Byte8F64<4>, U256, Vec<(Complex64, U256)>>),
    Bits512(CircuitExecutor<Byte8F64<8>, U512, Vec<(Complex64, U512)>>),
    Bits1024(CircuitExecutor<Byte8F64<16>, U1024, Vec<(Complex64, U1024)>>),
    Bits2048(CircuitExecutor<Byte8F64<32>, U2048, Vec<(Complex64, U2048)>>),
}

#[component(instruction = CircuitInstruction, message = CircuitMessage, effect = MeasurementEffect)]
impl Circuit {
    pub fn new(n_qubits: usize, coefficient_threshold: f64) -> Self {
        if n_qubits <= 64 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Self::Bits64(CircuitExecutor { tab })
        } else if n_qubits <= 128 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Self::Bits128(CircuitExecutor { tab })
        } else if n_qubits <= 256 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Self::Bits256(CircuitExecutor { tab })
        } else if n_qubits <= 512 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Self::Bits512(CircuitExecutor { tab })
        } else if n_qubits <= 1024 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Self::Bits1024(CircuitExecutor { tab })
        } else if n_qubits <= 2048 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Self::Bits2048(CircuitExecutor { tab })
        } else {
            panic!("No matching executor for {} qubits", n_qubits);
        }
    }

    fn execute(
        &mut self,
        inst: CircuitInstruction,
        msg: CircuitMessage,
    ) -> Result<Effects<MeasurementEffect>> {
        self.execute_instruction(&inst, &msg)
    }

    fn execute_instruction(
        &mut self,
        inst: &CircuitInstruction,
        msg: &CircuitMessage,
    ) -> Result<Effects<MeasurementEffect>> {
        match self {
            Self::Bits64(ex) => ex.execute_instruction(inst, msg),
            Self::Bits128(ex) => ex.execute_instruction(inst, msg),
            Self::Bits256(ex) => ex.execute_instruction(inst, msg),
            Self::Bits512(ex) => ex.execute_instruction(inst, msg),
            Self::Bits1024(ex) => ex.execute_instruction(inst, msg),
            Self::Bits2048(ex) => ex.execute_instruction(inst, msg),
        }
    }
}

pub struct CircuitEffect {
    pub inst: CircuitInstruction,
    pub msg: CircuitMessage,
}

#[observe(CircuitEffect, effect=MeasurementEffect)]
impl Circuit {
    fn observe_circuit_effect(
        &mut self,
        effect: &CircuitEffect,
    ) -> Result<Effects<MeasurementEffect>> {
        self.execute_instruction(&effect.inst, &effect.msg)
    }
}
