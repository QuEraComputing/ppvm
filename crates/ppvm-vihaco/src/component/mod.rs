// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod backend;
mod dispatch;
mod pauli_sum;
mod tableau;

use eyre::Result;
use vihaco::{Effects, component, observe};
use vihaco_circuit_isa::{CircuitEffect, CircuitInstruction, CircuitMessage};

use crate::device_info::PPVMDeviceInfo;
use crate::measurements::CircuitOutcomeEffect;

pub use pauli_sum::{
    LossyPauliSumCircuit, LossyPauliSumExecutor, PauliSumCircuit, PauliSumExecutor,
};
pub use tableau::{CircuitExecutor, TableauCircuit};

/// Largest qubit count supported by every circuit backend.
pub const MAX_QUBITS: usize = 2048;

/// Stable circuit-level backend selector. Its public behavior is independent of
/// whether the crate was compiled with `legacy` or `traits-2`.
pub enum Circuit {
    Tableau(TableauCircuit),
    PauliSum(PauliSumCircuit),
    LossyPauliSum(LossyPauliSumCircuit),
}

#[component(
    instruction = CircuitInstruction,
    message = CircuitMessage,
    effect = CircuitOutcomeEffect
)]
impl Circuit {
    pub fn tableau(info: &PPVMDeviceInfo) -> Result<Self> {
        Ok(Self::Tableau(TableauCircuit::new(
            info.n_qubits,
            info.coefficient_threshold,
        )?))
    }

    pub fn tableau_with_seed(info: &PPVMDeviceInfo, seed: u64) -> Result<Self> {
        Ok(Self::Tableau(TableauCircuit::new_with_seed(
            info.n_qubits,
            info.coefficient_threshold,
            seed,
        )?))
    }

    pub fn paulisum(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        Ok(Self::PauliSum(PauliSumCircuit::new(info, terms)?))
    }

    pub fn lossy_paulisum(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        Ok(Self::LossyPauliSum(LossyPauliSumCircuit::new(info, terms)?))
    }

    fn execute(
        &mut self,
        inst: CircuitInstruction,
        msg: CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        self.execute_instruction(&inst, &msg)
    }

    fn execute_instruction(
        &mut self,
        inst: &CircuitInstruction,
        msg: &CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        match self {
            Self::Tableau(circuit) => circuit.execute_instruction(inst, msg),
            Self::PauliSum(circuit) => circuit.execute_instruction(inst, msg),
            Self::LossyPauliSum(circuit) => circuit.execute_instruction(inst, msg),
        }
    }

    pub fn state_string(&self) -> String {
        match self {
            Self::Tableau(circuit) => circuit.state_string(),
            Self::PauliSum(circuit) => circuit.state_string(),
            Self::LossyPauliSum(circuit) => circuit.state_string(),
        }
    }
}

#[observe(CircuitEffect, effect = CircuitOutcomeEffect)]
impl Circuit {
    fn observe_circuit_effect(
        &mut self,
        effect: &CircuitEffect,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        self.execute_instruction(&effect.inst, &effect.msg)
    }
}

impl vihaco::Reset for Circuit {
    fn reset(&mut self) {
        match self {
            Self::Tableau(circuit) => circuit.reset(),
            Self::PauliSum(circuit) => circuit.reset(),
            Self::LossyPauliSum(circuit) => circuit.reset(),
        }
    }
}

impl Default for Circuit {
    fn default() -> Self {
        Self::tableau(&PPVMDeviceInfo::default()).expect("0-qubit tableau is always constructible")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn info(n_qubits: usize) -> PPVMDeviceInfo {
        PPVMDeviceInfo {
            n_qubits,
            ..Default::default()
        }
    }

    #[test]
    fn every_constructor_rejects_more_than_maximum() {
        let info = info(MAX_QUBITS + 1);
        assert!(Circuit::tableau(&info).is_err());
        assert!(Circuit::tableau_with_seed(&info, 0).is_err());
        assert!(Circuit::paulisum(&info, &[]).is_err());
        assert!(Circuit::lossy_paulisum(&info, &[]).is_err());
    }

    #[test]
    fn boundary_width_constructs() {
        assert!(Circuit::tableau(&info(MAX_QUBITS)).is_ok());
    }
}
