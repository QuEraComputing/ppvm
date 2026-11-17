// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::backend::CircuitBackend;
use super::buckets::{CircuitKind, LossyPauliSumBuckets, PauliSumBuckets, TableauBuckets};
use super::ops::CircuitOp;
use crate::device_info::{BackendKind, PPVMDeviceInfo};
use eyre::{Result, eyre};
use vihaco::{Execute, component};
use vihaco_circuit_isa::circuit;

/// Generic backend component - use this if you need one concrete variant in your machine
#[component]
#[instructions { circuit }]
pub struct CircuitExecutor<B: CircuitBackend> {
    pub backend: B,
}

impl<B: CircuitBackend> vihaco::Reset for CircuitExecutor<B> {
    fn reset(&mut self) {
        self.backend.reset();
    }
}

// The single execute impl for all CircuitOps
impl<B: CircuitBackend, I: CircuitOp> Execute<I> for CircuitExecutor<B> {
    type Effect = I::Effect;
    type Message = I::Message;
    type Error = eyre::Report;

    fn execute(
        &mut self,
        instruction: &I,
        message: Self::Message,
    ) -> Result<Self::Effect, Self::Error> {
        I::apply(instruction, &mut self.backend, message)
    }
}

/// circuit component that selects qubit size from a runtime value and can do all backends
#[component]
#[instructions { circuit }]
pub struct Circuit {
    inner: CircuitKind,
    n_qubits: usize,
}

impl<I: CircuitOp> Execute<I> for Circuit {
    type Message = I::Message;
    type Effect = I::Effect;
    type Error = eyre::Report;

    fn execute(&mut self, instruction: &I, msg: I::Message) -> Result<Self::Effect, eyre::Report> {
        match &mut self.inner {
            CircuitKind::Tableau(b) => b.apply(instruction, msg),
            CircuitKind::PauliSum(b) => b.apply(instruction, msg),
            CircuitKind::LossyPauliSum(b) => b.apply(instruction, msg),
        }
    }
}

/// Tableau-only component with runtime configured size
#[component]
#[instructions { circuit }]
pub struct TableauCircuit {
    inner: TableauBuckets,
}

impl<I: CircuitOp> Execute<I> for TableauCircuit {
    type Message = I::Message;
    type Effect = I::Effect;
    type Error = eyre::Report;

    fn execute(
        &mut self,
        instruction: &I,
        message: Self::Message,
    ) -> Result<Self::Effect, Self::Error> {
        self.inner.apply(instruction, message)
    }
}

impl vihaco::Reset for Circuit {
    fn reset(&mut self) {
        self.inner.reset()
    }
}

impl vihaco::Reset for TableauCircuit {
    fn reset(&mut self) {
        self.inner.reset()
    }
}

impl Circuit {
    /// Build the backend selected by `info.backend`. `seed` makes tableau runs
    /// reproducible; the PauliSum backends are deterministic and ignore it.
    pub fn new(info: &PPVMDeviceInfo, seed: Option<u64>) -> Result<Self> {
        match (info.backend, seed) {
            (BackendKind::Tableau, None) => Self::tableau(info),
            (BackendKind::Tableau, Some(seed)) => Self::tableau_with_seed(info, seed),
            (BackendKind::PauliSum, _) => Self::paulisum(info, &observable_terms(info)?),
            (BackendKind::LossyPauliSum, _) => Self::lossy_paulisum(info, &observable_terms(info)?),
        }
    }

    /// Number of qubits this circuit was built for.
    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    /// Human-readable rendering of the backend state (tableau or Pauli sum).
    pub fn state_string(&self) -> String {
        self.inner.state_string()
    }

    /// Build a Tableau-backed circuit. Tableau init only needs `n_qubits` and
    /// `coefficient_threshold` from `info`; no observable required.
    pub fn tableau(info: &PPVMDeviceInfo) -> Result<Self> {
        let tab = TableauBuckets::new(info.n_qubits, info.coefficient_threshold)?;
        Ok(Self {
            inner: CircuitKind::Tableau(tab),
            n_qubits: info.n_qubits,
        })
    }

    /// Same as [`Circuit.tableau`], but seed the tableau RNG deterministically
    /// so a shot is reproducible.
    pub fn tableau_with_seed(info: &PPVMDeviceInfo, seed: u64) -> Result<Self> {
        let tab = TableauBuckets::new_with_seed(info.n_qubits, info.coefficient_threshold, seed)?;
        Ok(Self {
            inner: CircuitKind::Tableau(tab),
            n_qubits: info.n_qubits,
        })
    }

    /// Build a PauliSum-backed circuit, seeding the state with every term in
    /// `terms`. Each `(word, coef)` is added via `state += (word, coef)`; the
    /// caller is responsible for having parsed/validated the words against
    /// `info.n_qubits` (see `parse_observable_terms` in `composite.rs`).
    pub fn paulisum(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        let psum = PauliSumBuckets::new(info, terms)?;
        Ok(Self {
            inner: CircuitKind::PauliSum(psum),
            n_qubits: info.n_qubits,
        })
    }

    /// Build a LossyPauliSum-backed circuit. Same contract as
    /// [`Circuit::paulisum`].
    pub fn lossy_paulisum(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        let psum = LossyPauliSumBuckets::new(info, terms)?;
        Ok(Self {
            inner: CircuitKind::LossyPauliSum(psum),
            n_qubits: info.n_qubits,
        })
    }
}

impl Default for Circuit {
    fn default() -> Self {
        // Default backend is Tableau with 0 qubits, which always fits the
        // smallest bucket, so construction here is infallible.
        Self::tableau(&PPVMDeviceInfo::default()).expect("0-qubit tableau is always constructible")
    }
}

impl TableauCircuit {
    /// Build a Tableau-backed circuit. Tableau init only needs `n_qubits` and
    /// `coefficient_threshold` from `info`; no observable required.
    pub fn new(info: &PPVMDeviceInfo) -> Result<Self> {
        let tab = TableauBuckets::new(info.n_qubits, info.coefficient_threshold)?;
        Ok(Self { inner: tab })
    }

    /// Same as [`Circuit.tableau`], but seed the tableau RNG deterministically
    /// so a shot is reproducible.
    pub fn new_with_seed(info: &PPVMDeviceInfo, seed: u64) -> Result<Self> {
        let tab = TableauBuckets::new_with_seed(info.n_qubits, info.coefficient_threshold, seed)?;
        Ok(Self { inner: tab })
    }
}

impl Default for TableauCircuit {
    fn default() -> Self {
        Self::new(&PPVMDeviceInfo::default()).expect("0-qubit tableau is always constructible")
    }
}

/// The PauliSum backends' initial observable, parsed from `circuit.observable`.
fn observable_terms(info: &PPVMDeviceInfo) -> Result<Vec<(String, f64)>> {
    let observable = info.observable.as_deref().ok_or_else(|| {
        eyre!(
            "the {:?} backend requires `circuit.observable` to be set",
            info.backend
        )
    })?;
    crate::observable::parse_pauli_sum_terms(observable, info.n_qubits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::component::MAX_QUBITS;
    use crate::measurements::{CircuitOutcomeEffect, MeasurementOutcome};
    use vihaco::PopFromStack;

    fn info(n_qubits: usize) -> PPVMDeviceInfo {
        PPVMDeviceInfo {
            n_qubits,
            ..Default::default()
        }
    }

    /// Measure `addr` and read the outcome out of the returned effect.
    fn measure(circuit: &mut Circuit, addr: u64) -> MeasurementOutcome {
        let effects = circuit
            .execute(&circuit::Measure, PopFromStack((addr,)))
            .unwrap();
        match effects.into_iter().next() {
            Some(CircuitOutcomeEffect::Measurement(m)) => m.measurement_results[0],
            other => panic!("expected a measurement effect, got {other:?}"),
        }
    }

    #[test]
    fn tableau_backend_x_flips_qubit() {
        let mut circuit = Circuit::tableau(&info(1)).unwrap();
        circuit.execute(&circuit::X, PopFromStack((0,))).unwrap();
        assert_eq!(measure(&mut circuit, 0), MeasurementOutcome::One);
    }

    #[test]
    fn tableau_backend_cnot_propagates_flip() {
        // X(0); CNOT(0, 1) leaves both qubits in |1>.
        let mut circuit = Circuit::tableau(&info(2)).unwrap();
        circuit.execute(&circuit::X, PopFromStack((0,))).unwrap();
        circuit
            .execute(&circuit::CNOT, PopFromStack((0, 1)))
            .unwrap();
        assert_eq!(measure(&mut circuit, 0), MeasurementOutcome::One);
        assert_eq!(measure(&mut circuit, 1), MeasurementOutcome::One);
    }

    #[test]
    fn tableau_t_gate_doubles_coefficients() {
        // T on |+> is non-Clifford, so the tableau splits into two terms; the
        // following CNOT is Clifford and keeps that count.
        let mut circuit = Circuit::tableau(&info(2)).unwrap();
        circuit.execute(&circuit::H, PopFromStack((0,))).unwrap();
        circuit.execute(&circuit::T, PopFromStack((0,))).unwrap();
        circuit
            .execute(&circuit::CNOT, PopFromStack((0, 1)))
            .unwrap();
        let CircuitKind::Tableau(TableauBuckets::Bits64(executor)) = &circuit.inner else {
            panic!("a 2-qubit tableau should use the 64-bit bucket");
        };
        assert_eq!(executor.backend.coefficients.len(), 2);
    }

    // The widest bucket is 2048 qubits; beyond it constructors must error
    // rather than panic, since a panic would tear down the TUI.

    #[test]
    fn tableau_rejects_more_than_2048_qubits() {
        assert!(Circuit::tableau(&info(MAX_QUBITS + 1)).is_err());
    }

    #[test]
    fn tableau_with_seed_rejects_more_than_2048_qubits() {
        assert!(Circuit::tableau_with_seed(&info(MAX_QUBITS + 1), 0).is_err());
    }

    #[test]
    fn paulisum_rejects_more_than_2048_qubits() {
        assert!(Circuit::paulisum(&info(MAX_QUBITS + 1), &[]).is_err());
    }

    #[test]
    fn lossy_paulisum_rejects_more_than_2048_qubits() {
        assert!(Circuit::lossy_paulisum(&info(MAX_QUBITS + 1), &[]).is_err());
    }

    #[test]
    fn constructs_at_the_2048_qubit_boundary() {
        assert!(Circuit::tableau(&info(MAX_QUBITS)).is_ok());
    }
}
