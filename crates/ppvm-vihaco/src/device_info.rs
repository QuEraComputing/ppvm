// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Device configuration for the circuit component. Split into its own module
//! so [`crate::component`] can depend on [`PPVMDeviceInfo`] without pulling in
//! the composite machine (which is added in a later PR).

pub const PPVM_MAGIC: u32 = 0x5050564D;

/// Which execution backend the circuit runs on. Selected via the
/// `device circuit.backend` header; defaults to `Tableau` so existing
/// programs that don't declare a backend keep working.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, vihaco_parser::Parse)]
pub enum BackendKind {
    #[default]
    Tableau,
    PauliSum,
    #[token = "lossy_paulisum"]
    LossyPauliSum,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PPVMDeviceInfo {
    pub magic: u32,
    pub n_qubits: usize,
    pub coefficient_threshold: f64,
    pub backend: BackendKind,
    pub observable: Option<String>,
    pub max_pauli_weight: Option<usize>,
}

impl Default for PPVMDeviceInfo {
    fn default() -> Self {
        Self {
            magic: PPVM_MAGIC,
            n_qubits: 0,
            coefficient_threshold: 1e-10,
            backend: BackendKind::default(),
            observable: None,
            max_pauli_weight: None,
        }
    }
}
