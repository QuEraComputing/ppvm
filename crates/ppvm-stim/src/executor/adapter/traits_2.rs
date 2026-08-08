// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_tableau_2::prelude::{
    Bitstring, Clifford, CliffordBatch, CliffordExtensions, CliffordExtensionsBatch,
    CorrelatedLossChannel, Depolarizing, Depolarizing2, GeneralizedTableau, LossChannel, Measure,
    PauliError, Reset, RotationOne, RowStorage, TGate, TwoQubitPauliError, U3Gate,
};

use super::StimTableau;

macro_rules! unary {
    ($name:ident, $trait:ident) => {
        fn $name(&mut self, q: usize) {
            $trait::$name(self, q);
        }
    };
}
macro_rules! binary {
    ($name:ident, $trait:ident) => {
        fn $name(&mut self, a: usize, b: usize) {
            $trait::$name(self, a, b);
        }
    };
}
macro_rules! batch {
    ($name:ident, $trait:ident, $ty:ty) => {
        fn $name(&mut self, q: &[$ty]) {
            $trait::$name(self, q);
        }
    };
}

impl<A, I, H> StimTableau for GeneralizedTableau<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
{
    type Config = A;
    type Index = I;
    type Coefficients = H;

    unary!(reset, Reset);
    unary!(x, Clifford);
    unary!(y, Clifford);
    unary!(z, Clifford);
    unary!(h, Clifford);
    unary!(s, Clifford);
    unary!(s_dag, CliffordExtensions);
    unary!(sqrt_x, CliffordExtensions);
    unary!(sqrt_x_dag, CliffordExtensions);
    unary!(sqrt_y, CliffordExtensions);
    unary!(sqrt_y_dag, CliffordExtensions);
    binary!(cnot, Clifford);
    binary!(cy, CliffordExtensions);
    binary!(cz, Clifford);
    batch!(x_many, CliffordBatch, usize);
    batch!(y_many, CliffordBatch, usize);
    batch!(z_many, CliffordBatch, usize);
    batch!(h_many, CliffordBatch, usize);
    batch!(s_many, CliffordBatch, usize);
    batch!(s_dag_many, CliffordExtensionsBatch, usize);
    batch!(sqrt_x_many, CliffordExtensionsBatch, usize);
    batch!(sqrt_x_dag_many, CliffordExtensionsBatch, usize);
    batch!(sqrt_y_many, CliffordExtensionsBatch, usize);
    batch!(sqrt_y_dag_many, CliffordExtensionsBatch, usize);
    batch!(cnot_many, CliffordBatch, (usize, usize));
    batch!(cy_many, CliffordExtensionsBatch, (usize, usize));
    batch!(cz_many, CliffordBatch, (usize, usize));
    unary!(t, TGate);
    unary!(t_dag, TGate);

    fn rx(&mut self, q: usize, theta: f64) {
        RotationOne::rx(self, q, theta);
    }
    fn ry(&mut self, q: usize, theta: f64) {
        RotationOne::ry(self, q, theta);
    }
    fn rz(&mut self, q: usize, theta: f64) {
        RotationOne::rz(self, q, theta);
    }
    fn u3(&mut self, q: usize, theta: f64, phi: f64, lambda: f64) {
        U3Gate::u3(self, q, theta, phi, lambda);
    }
    fn depolarize1(&mut self, q: usize, p: f64) {
        Depolarizing::depolarize1(self, q, p);
    }
    fn depolarize2(&mut self, a: usize, b: usize, p: f64) {
        Depolarizing2::depolarize2(self, a, b, p);
    }
    fn pauli_error(&mut self, q: usize, p: [f64; 3]) {
        PauliError::pauli_error(self, q, p);
    }
    fn two_qubit_pauli_error(&mut self, a: usize, b: usize, p: [f64; 15]) {
        TwoQubitPauliError::two_qubit_pauli_error(self, a, b, p);
    }
    fn loss_channel(&mut self, q: usize, p: f64) {
        LossChannel::loss_channel(self, q, p);
    }
    fn correlated_loss_channel(&mut self, a: usize, b: usize, p: [f64; 3]) {
        CorrelatedLossChannel::correlated_loss_channel(self, a, b, p);
    }
    fn measure(&mut self, q: usize) -> Option<bool> {
        Measure::measure(self, q)
    }
    fn measure_many(&mut self, q: &[usize]) -> Vec<Option<bool>> {
        Measure::measure_many(self, q)
    }
    fn measure_noisy(&mut self, q: usize, p: f64) -> Option<bool> {
        GeneralizedTableau::measure_noisy(self, q, p)
    }
    fn flip_with_prob(&mut self, bit: bool, p: f64) -> bool {
        GeneralizedTableau::flip_with_prob(self, bit, p)
    }
    fn measurement_record(&self) -> &[Option<bool>] {
        self.current_measurement_record()
    }
    fn append_measurement_record(&mut self, result: Option<bool>) {
        GeneralizedTableau::append_measurement_record(self, result);
    }
    fn overwrite_last_measurement_record(&mut self, result: Option<bool>) {
        GeneralizedTableau::overwrite_last_measurement_record(self, result);
    }
}
