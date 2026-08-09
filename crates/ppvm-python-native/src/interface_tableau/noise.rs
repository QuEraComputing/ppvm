// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

macro_rules! create_tableau_noise {
    ($name: ident, $type: ident, $storage:expr, $indexType: ident) => {
        #[pymethods]
        impl $name {
            // noise
            pub fn x_error(&mut self, targets: Vec<usize>, p: f64) {
                draw!(self.inner.x_error_many(targets.as_slice(), p));
            }

            pub fn y_error(&mut self, targets: Vec<usize>, p: f64) {
                draw!(self.inner.y_error_many(targets.as_slice(), p));
            }

            pub fn z_error(&mut self, targets: Vec<usize>, p: f64) {
                draw!(self.inner.z_error_many(targets.as_slice(), p));
            }

            pub fn pauli_error(&mut self, targets: Vec<usize>, p: [f64; 3]) {
                draw!(self.inner.pauli_error_many(targets.as_slice(), p));
            }

            pub fn depolarize1(&mut self, targets: Vec<usize>, p: f64) {
                draw!(self.inner.depolarize1_many(targets.as_slice(), p));
            }

            pub fn depolarize2(&mut self, targets: Vec<usize>, p: f64) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                draw!(self.inner.depolarize2_many(&pairs, p));
                Ok(())
            }

            pub fn two_qubit_pauli_error(
                &mut self,
                targets: Vec<usize>,
                p: [f64; 15],
            ) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                draw!(self.inner.two_qubit_pauli_error_many(&pairs, p));
                Ok(())
            }

            pub fn loss_channel(&mut self, addr0: usize, p: f64) {
                draw!(self.inner.loss_channel(addr0, p));
            }

            pub fn correlated_loss_channel(&mut self, addr0: usize, addr1: usize, p: [f64; 3]) {
                draw!(self.inner.correlated_loss_channel(addr0, addr1, p));
            }

            pub fn reset_loss_channel(&mut self, addr0: usize) {
                self.inner.reset_loss_channel(addr0);
            }

            pub fn asymmetric_loss_channel(&mut self, addr0: usize, p0: f64, p1: f64) {
                draw!(self.inner.asymmetric_loss_channel(addr0, p0, p1));
            }

            pub fn reset(&mut self, targets: Vec<usize>) {
                draw!(self.inner.reset_many(targets.as_slice()));
            }

            pub fn reset_x(&mut self, targets: Vec<usize>) {
                draw!(self.inner.reset_x_many(targets.as_slice()));
            }

            pub fn reset_y(&mut self, targets: Vec<usize>) {
                draw!(self.inner.reset_y_many(targets.as_slice()));
            }

            pub fn reset_z(&mut self, targets: Vec<usize>) {
                draw!(self.inner.reset_z_many(targets.as_slice()));
            }

            pub fn is_lost(&self, addr0: usize) -> bool {
                self.inner.is_lost[addr0]
            }

            pub fn loss_values(&self) -> Vec<bool> {
                self.inner.is_lost.clone()
            }
        }
    };
}
