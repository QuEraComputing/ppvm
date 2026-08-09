// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

macro_rules! create_interface_noise {
    ($name: ident, $type: ident, $loss: tt) => {
        #[pymethods]
        impl $name {
            // noise
            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn x_error(&mut self, targets: Vec<usize>, p: f64, truncate: bool) {
                draw!(self.inner.x_error_many(targets.as_slice(), p));
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn y_error(&mut self, targets: Vec<usize>, p: f64, truncate: bool) {
                draw!(self.inner.y_error_many(targets.as_slice(), p));
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn z_error(&mut self, targets: Vec<usize>, p: f64, truncate: bool) {
                draw!(self.inner.z_error_many(targets.as_slice(), p));
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn pauli_error(&mut self, targets: Vec<usize>, p: [f64; 3], truncate: bool) {
                draw!(self.inner.pauli_error_many(targets.as_slice(), p));
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn two_qubit_pauli_error(
                &mut self,
                targets: Vec<usize>,
                p: [f64; 15],
                truncate: bool,
            ) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                draw!(self.inner.two_qubit_pauli_error_many(&pairs, p));
                if truncate {
                    self.inner.truncate();
                }
                Ok(())
            }

            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn depolarize1(&mut self, targets: Vec<usize>, p: f64, truncate: bool) {
                draw!(self.inner.depolarize1_many(targets.as_slice(), p));
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn depolarize2(
                &mut self,
                targets: Vec<usize>,
                p: f64,
                truncate: bool,
            ) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                draw!(self.inner.depolarize2_many(&pairs, p));
                if truncate {
                    self.inner.truncate();
                }
                Ok(())
            }

            #[pyo3(signature = (addr0, gamma, truncate = true))]
            pub fn amplitude_damping(&mut self, addr0: usize, gamma: f64, truncate: bool) {
                self.inner.amplitude_damping(addr0, gamma);
                if truncate {
                    self.inner.truncate();
                }
            }
        }
    };
}
