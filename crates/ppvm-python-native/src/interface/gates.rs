// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

macro_rules! create_interface_gates {
    ($name: ident, $type: ident, $loss: tt) => {
        #[pymethods]
        impl $name {
            // clifford
            #[pyo3(signature = (targets, truncate = true))]
            pub fn x(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.x_many(targets.as_slice());
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn y(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.y_many(targets.as_slice());
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn z(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.z_many(targets.as_slice());
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn h(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.h_many(targets.as_slice());
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn s(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.s_many(targets.as_slice());
                if truncate {
                    self.inner.truncate();
                }
            }

            // two-qubit clifford (+ stim aliases)
            #[pyo3(signature = (targets, truncate = true))]
            pub fn cnot(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.cnot_many(&pairs);
                if truncate {
                    self.inner.truncate();
                }
                Ok(())
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn cx(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                self.cnot(targets, truncate)
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn zcx(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                self.cnot(targets, truncate)
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn cz(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.cz_many(&pairs);
                if truncate {
                    self.inner.truncate();
                }
                Ok(())
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn zcz(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                self.cz(targets, truncate)
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn cy(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.cy_many(&pairs);
                if truncate {
                    self.inner.truncate();
                }
                Ok(())
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn zcy(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                self.cy(targets, truncate)
            }

            // clifford extensions
            #[pyo3(signature = (targets, truncate = true))]
            pub fn s_dag(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.s_dag_many(targets.as_slice());
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn sqrt_x(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.sqrt_x_many(targets.as_slice());
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn sqrt_y(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.sqrt_y_many(targets.as_slice());
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn sqrt_x_dag(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.sqrt_x_dag_many(targets.as_slice());
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn sqrt_y_dag(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.sqrt_y_dag_many(targets.as_slice());
                if truncate {
                    self.inner.truncate();
                }
            }
        }
    };
}
