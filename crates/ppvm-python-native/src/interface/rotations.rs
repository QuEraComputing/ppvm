// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

macro_rules! create_interface_rotations {
    ($name: ident, $type: ident, $loss: tt) => {
        #[pymethods]
        impl $name {
            // rot1
            #[pyo3(signature = (targets, theta, truncate = true))]
            pub fn rx(&mut self, targets: Vec<usize>, theta: f64, truncate: bool) {
                self.inner.rx_many(targets.as_slice(), theta);
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, theta, truncate = true))]
            pub fn ry(&mut self, targets: Vec<usize>, theta: f64, truncate: bool) {
                self.inner.ry_many(targets.as_slice(), theta);
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (targets, theta, truncate = true))]
            pub fn rz(&mut self, targets: Vec<usize>, theta: f64, truncate: bool) {
                self.inner.rz_many(targets.as_slice(), theta);
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (addr0, axis_angle, theta, truncate = true))]
            pub fn r(&mut self, addr0: usize, axis_angle: f64, theta: f64, truncate: bool) {
                self.inner.r(addr0, axis_angle, theta);
                if truncate {
                    self.inner.truncate();
                }
            }

            // rot2
            #[pyo3(signature = (targets, theta, truncate = true))]
            pub fn rxx(&mut self, targets: Vec<usize>, theta: f64, truncate: bool) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.rxx_many(&pairs, theta);
                if truncate {
                    self.inner.truncate();
                }
                Ok(())
            }

            #[pyo3(signature = (targets, theta, truncate = true))]
            pub fn ryy(&mut self, targets: Vec<usize>, theta: f64, truncate: bool) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.ryy_many(&pairs, theta);
                if truncate {
                    self.inner.truncate();
                }
                Ok(())
            }

            #[pyo3(signature = (targets, theta, truncate = true))]
            pub fn rzz(&mut self, targets: Vec<usize>, theta: f64, truncate: bool) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.rzz_many(&pairs, theta);
                if truncate {
                    self.inner.truncate();
                }
                Ok(())
            }
        }
    };
}
