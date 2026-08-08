// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

macro_rules! create_tableau_gates {
    ($name: ident, $type: ident, $storage:expr, $indexType: ident) => {
        #[pymethods]
        impl $name {
            // clifford
            pub fn x(&mut self, targets: Vec<usize>) {
                self.inner.x_many(targets.as_slice());
            }

            pub fn y(&mut self, targets: Vec<usize>) {
                self.inner.y_many(targets.as_slice());
            }

            pub fn z(&mut self, targets: Vec<usize>) {
                self.inner.z_many(targets.as_slice());
            }

            pub fn h(&mut self, targets: Vec<usize>) {
                self.inner.h_many(targets.as_slice());
            }

            pub fn s(&mut self, targets: Vec<usize>) {
                self.inner.s_many(targets.as_slice());
            }

            pub fn s_dag(&mut self, targets: Vec<usize>) {
                self.inner.s_dag_many(targets.as_slice());
            }

            // clifford extensions
            pub fn sqrt_x(&mut self, targets: Vec<usize>) {
                self.inner.sqrt_x_many(targets.as_slice());
            }

            pub fn sqrt_x_dag(&mut self, targets: Vec<usize>) {
                self.inner.sqrt_x_dag_many(targets.as_slice());
            }

            pub fn sqrt_y(&mut self, targets: Vec<usize>) {
                self.inner.sqrt_y_many(targets.as_slice());
            }

            pub fn sqrt_y_dag(&mut self, targets: Vec<usize>) {
                self.inner.sqrt_y_dag_many(targets.as_slice());
            }

            pub fn t(&mut self, targets: Vec<usize>) {
                self.inner.t_many(targets.as_slice());
            }

            pub fn t_dag(&mut self, targets: Vec<usize>) {
                self.inner.t_dag_many(targets.as_slice());
            }

            // two-qubit clifford (+ stim aliases)
            pub fn cnot(&mut self, targets: Vec<usize>) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.cnot_many(&pairs);
                Ok(())
            }

            pub fn cx(&mut self, targets: Vec<usize>) -> PyResult<()> {
                self.cnot(targets)
            }

            pub fn zcx(&mut self, targets: Vec<usize>) -> PyResult<()> {
                self.cnot(targets)
            }

            pub fn cy(&mut self, targets: Vec<usize>) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.cy_many(&pairs);
                Ok(())
            }

            pub fn zcy(&mut self, targets: Vec<usize>) -> PyResult<()> {
                self.cy(targets)
            }

            pub fn cz(&mut self, targets: Vec<usize>) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.cz_many(&pairs);
                Ok(())
            }

            pub fn zcz(&mut self, targets: Vec<usize>) -> PyResult<()> {
                self.cz(targets)
            }

            pub fn cz_block(&mut self, control_base: usize, target_base: usize, count: usize) {
                self.inner.cz_block(control_base, target_base, count);
            }

            // rot1
            pub fn rx(&mut self, targets: Vec<usize>, theta: f64) {
                self.inner.rx_many(targets.as_slice(), theta);
            }

            pub fn ry(&mut self, targets: Vec<usize>, theta: f64) {
                self.inner.ry_many(targets.as_slice(), theta);
            }

            pub fn rz(&mut self, targets: Vec<usize>, theta: f64) {
                self.inner.rz_many(targets.as_slice(), theta);
            }

            pub fn u3(&mut self, addr0: usize, theta: f64, phi: f64, lam: f64) {
                self.inner.u3(addr0, theta, phi, lam);
            }

            pub fn r(&mut self, addr0: usize, axis_angle: f64, theta: f64) {
                self.inner.r(addr0, axis_angle, theta);
            }

            // rot2
            pub fn rxx(&mut self, targets: Vec<usize>, theta: f64) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.rxx_many(&pairs, theta);
                Ok(())
            }

            pub fn ryy(&mut self, targets: Vec<usize>, theta: f64) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.ryy_many(&pairs, theta);
                Ok(())
            }

            pub fn rzz(&mut self, targets: Vec<usize>, theta: f64) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.rzz_many(&pairs, theta);
                Ok(())
            }
        }
    };
}
