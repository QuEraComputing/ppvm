// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

macro_rules! create_interface_loss_methods {
    ($name: ident, $type: ident, false) => {};
    ($name: ident, $type: ident, true) => {
        #[pymethods]
        impl $name {
            #[pyo3(signature = (addr0, p, truncate = true))]
            pub fn loss_channel(&mut self, addr0: usize, p: f64, truncate: bool) {
                draw!(self.inner.loss_channel(addr0, p));
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (addr0, addr1, p, truncate = true))]
            pub fn correlated_loss_channel(
                &mut self,
                addr0: usize,
                addr1: usize,
                p: [f64; 3],
                truncate: bool,
            ) {
                draw!(self.inner.correlated_loss_channel(addr0, addr1, p));
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (addr0, truncate = true))]
            pub fn reset_loss_channel(&mut self, addr0: usize, truncate: bool) {
                self.inner.reset_loss_channel(addr0);
                if truncate {
                    self.inner.truncate();
                }
            }
        }
    };
}
