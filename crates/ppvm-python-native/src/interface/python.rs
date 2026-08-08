// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

macro_rules! create_interface_python {
    ($name: ident, $type: ident, $loss: tt) => {
        #[pymethods]
        impl $name {
            // some python niceties

            fn __copy__(&self) -> Self {
                Self {
                    inner: self.inner.clone(),
                }
            }

            fn __richcmp__(
                &self,
                other: PyRef<$name>,
                op: pyo3::basic::CompareOp,
            ) -> PyResult<bool> {
                match op {
                    pyo3::basic::CompareOp::Eq => Ok(self.inner == other.inner),
                    pyo3::basic::CompareOp::Ne => Ok(self.inner != other.inner),
                    _ => Err(pyo3::exceptions::PyNotImplementedError::new_err(
                        "Only equality and inequality comparisons are supported for PauliSum.",
                    )),
                }
            }

            fn __deepcopy__(&self, _memo: &Bound<'_, PyAny>) -> Self {
                Self {
                    inner: self.inner.clone(),
                }
            }

            fn __len__(&self) -> usize {
                self.inner.len()
            }

            pub fn terms(&self) -> Vec<(String, f64)> {
                sum_iter!(self.inner)
                    .map(|(k, v)| (k.to_string(), v))
                    .collect()
            }

            pub fn weights(&self) -> Vec<(String, usize)> {
                sum_iter!(self.inner)
                    .map(|(k, _v)| (k.to_string(), k.weight()))
                    .collect()
            }

            pub fn current_max_weight(&self) -> usize {
                sum_iter!(self.inner)
                    .map(|(k, _v)| k.weight())
                    .max()
                    .unwrap_or(0)
            }
        }
    };
}
