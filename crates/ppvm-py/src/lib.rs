use pyo3::prelude::*;

pub mod interface;

#[pymodule]
pub mod ppvm_py {
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash10;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash100;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash500;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash1000;
}
