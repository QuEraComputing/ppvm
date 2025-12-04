use pyo3::prelude::*;

pub mod interface;

#[pymodule]
pub mod ppvm_python {
    // NOTE: it's not possible to use #[pymodule_export] inside a macro_rules!
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash0;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash1;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash2;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash3;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash4;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash5;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash6;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash7;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash8;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash9;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash10;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash11;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash12;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash13;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash14;
    #[pymodule_export]
    pub use crate::interface::PauliSumIndexMapFxHash15;
}
