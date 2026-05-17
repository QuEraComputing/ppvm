use pyo3::prelude::*;

pub mod interface;
pub mod interface_tableau;
pub mod stim_program;

#[pymodule]
pub mod ppvm_python_native {
    // NOTE: it's not possible to use #[pymodule_export] inside a macro_rules!

    // PauliSum
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

    // PauliSum with Loss
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash0;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash1;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash2;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash3;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash4;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash5;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash6;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash7;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash8;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash9;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash10;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash11;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash12;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash13;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash14;
    #[pymodule_export]
    pub use crate::interface::PauliSumLossIndexMapFxHash15;

    // Generalized Tableau
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau1;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau2;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau3;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau4;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau5;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau6;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau7;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau8;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau9;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau10;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau11;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau12;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau13;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau14;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau15;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau16;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau17;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau18;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau19;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau20;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau21;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau22;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau23;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau24;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau25;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau26;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau27;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau28;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau29;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau30;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau31;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau32;

    // Stim
    #[pymodule_export]
    pub use crate::stim_program::PyStimProgram;
}
