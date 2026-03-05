use pyo3::prelude::*;

pub mod interface;
pub mod interface_tableau;

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
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau33;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau34;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau35;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau36;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau37;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau38;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau39;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau40;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau41;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau42;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau43;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau44;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau45;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau46;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau47;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau48;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau49;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau50;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau51;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau52;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau53;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau54;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau55;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau56;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau57;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau58;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau59;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau60;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau61;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau62;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau63;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau64;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau65;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau66;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau67;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau68;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau69;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau70;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau71;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau72;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau73;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau74;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau75;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau76;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau77;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau78;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau79;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau80;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau81;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau82;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau83;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau84;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau85;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau86;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau87;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau88;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau89;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau90;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau91;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau92;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau93;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau94;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau95;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau96;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau97;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau98;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau99;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau100;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau101;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau102;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau103;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau104;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau105;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau106;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau107;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau108;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau109;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau110;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau111;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau112;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau113;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau114;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau115;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau116;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau117;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau118;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau119;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau120;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau121;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau122;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau123;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau124;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau125;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau126;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau127;
    #[pymodule_export]
    pub use crate::interface_tableau::GeneralizedTableau128;
}
