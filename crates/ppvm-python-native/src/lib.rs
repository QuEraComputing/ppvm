// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

pub mod interface;
pub mod interface_tableau;
pub mod interface_tableau_sum;
pub mod stim_program;

pub(crate) fn flat_pairs(targets: &[usize]) -> PyResult<Vec<(usize, usize)>> {
    if !targets.len().is_multiple_of(2) {
        return Err(PyValueError::new_err(
            "two-qubit operations require an even number of targets",
        ));
    }
    Ok(targets
        .chunks_exact(2)
        .map(|pair| (pair[0], pair[1]))
        .collect())
}

// Imported by Python only as the private `ppvm._core` submodule; maturin's
// `module-name = "ppvm._core"` maps this `PyInit__core` symbol into the wheel.
#[pymodule]
pub mod _core {
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

    // Generalized Tableau Sum
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum1;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum2;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum3;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum4;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum5;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum6;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum7;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum8;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum9;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum10;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum11;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum12;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum13;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum14;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum15;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum16;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum17;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum18;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum19;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum20;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum21;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum22;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum23;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum24;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum25;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum26;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum27;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum28;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum29;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum30;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum31;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::GeneralizedTableauSum32;

    // Tableau Sum Sampler
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler1;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler2;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler3;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler4;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler5;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler6;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler7;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler8;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler9;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler10;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler11;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler12;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler13;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler14;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler15;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler16;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler17;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler18;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler19;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler20;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler21;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler22;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler23;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler24;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler25;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler26;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler27;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler28;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler29;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler30;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler31;
    #[pymodule_export]
    pub use crate::interface_tableau_sum::TableauSumSampler32;

    // Stim
    #[pymodule_export]
    pub use crate::stim_program::PyStimProgram;
}
