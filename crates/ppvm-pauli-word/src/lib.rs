// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Packed Pauli-word types: `PauliWord`, `PhasedPauliWord`, `LossyPauliWord`,
//! and `PauliPattern`. Built on the `ppvm-traits` foundation.
pub mod loss;
pub mod pattern;
pub mod phase;
pub mod word;

pub mod prelude {
    pub use crate::loss::LossyPauliWord;
    pub use crate::pattern::PauliPattern;
    pub use crate::phase::PhasedPauliWord;
    pub use crate::word::PauliWord;
    pub use ppvm_traits::prelude::*;
}
