// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Declarations in a PPVM program's `.header(...)` section, collected into the
//! program metadata and a [`PPVMDeviceInfo`].

use crate::device_info::{BackendKind, PPVMDeviceInfo};
use vihaco::module::ProgramMetadata;
use vihaco::{Parse, Type, Value};
use vihaco_parser::Ident;

pub type Metadata = ProgramMetadata<Value, Type>;

#[derive(Debug, Clone, PartialEq, Parse)]
#[syntax_class(value)]
pub enum BackendKindSyntax {
    #[pattern = "`tableau`"]
    Tableau,
    #[pattern = "`paulisum`"]
    PauliSum,
    #[pattern = "`lossy_paulisum`"]
    LossyPauliSum,
}

impl From<BackendKindSyntax> for BackendKind {
    fn from(value: BackendKindSyntax) -> Self {
        match value {
            BackendKindSyntax::Tableau => Self::Tableau,
            BackendKindSyntax::PauliSum => Self::PauliSum,
            BackendKindSyntax::LossyPauliSum => Self::LossyPauliSum,
        }
    }
}

/// One `;`-terminated declaration in a PPVM header section.
#[derive(Debug, Parse)]
#[syntax_class(metadata)]
pub enum HeaderDeclaration {
    #[pattern = "$0"]
    Metadata(Metadata),
    #[pattern = "circuit.n_qubits $0"]
    NumQubits(usize),
    #[pattern = "circuit.coefficient_threshold $0"]
    CoefficientThreshold(f64),
    #[pattern = "circuit.backend $0"]
    Backend(BackendKindSyntax),
    #[pattern = "circuit.observable $0"]
    Observable(Ident),
    #[pattern = "circuit.max_pauli_weight $0"]
    MaxPauliWeight(usize),
}

/// Collect header declarations into the program metadata and device info.
///
/// Exactly one metadata block is required. A repeated `circuit.*` declaration
/// overwrites the earlier one, as in the previous header syntax.
pub fn collect_header(
    declarations: impl IntoIterator<Item = HeaderDeclaration>,
) -> eyre::Result<(Metadata, PPVMDeviceInfo)> {
    let mut metadata = None;
    let mut info = PPVMDeviceInfo::default();
    for declaration in declarations {
        match declaration {
            HeaderDeclaration::Metadata(m) => {
                eyre::ensure!(metadata.is_none(), "duplicate metadata declaration");
                metadata = Some(m);
            }
            HeaderDeclaration::NumQubits(n) => info.n_qubits = n,
            HeaderDeclaration::CoefficientThreshold(t) => info.coefficient_threshold = t,
            HeaderDeclaration::Backend(b) => info.backend = b.into(),
            HeaderDeclaration::Observable(s) => info.observable = Some(s.0),
            HeaderDeclaration::MaxPauliWeight(w) => info.max_pauli_weight = Some(w),
        }
    }
    let metadata = metadata.ok_or_else(|| eyre::eyre!("missing program metadata"))?;
    Ok((metadata, info))
}
