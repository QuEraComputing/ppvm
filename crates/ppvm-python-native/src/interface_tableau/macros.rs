// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

macro_rules! create_interface {
    ($name: ident, $type: ident, $storage:expr, $indexType: ident) => {
        create_tableau_state!($name, $type, $storage, $indexType);
        create_tableau_gates!($name, $type, $storage, $indexType);
        create_tableau_noise!($name, $type, $storage, $indexType);
        create_tableau_stim!($name, $type, $storage, $indexType);
    };
}

macro_rules! create_interface_range {
    ($name: ident, $indexType: ident, $( $n: expr),* ) => {
        paste! {
        $(
            type [<$name$n>] = crate::backend::GeneralizedTableau<$n, $indexType>;
            create_interface!([<GeneralizedTableau$n>], [<$name$n>], $n, $indexType);
        )*
    }
    };
}
