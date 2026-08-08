// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

// Adapted from https://pyo3.rs/v0.27.1/class.html#no-generic-parameters
macro_rules! create_interface {
    ($name: ident, $type: ident, $loss: tt) => {
        create_interface_state!($name, $type, $loss);
        create_interface_gates!($name, $type, $loss);
        create_interface_noise!($name, $type, $loss);
        create_interface_rotations!($name, $type, $loss);
        create_interface_python!($name, $type, $loss);
        create_interface_loss_methods!($name, $type, $loss);
    };
}

macro_rules! create_interface_range {
    ($name: ident, false, $( $n: expr),* ) => {
        paste! {
        $(
            type [<$name$n>] = crate::backend::OrdinaryPauliSum<{(2 as usize).pow($n)}>;
            create_interface!([<PauliSum$name$n>], [<$name$n>], false);
        )*
    }
    };

    ($name: ident, true, $( $n: expr),* ) => {
        paste! {
        $(
            type [<Loss$name$n>] = crate::backend::LossyPauliSum<{(2 as usize).pow($n)}>;
            create_interface!([<PauliSumLoss$name$n>], [<Loss$name$n>], true);
        )*
    }
    };
}
