// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

macro_rules! create_sum_interface {
    ($tab_name: ident, $sampler_name: ident, $type: ident, $storage:expr, $indexType: ident) => {
        create_sum_state!($tab_name, $sampler_name, $type, $storage, $indexType);
        create_sum_gates!($tab_name, $sampler_name, $type, $storage, $indexType);
        create_sum_noise!($tab_name, $sampler_name, $type, $storage, $indexType);
        create_sum_sampler!($tab_name, $sampler_name, $type, $storage, $indexType);
    };
}

macro_rules! create_sum_interface_range {
    ($indexType: ident, $( $n: expr),* ) => {
        paste! {
        $(
            type [<SumConfig$n>] = crate::backend::GeneralizedTableauSum<$n, $indexType>;
            create_sum_interface!(
                [<GeneralizedTableauSum$n>],
                [<TableauSumSampler$n>],
                [<SumConfig$n>],
                $n,
                $indexType
            );
        )*
        }
    };
}
