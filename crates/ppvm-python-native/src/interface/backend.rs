// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

#[cfg(feature = "legacy")]
macro_rules! construct_pauli_sum {
    ($type:ty, $strategy:expr, $n:expr, $preserve:expr, $terms:expr, $coeffs:expr) => {{
        let mut ps = <$type>::builder()
            .n_qubits($n)
            .strategy($strategy)
            .capacity($n)
            .preserve_strings($preserve.into_iter().map(Into::into).collect())
            .build();
        for (term, coefficient) in $terms.iter().zip($coeffs.iter()) {
            ps += (term.to_owned(), *coefficient);
        }
        ps
    }};
}

#[cfg(feature = "traits-2")]
macro_rules! construct_pauli_sum {
    ($type:ty, $strategy:expr, $n:expr, $preserve:expr, $terms:expr, $coeffs:expr) => {{
        let mut ps = <$type>::with_capacity($n, $strategy, $n)
            .preserving($preserve.into_iter().map(Into::into));
        for (term, coefficient) in $terms.iter().zip($coeffs.iter()) {
            ps += (term.to_owned(), *coefficient);
        }
        ps
    }};
}

#[cfg(feature = "legacy")]
macro_rules! sum_iter {
    ($sum:expr) => {
        $sum.data().iter().map(|(key, value)| (key.clone(), *value))
    };
}

#[cfg(feature = "traits-2")]
macro_rules! sum_iter {
    ($sum:expr) => {
        $sum.iter()
    };
}
#[cfg(feature = "legacy")]
macro_rules! create_strategy {
    (false, $min_abs_coeff:ident, $max_pauli_weight:ident, $_max_loss_weight:ident) => {
        CombinedStrategy(
            CoefficientThreshold($min_abs_coeff),
            MaxPauliWeight($max_pauli_weight),
        )
    };
    (true, $min_abs_coeff:ident, $max_pauli_weight:ident, $max_loss_weight:ident) => {
        CombinedStrategy(
            CombinedStrategy(
                CoefficientThreshold($min_abs_coeff),
                MaxPauliWeight($max_pauli_weight),
            ),
            MaxLossWeight($max_loss_weight),
        )
    };
}

#[cfg(feature = "traits-2")]
macro_rules! create_strategy {
    (false, $min_abs_coeff:ident, $max_pauli_weight:ident, $_max_loss_weight:ident) => {
        CombinedPolicy(
            CoefficientThreshold {
                threshold: $min_abs_coeff,
            },
            MaxPauliWeight($max_pauli_weight),
        )
    };
    (true, $min_abs_coeff:ident, $max_pauli_weight:ident, $max_loss_weight:ident) => {
        CombinedPolicy(
            CombinedPolicy(
                CoefficientThreshold {
                    threshold: $min_abs_coeff,
                },
                MaxPauliWeight($max_pauli_weight),
            ),
            MaxLossWeight($max_loss_weight),
        )
    };
}
