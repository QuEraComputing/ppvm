// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Comparisons include width and explicit zero terms, but ignore configuration.

use approx::{AbsDiffEq, RelativeEq};
use ppvm_traits_2::Support;

use crate::Sum;

impl<S: Support, P> Sum<S, P> {
    fn compare_coefficients(
        &self,
        other: &Self,
        equal: impl Fn(&S::Coeff, &S::Coeff) -> bool,
    ) -> bool {
        self.n_sites() == other.n_sites()
            && self.len() == other.len()
            && self.iter().all(|(key, coeff)| {
                other
                    .get(&key)
                    .is_some_and(|other_coeff| equal(&coeff, &other_coeff))
            })
    }
}

impl<S: Support, P> PartialEq for Sum<S, P> {
    fn eq(&self, other: &Self) -> bool {
        self.compare_coefficients(other, |a, b| a == b)
    }
}

impl<S: Support, P> Eq for Sum<S, P> where S::Coeff: Eq {}

impl<S: Support, P> AbsDiffEq for Sum<S, P>
where
    S::Coeff: AbsDiffEq,
    <S::Coeff as AbsDiffEq>::Epsilon: Copy,
{
    type Epsilon = <S::Coeff as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        S::Coeff::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.compare_coefficients(other, |a, b| a.abs_diff_eq(b, epsilon))
    }
}

impl<S: Support, P> RelativeEq for Sum<S, P>
where
    S::Coeff: RelativeEq,
    <S::Coeff as AbsDiffEq>::Epsilon: Copy,
{
    fn default_max_relative() -> Self::Epsilon {
        S::Coeff::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        self.compare_coefficients(other, |a, b| a.relative_eq(b, epsilon, max_relative))
    }
}
