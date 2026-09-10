// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::arithmetic::Coefficient;

/// A coefficient ring in which halving (`0.5·x`) is total and exact: the
/// capability the projective computational-basis measurement kernel needs to
/// apply the `(I ± Z)/2` projectors.
pub trait Halvable: Coefficient {
    /// Divide by two. Impls must be exact: `x.half() + x.half() == x`.
    fn half(&self) -> Self;
}

impl Halvable for f64 {
    #[inline]
    fn half(&self) -> Self {
        *self / 2.0
    }
}

impl Halvable for num::Complex<f64> {
    #[inline]
    fn half(&self) -> Self {
        *self / 2.0
    }
}
