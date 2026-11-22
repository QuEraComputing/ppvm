// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Operator boilerplate shared by the arithmetic modules.

/// `impl Op<R> for Term` in terms of `OpAssign<R>`, for each listed `R`.
macro_rules! binary_ops {
    ($op:ident::$f:ident from $assign:ident::$af:ident for $($rhs:ty),*) => {$(
        impl std::ops::$op<$rhs> for crate::Term {
            type Output = crate::Term;

            #[inline]
            fn $f(mut self, rhs: $rhs) -> crate::Term {
                std::ops::$assign::$af(&mut self, rhs);
                self
            }
        }
    )*};
}
