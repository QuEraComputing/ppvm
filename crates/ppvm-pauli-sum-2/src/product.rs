// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Distributive products: `(Σ a_p P)(Σ b_q Q) = Σ a_p b_q i^k (PQ)`, where the word
//! product's residual phase `i^k` is absorbed into the coefficient. This needs a
//! coefficient ring with an imaginary unit, so real-coefficient sums and lossy
//! words have no product. Exact zeros stay stored; reduction and truncation are
//! explicit. Results keep the left operand's policy and preserved keys.

use std::ops::{Mul, MulAssign};

use ppvm_traits_2::{ImaginaryUnit, KeyProduct, Multiply};

use crate::{StoreAlloc, Sum};

/// Multiplies every pair of terms, so this costs `|self| · |rhs|` word products.
/// Order matters, since Pauli words need not commute. Policies may differ; the
/// result keeps the left one. Panics if the widths differ.
impl<S, P, P2> Mul<&Sum<S, P2>> for &Sum<S, P>
where
    S: Multiply + StoreAlloc,
    S::Key: KeyProduct,
    S::Coeff: ImaginaryUnit,
    P: Clone,
{
    type Output = Sum<S, P>;

    fn mul(self, rhs: &Sum<S, P2>) -> Self::Output {
        assert_eq!(
            self.n_sites(),
            rhs.n_sites(),
            "product requires equal widths"
        );
        let mut out = self.empty_like(self.len().max(rhs.len()));
        self.storage.multiply_into(&rhs.storage, &mut out.storage);
        out
    }
}

impl<S, P, P2> MulAssign<&Sum<S, P2>> for Sum<S, P>
where
    S: Multiply + StoreAlloc,
    S::Key: KeyProduct,
    S::Coeff: ImaginaryUnit,
    P: Clone,
{
    fn mul_assign(&mut self, rhs: &Sum<S, P2>) {
        *self = &*self * rhs;
    }
}
