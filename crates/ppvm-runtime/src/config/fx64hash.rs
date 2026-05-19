// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::collections::HashMap;
use std::marker::PhantomData;

use crate::traits::{Coefficient, NoStrategy, PauliWordTrait, Strategy};
use crate::{config::Config, word::PauliWord};

/// `HashMap`-backed [`Config`] with `[u64; N]` storage and `FxHasher`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Byte8<
    const N: usize,
    C: Coefficient,
    St: Strategy = NoStrategy,
    W: PauliWordTrait = PauliWord<[u64; N], fxhash::FxBuildHasher>,
>(PhantomData<(C, St, W)>);

impl<const N: usize, C: Coefficient, St: Strategy, W: PauliWordTrait> Config
    for Byte8<N, C, St, W>
{
    type Storage = [u64; N];
    type Coeff = C;
    type BuildHasher = fxhash::FxBuildHasher;
    type PauliWordType = W;
    type Map = HashMap<W, C, Self::BuildHasher>;
    type Strategy = St;
}

/// [`Byte8`] specialised to `f64` coefficients.
pub type Byte8F64<const N: usize, St = NoStrategy> = Byte8<N, f64, St>;
