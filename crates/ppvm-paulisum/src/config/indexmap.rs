// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::marker::PhantomData;

use ppvm_runtime::traits::{Coefficient, NoStrategy, PauliWordTrait, Strategy};
use ppvm_runtime::config::Config;
use ppvm_runtime::word::PauliWord;

/// `IndexMap`-backed [`Config`] with `[u8; N]` storage and `FxHasher`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ByteFxHash<
    const N: usize,
    C: Coefficient,
    St: Strategy = NoStrategy,
    W: PauliWordTrait = PauliWord<[u8; N], fxhash::FxBuildHasher>,
>(PhantomData<(C, St, W)>);

impl<const N: usize, C: Coefficient, St: Strategy, W: PauliWordTrait> Config
    for ByteFxHash<N, C, St, W>
{
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = fxhash::FxBuildHasher;
    type PauliWordType = W;
    type Map = indexmap::IndexMap<Self::PauliWordType, Self::Coeff, Self::BuildHasher>;
    type Strategy = St;
}

/// `IndexMap`-backed [`Config`] with `[u8; N]` storage and `gxhash`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ByteGxHash<
    const N: usize,
    C: Coefficient,
    St: Strategy = NoStrategy,
    W: PauliWordTrait = PauliWord<[u8; N], gxhash::GxBuildHasher>,
>(PhantomData<(C, St, W)>);

impl<const N: usize, C: Coefficient, St: Strategy, W: PauliWordTrait> Config
    for ByteGxHash<N, C, St, W>
{
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = gxhash::GxBuildHasher;
    type PauliWordType = W;
    type Map = indexmap::IndexMap<W, Self::Coeff, Self::BuildHasher>;
    type Strategy = St;
}

/// [`ByteFxHash`] specialised to `f64` coefficients.
pub type ByteFxHashF64<
    const N: usize,
    St = NoStrategy,
    Wd = PauliWord<[u8; N], fxhash::FxBuildHasher>,
> = ByteFxHash<N, f64, St, Wd>;
/// [`ByteGxHash`] specialised to `f64` coefficients.
pub type ByteGxHashF64<
    const N: usize,
    St = NoStrategy,
    Wd = PauliWord<[u8; N], gxhash::GxBuildHasher>,
> = ByteGxHash<N, f64, St, Wd>;
