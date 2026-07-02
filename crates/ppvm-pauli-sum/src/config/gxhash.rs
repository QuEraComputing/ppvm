// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::collections::HashMap;
use std::marker::PhantomData;

use ppvm_pauli_word::word::PauliWord;
use ppvm_traits::config::Config;
use ppvm_traits::traits::{Coefficient, NoStrategy, PauliWordTrait, Strategy};

/// `HashMap`-backed [`Config`] with `[u8; N]` storage and `gxhash`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Byte<
    const N: usize,
    C: Coefficient,
    St: Strategy = NoStrategy,
    W: PauliWordTrait = PauliWord<[u8; N], gxhash::GxBuildHasher>,
>(PhantomData<(C, St, W)>);

impl<const N: usize, C: Coefficient, St: Strategy, W: PauliWordTrait + Send + Sync> Config
    for Byte<N, C, St, W>
{
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = gxhash::GxBuildHasher;
    type PauliWordType = W;
    type Map = HashMap<W, C, Self::BuildHasher>;
    type Strategy = St;
}

/// [`Byte`] specialised to `f64` coefficients.
pub type ByteF64<const N: usize, St = NoStrategy, W = PauliWord<[u8; N], gxhash::GxBuildHasher>> =
    Byte<N, f64, St, W>;
