// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::collections::HashMap;
use std::marker::PhantomData;

use ppvm_pauli_word::word::PauliWord;
use ppvm_traits::config::Config;
use ppvm_traits::traits::{Coefficient, NoStrategy, PauliWordTrait, Strategy};

/// `HashMap`-backed [`Config`] with native-word (`[usize; N]`) storage and
/// `FxHasher`.
///
/// The storage element is `usize`, i.e. the target's native machine word:
/// `u64` on 64-bit targets (identical layout and performance to a hardcoded
/// `[u64; N]`) and `u32` on 32-bit targets such as `wasm32`. Using `usize`
/// rather than `u64` keeps this config available on every target — `bitvec`
/// only implements `BitStore` for `u64` on 64-bit pointer widths, but always
/// implements it for `usize`. (The `Byte8` name refers to the 64-bit word
/// this config packs into on native targets.)
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Byte8<
    const N: usize,
    C: Coefficient,
    St: Strategy = NoStrategy,
    W: PauliWordTrait = PauliWord<[usize; N], fxhash::FxBuildHasher>,
>(PhantomData<(C, St, W)>);

impl<const N: usize, C: Coefficient, St: Strategy, W: PauliWordTrait> Config
    for Byte8<N, C, St, W>
{
    type Storage = [usize; N];
    type Coeff = C;
    type BuildHasher = fxhash::FxBuildHasher;
    type PauliWordType = W;
    type Map = HashMap<W, C, Self::BuildHasher>;
    type Strategy = St;
}

/// [`Byte8`] specialised to `f64` coefficients.
pub type Byte8F64<const N: usize, St = NoStrategy> = Byte8<N, f64, St>;
