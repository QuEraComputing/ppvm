use std::marker::PhantomData;

use crate::traits::{Coefficient, NoStrategy, PauliWordTrait, Strategy};
use crate::{config::Config, word::PauliWord};

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

pub type ByteFxHashF64<const N: usize, St = NoStrategy> = ByteFxHash<N, f64, St>;
pub type ByteGxHashF64<const N: usize, St = NoStrategy> = ByteGxHash<N, f64, St>;
