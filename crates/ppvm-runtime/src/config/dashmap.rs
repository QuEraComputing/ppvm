use std::marker::PhantomData;

use crate::traits::{Coefficient, NoStrategy, PauliWordTrait, Strategy};
use crate::{config::Config, word::PauliWord};

/// `DashMap`-backed concurrent [`Config`] with `[u8; N]` storage and `FxHasher`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ByteFxHash<
    const N: usize,
    C: Coefficient + Sync + Send,
    St: Strategy = NoStrategy,
    W: PauliWordTrait = PauliWord<[u8; N], fxhash::FxBuildHasher>,
>(PhantomData<(C, St, W)>);

impl<const N: usize, C: Coefficient + Sync + Send, St: Strategy, W: PauliWordTrait + Sync + Send>
    Config for ByteFxHash<N, C, St, W>
{
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = fxhash::FxBuildHasher;
    type PauliWordType = W;
    type Map = dashmap::DashMap<W, C, Self::BuildHasher>;
    type Strategy = St;
}

/// `DashMap`-backed concurrent [`Config`] with `[u8; N]` storage and `gxhash`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ByteGxHash<
    const N: usize,
    C: Coefficient + Sync + Send,
    St: Strategy = NoStrategy,
    W: PauliWordTrait + Sync + Send = PauliWord<[u8; N], gxhash::GxBuildHasher>,
>(PhantomData<(C, St, W)>);

impl<const N: usize, C: Coefficient + Sync + Send, St: Strategy, W: PauliWordTrait + Sync + Send>
    Config for ByteGxHash<N, C, St, W>
{
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = gxhash::GxBuildHasher;
    type PauliWordType = W;
    type Map = dashmap::DashMap<W, C, Self::BuildHasher>;
    type Strategy = St;
}

/// [`ByteFxHash`] specialised to `f64` coefficients.
pub type ByteFxHashF64<const N: usize, St = NoStrategy> = ByteFxHash<N, f64, St>;
/// [`ByteGxHash`] specialised to `f64` coefficients.
pub type ByteGxHashF64<const N: usize, St = NoStrategy> = ByteGxHash<N, f64, St>;
