use std::collections::HashMap;
use std::marker::PhantomData;

use crate::traits::{Coefficient, NoStrategy, PauliWordTrait, Strategy};
use crate::{config::Config, word::PauliWord};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Byte<
    const N: usize,
    C: Coefficient,
    St: Strategy = NoStrategy,
    W: PauliWordTrait<[u8; N], fxhash::FxBuildHasher> = PauliWord<[u8; N], fxhash::FxBuildHasher>,
>(PhantomData<(C, St, W)>);

impl<
    const N: usize,
    C: Coefficient,
    St: Strategy,
    W: PauliWordTrait<[u8; N], fxhash::FxBuildHasher>,
> Config for Byte<N, C, St, W>
{
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = fxhash::FxBuildHasher;
    type PauliWordType = W;
    type Map = HashMap<W, C, Self::BuildHasher>;
    type Strategy = St;
}

pub type ByteF64<const N: usize, St = NoStrategy> = Byte<N, f64, St>;
