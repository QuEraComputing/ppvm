use std::marker::PhantomData;

use crate::traits::{Coefficient, NoStrategy, Strategy};
use crate::{config::Config, word::PauliWord};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ByteFxHash<const N: usize, C: Coefficient, St: Strategy = NoStrategy>(
    PhantomData<(C, St)>,
);

impl<const N: usize, C: Coefficient, St: Strategy> Config for ByteFxHash<N, C, St> {
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = fxhash::FxBuildHasher;
    type Map =
        indexmap::IndexMap<PauliWord<[u8; N], Self::BuildHasher>, Self::Coeff, Self::BuildHasher>;
    type Strategy = St;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ByteGxHash<const N: usize, C: Coefficient, St: Strategy = NoStrategy>(
    PhantomData<(C, St)>,
);

impl<const N: usize, C: Coefficient, St: Strategy> Config for ByteGxHash<N, C, St> {
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = gxhash::GxBuildHasher;
    type Map =
        indexmap::IndexMap<PauliWord<[u8; N], Self::BuildHasher>, Self::Coeff, Self::BuildHasher>;
    type Strategy = St;
}

pub type ByteFxHashF64<const N: usize, St = NoStrategy> = ByteFxHash<N, f64, St>;
pub type ByteGxHashF64<const N: usize, St = NoStrategy> = ByteGxHash<N, f64, St>;
