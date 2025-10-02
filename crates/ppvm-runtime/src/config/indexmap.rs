use std::marker::PhantomData;

use crate::{config::Config, word::PauliWord};
use crate::traits::Coefficient;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ByteFxHash<const N: usize, C: Coefficient>(PhantomData<C>);

impl<const N: usize, C: Coefficient> Config for ByteFxHash<N, C> {
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = fxhash::FxBuildHasher;
    type Map = indexmap::IndexMap<PauliWord<[u8; N], Self::BuildHasher>, Self::Coeff, Self::BuildHasher>;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ByteGxHash<const N: usize, C: Coefficient>(PhantomData<C>);

impl<const N: usize, C: Coefficient> Config for ByteGxHash<N, C> {
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = gxhash::GxBuildHasher;
    type Map = indexmap::IndexMap<PauliWord<[u8; N], Self::BuildHasher>, Self::Coeff, Self::BuildHasher>;
}

pub type ByteFxHashF64<const N: usize> = ByteFxHash<N, f64>;
pub type ByteGxHashF64<const N: usize> = ByteGxHash<N, f64>;
