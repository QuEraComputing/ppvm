use std::marker::PhantomData;

use crate::traits::Coefficient;
use crate::{config::Config, word::PauliWord};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ByteFxHash<const N: usize, C: Coefficient + Sync + Send>(PhantomData<C>);

impl<const N: usize, C: Coefficient + Sync + Send> Config for ByteFxHash<N, C> {
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = fxhash::FxBuildHasher;
    type Map = dashmap::DashMap<PauliWord<[u8; N], Self::BuildHasher>, C, Self::BuildHasher>;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ByteGxHash<const N: usize, C: Coefficient + Sync + Send>(PhantomData<C>);

impl<const N: usize, C: Coefficient + Sync + Send> Config for ByteGxHash<N, C> {
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = gxhash::GxBuildHasher;
    type Map = dashmap::DashMap<PauliWord<[u8; N], Self::BuildHasher>, C, Self::BuildHasher>;
}

pub type ByteFxHashF64<const N: usize> = ByteFxHash<N, f64>;
pub type ByteGxHashF64<const N: usize> = ByteGxHash<N, f64>;
