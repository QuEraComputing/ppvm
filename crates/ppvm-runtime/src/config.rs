use std::marker::PhantomData;

use crate::traits::{ACMap, Coefficient, PauliStorage};
#[cfg(feature = "dashmap")]
use crate::word::PauliWord;

pub trait Config {
    type Storage: PauliStorage;
    type Coeff: Coefficient;
    type Map: ACMap<Self::Storage, Self::Coeff>;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Byte<const N: usize, C: Coefficient, M: ACMap<[u8; N], C>>(PhantomData<(C, M)>);
pub type ByteF64<const N: usize, M> = Byte<N, f64, M>;

impl<const N: usize, C: Coefficient, M: ACMap<[u8; N], C>> Config for Byte<N, C, M> {
    type Storage = [u8; N];
    type Coeff = C;
    type Map = M;
}

#[cfg(feature = "dashmap")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ByteDashMap<const N: usize, C: Coefficient>(PhantomData<C>);

#[cfg(feature = "dashmap")]
impl<const N: usize, C: Coefficient> Config for ByteDashMap<N, C> {
    type Storage = [u8; N];
    type Coeff = C;
    type Map = dashmap::DashMap<PauliWord<[u8; N]>, C>;
}

pub type ByteDashMapF64<const N: usize> = ByteDashMap<N, f64>;
