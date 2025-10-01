use std::marker::PhantomData;

use crate::traits::{Coefficient, ACMap, PauliStorage};

pub trait Config {
    type Storage: PauliStorage;
    type Coeff: Coefficient;
    type Map: ACMap<Self::Storage, Self::Coeff>;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Byte<const N: usize, C: Coefficient, M: ACMap<[u8; N], C>>(PhantomData<(C, M)>);

impl<const N: usize, C: Coefficient, M: ACMap<[u8; N], C>> Config for Byte<N, C, M> {
    type Storage = [u8; N];
    type Coeff = C;
    type Map = M;
}

pub type ByteF64<const N: usize, M> = Byte<N, f64, M>;
