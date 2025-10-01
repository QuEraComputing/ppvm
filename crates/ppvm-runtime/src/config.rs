use std::marker::PhantomData;

use crate::traits::{Coefficient, Map, PauliStorage};

pub trait Config {
    type Storage: PauliStorage;
    type Value: Coefficient;
    type MapType: Map<Self::Storage, Self::Value>;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Byte<const N: usize, C: Coefficient, M: Map<[u8; N], C>>(PhantomData<(C, M)>);

impl<const N: usize, C: Coefficient, M: Map<[u8; N], C>> Config for Byte<N, C, M> {
    type Storage = [u8; N];
    type Value = C;
    type MapType = M;
}

pub type ByteF64<const N: usize, M> = Byte<N, f64, M>;
