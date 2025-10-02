use std::collections::HashMap;
use std::marker::PhantomData;

use crate::traits::Coefficient;
use crate::{config::Config, word::PauliWord};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Byte<const N: usize, C: Coefficient>(PhantomData<C>);

impl<const N: usize, C: Coefficient> Config for Byte<N, C> {
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = fxhash::FxBuildHasher;
    type Map = HashMap<PauliWord<[u8; N], Self::BuildHasher>, C, Self::BuildHasher>;
}

pub type ByteF64<const N: usize> = Byte<N, f64>;
