use std::marker::PhantomData;

use crate::{config::Config, word::PauliWord};
use crate::traits::Coefficient;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Byte<const N: usize, C: Coefficient + Sync + Send>(PhantomData<C>);

impl<const N: usize, C: Coefficient + Sync + Send> Config for Byte<N, C> {
    type Storage = [u8; N];
    type Coeff = C;
    type Map = dashmap::DashMap<PauliWord<[u8; N]>, C>;
}

pub type ByteF64<const N: usize> = Byte<N, f64>;
