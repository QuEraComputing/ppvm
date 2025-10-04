use std::collections::HashMap;
use std::marker::PhantomData;

use crate::traits::{Coefficient, NoStrategy, Strategy};
use crate::{config::Config, word::PauliWord};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Byte<const N: usize, C: Coefficient, St: Strategy = NoStrategy>(PhantomData<(C, St)>);

impl<const N: usize, C: Coefficient, St: Strategy> Config for Byte<N, C, St> {
    type Storage = [u8; N];
    type Coeff = C;
    type BuildHasher = gxhash::GxBuildHasher;
    type Map = HashMap<PauliWord<[u8; N], Self::BuildHasher>, C, Self::BuildHasher>;
    type Strategy = St;
}

pub type ByteF64<const N: usize, St = NoStrategy> = Byte<N, f64, St>;
